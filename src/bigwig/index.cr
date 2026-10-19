module BigWig::Internal
  # :nodoc:
  # Only the current path is resident. A whole-chromosome query does not first
  # collect every matching block or expand the complete R-tree.
  class BlockIterator
    include Iterator(BlockRef)

    private class Frame
      getter offset : Int64
      getter? leaf : Bool
      getter count : Int32
      property item : Int32 = 0
      property cursor : Cursor

      def initialize(@offset, @leaf, @count, @cursor)
      end
    end

    @stack = Array(Frame).new
    @started = false
    @finished = false
    @fanout : Int32 = 0

    def initialize(@source : Source, @index_offset : Int64, @region : Region,
                   @big_endian : Bool, @limit : Int32)
    end

    private def load_node(offset : Int64) : Frame
      if @stack.size >= MAX_DEPTH || @stack.any? { |node| node.offset == offset }
        raise BigWig::FormatError.new("cyclic or excessively deep R-tree")
      end
      c = Cursor.new(@source.read(offset, 4, @limit), @big_endian)
      flag = c.u8
      raise BigWig::FormatError.new("invalid R-tree node flag") unless flag <= 1
      c.skip(1)
      count = c.u16.to_i
      raise BigWig::FormatError.new("R-tree node exceeds fanout") if count > @fanout
      width = flag == 1 ? 32 : 24
      data = @source.read(offset + 4, count * width, @limit)
      Frame.new(offset, flag == 1, count, Cursor.new(data, @big_endian))
    end

    private def start_iteration : Nil
      @started = true
      if @region.start == @region.stop || @index_offset == 0
        @finished = true
        return
      end
      cursor = Cursor.new(@source.read(@index_offset, 48, @limit), @big_endian)
      raise BigWig::FormatError.new("invalid R-tree signature") unless cursor.u32 == RTREE_MAGIC
      fanout = cursor.u32
      raise BigWig::FormatError.new("invalid R-tree fanout") unless 1 <= fanout <= 65_535
      @fanout = fanout.to_i
      cursor.skip(40)
      @stack << load_node(@index_offset + 48)
    end

    private def read_ref(frame : Frame) : BlockRef
      cursor = frame.cursor
      start_chrom = cursor.u32
      first = cursor.u32.to_i64
      end_chrom = cursor.u32
      last = cursor.u32.to_i64
      offset = cursor.offset
      size = frame.leaf? ? cursor.offset : 0_i64
      frame.cursor = cursor
      frame.item += 1
      if start_chrom > end_chrom || (start_chrom == end_chrom && first > last)
        raise BigWig::FormatError.new("invalid R-tree bounding interval")
      end
      BlockRef.new(start_chrom, first, end_chrom, last, offset, size)
    end

    def next : BlockRef | Iterator::Stop
      raise BigWig::ClosedError.new("reader is closed") if @source.closed?
      return stop if @finished
      start_iteration unless @started
      until @stack.empty?
        frame = @stack.last
        if frame.item == frame.count
          @stack.pop
          next
        end
        ref = read_ref(frame)
        next unless ref.overlaps?(@region)
        if frame.leaf?
          raise BigWig::FormatError.new("empty indexed data block") unless ref.size > 0
          @source.check_range(ref.offset, ref.size)
          return ref
        else
          @stack << load_node(ref.offset)
        end
      end
      @finished = true
      stop
    end
  end

  # :nodoc:
  class BlockLoader
    @compressed = Bytes.empty
    @decoded = Bytes.empty

    def initialize(@source : Source, @buffer_size : Int32, @limit : Int32)
    end

    def load(ref : BlockRef) : Bytes
      raise BigWig::LimitError.new("compressed block exceeds max_block_bytes") if ref.size > @limit
      count = ref.size.to_i
      @compressed = Bytes.new(count) if @compressed.size < count
      data = @compressed[0, count]
      @source.read_exact_at(ref.offset, data)
      return data if @buffer_size == 0
      @decoded = Bytes.new(@buffer_size) if @decoded.size < @buffer_size
      stream = IO::Memory.new(data)
      total = 0
      Compress::Zlib::Reader.open(stream) do |zlib|
        while total < @buffer_size
          n = zlib.read(@decoded[total, @buffer_size - total])
          break if n == 0
          total += n
        end
        # Reading through EOF both checks the Adler checksum and catches
        # decompression exceeding the declared bound, without growing memory.
        unless zlib.read_byte.nil?
          raise BigWig::FormatError.new("decompressed block exceeds header buffer size")
        end
      end
      @decoded[0, total]
    rescue ex : Compress::Zlib::Error
      raise BigWig::FormatError.new("invalid zlib block: #{ex.message}")
    rescue ex : Compress::Deflate::Error
      raise BigWig::FormatError.new("invalid deflate block: #{ex.message}")
    rescue ex : IO::EOFError
      raise BigWig::FormatError.new("truncated compressed block: #{ex.message}")
    end
  end
end

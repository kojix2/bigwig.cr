module BigWig
  # Reads local/remote bigWig and bigBed files. Prefer BigWig.open with a block.
  # Each iterator has its own traversal and decode buffers. Advancing one
  # iterator cannot invalidate values returned by another.
  class Reader
    @file_type = FileType::BigWig
    @big_endian = false
    @buffer_size = 0
    @index_offset = 0_i64
    @sql_offset = 0_i64
    @sql_cache : String? = nil
    @sql_loaded = false
    @zoom_levels = Array(Internal::ZoomLevel).new
    @chroms = Hash(String, Internal::Chromosome).new
    @chrom_sizes = Hash(UInt32, Int64).new
    @header = Header.new(0, 0, 0_u64, Float64::NAN, Float64::NAN,
      Float64::NAN, Float64::NAN, 0, 0, false)

    # :nodoc:
    def initialize(@source : Internal::Source, @options : ReadOptions = ReadOptions.new)
      @options.validate
      begin
        read_header
      rescue ex
        @source.close
        raise ex
      end
    end

    def self.new(path : String | Path, *, options : ReadOptions = ReadOptions.new,
                 headers : HTTP::Headers = HTTP::Headers.new)
      BigWig.open(path, options: options, headers: headers)
    end

    def self.new(io : IO, *, sync_close : Bool = false, options : ReadOptions = ReadOptions.new)
      BigWig.open(io, sync_close: sync_close, options: options)
    end

    def close : Nil
      @source.close
    end

    def closed? : Bool
      @source.closed?
    end

    private def check_open : Nil
      raise ClosedError.new("reader is closed") if closed?
    end

    def file_type : FileType
      check_open
      @file_type
    end

    def bigwig? : Bool
      file_type == FileType::BigWig
    end

    def bigbed? : Bool
      file_type == FileType::BigBed
    end

    def header : Header
      check_open
      @header
    end

    # A copy in chromosome-ID order, not a mutable view of internal metadata.
    def chromosomes : Hash(String, Int64)
      check_open
      result = Hash(String, Int64).new
      @chroms.values.sort_by!(&.id).each { |chrom| result[chrom.name] = chrom.size }
      result
    end

    def chromosome_size(chrom : String) : Int64?
      check_open
      @chroms[chrom]?.try(&.size)
    end

    # Reduction factors actually stored in the file, from fine to coarse.
    def zoom_levels : Array(Int64)
      check_open
      @zoom_levels.map(&.reduction)
    end

    # Unparsed AutoSQL, or nil when absent. Only available for bigBed.
    def sql : String?
      require_bigbed
      return @sql_cache if @sql_loaded
      if @sql_offset != 0
        data = IO::Memory.new
        offset = @sql_offset
        loop do
          size = Math.min(1024_i64, @source.size - offset).to_i
          raise FormatError.new("unterminated AutoSQL") if size <= 0
          bytes = @source.read(offset, size, @options.max_block_bytes)
          length = bytes.index(0_u8)
          count = length || size
          if data.size.to_i64 + count > @options.max_block_bytes
            raise LimitError.new("AutoSQL exceeds max_block_bytes")
          end
          data.write(bytes[0, count])
          if length
            @sql_cache = data.to_s
            break
          end
          offset += size
        end
      end
      @sql_loaded = true
      @sql_cache
    end

    private def require_bigwig : Nil
      check_open
      raise UnsupportedError.new("this operation requires bigWig") unless @file_type == FileType::BigWig
    end

    private def require_bigbed : Nil
      check_open
      raise UnsupportedError.new("this operation requires bigBed") unless @file_type == FileType::BigBed
    end

    private def region(chrom : String, start : Int, stop : Int?) : Internal::Region
      check_open
      chromosome = @chroms[chrom]? || raise ArgumentError.new("unknown chromosome: #{chrom}")
      first = Internal.coordinate(start, "start")
      last = stop ? Internal.coordinate(stop, "stop") : chromosome.size
      unless first <= last <= chromosome.size
        raise ArgumentError.new("expected 0 <= start <= stop <= #{chromosome.size} for #{chrom}")
      end
      Internal::Region.new(chromosome.id, first, last)
    end

    private def blocks(query : Internal::Region, offset : Int64 = @index_offset)
      Internal::BlockIterator.new(@source, offset, query, @big_endian, @options.max_block_bytes)
    end

    private def loader
      Internal::BlockLoader.new(@source, @buffer_size, @options.max_block_bytes)
    end

    private def interval_iterator(query : Internal::Region)
      Internal::IntervalIterator.new(@source, blocks(query), loader, query, @big_endian, @chrom_sizes)
    end

    # Lazily yields whole overlapping records. Empty ranges yield nothing.
    def each_interval(chrom : String, start : Int = 0, stop : Int? = nil)
      require_bigwig
      interval_iterator(region(chrom, start, stop))
    end

    def each_interval(chrom : String, start : Int = 0, stop : Int? = nil, & : Interval ->) : Nil
      each_interval(chrom, start, stop).each { |value| yield value }
    end

    def each_entry(chrom : String, start : Int = 0, stop : Int? = nil)
      require_bigbed
      query = region(chrom, start, stop)
      Internal::EntryIterator.new(@source, blocks(query), loader, query, @big_endian, @chrom_sizes)
    end

    def each_entry(chrom : String, start : Int = 0, stop : Int? = nil, & : BedEntry ->) : Nil
      each_entry(chrom, start, stop).each { |value| yield value }
    end

    def intervals(chrom : String, start : Int = 0, stop : Int? = nil) : Array(Interval)
      result = Array(Interval).new
      each_interval(chrom, start, stop) do |value|
        check_allocation((result.size.to_i64 + 1) * sizeof(Interval))
        result << value
      end
      result
    end

    def entries(chrom : String, start : Int = 0, stop : Int? = nil) : Array(BedEntry)
      result = Array(BedEntry).new
      bytes = 0_i64
      each_entry(chrom, start, stop) do |value|
        bytes += sizeof(BedEntry) + value.rest.bytesize
        check_allocation(bytes)
        result << value
      end
      result
    end

    private def check_allocation(bytes : Int64) : Nil
      if bytes > @options.max_array_bytes
        raise LimitError.new("result exceeds max_array_bytes; use an iterator or a reusable buffer")
      end
    end

    # Dense values, with NaN for uncovered bases. Allocation is proportional to
    # stop-start; use read_values/each_interval for unbounded regions.
    def values(chrom : String, start : Int, stop : Int) : Array(Float32)
      require_bigwig
      query = region(chrom, start, stop)
      count = query.stop - query.start
      check_allocation(count * sizeof(Float32))
      raise LimitError.new("dense result exceeds Crystal's array size limit") if count > Int32::MAX
      result = Array(Float32).new(count.to_i, Float32::NAN)
      fill_values(query, result.to_unsafe.to_slice(result.size))
      result
    end

    # Fills every element, including NaN gaps. Returns buffer.size; never a
    # short read. A range extending past the chromosome raises ArgumentError.
    def read_values(chrom : String, start : Int, buffer : Slice(Float32)) : Int32
      require_bigwig
      first = Internal.coordinate(start, "start")
      query = region(chrom, first, first + buffer.size)
      fill_values(query, buffer)
      buffer.size
    end

    private def fill_values(query : Internal::Region, buffer : Slice(Float32)) : Nil
      buffer.fill(Float32::NAN)
      interval_iterator(query).each do |item|
        first = Math.max(item.start, query.start)
        last = Math.min(item.stop, query.stop)
        buffer[(first - query.start).to_i, (last - first).to_i].fill(item.value)
      end
    end

    def value(chrom : String, position : Int) : Float32
      buffer = StaticArray(Float32, 1).new(Float32::NAN)
      read_values(chrom, position, buffer.to_slice)
      buffer[0]
    end

    def each_value(chrom : String, start : Int, stop : Int)
      require_bigwig
      query = region(chrom, start, stop)
      Internal::ValueIterator.new(@source, interval_iterator(query), query.start, query.stop)
    end

    def each_value(chrom : String, start : Int, stop : Int, & : Float32 ->) : Nil
      each_value(chrom, start, stop).each { |value| yield value }
    end

    {% for method in %w[each_interval each_entry each_value] %}
      def {{ method.id }}(chrom : String, range : Range(B, E)) forall B, E
        first, last = Internal.range_bounds(range)
        {{ method.id }}(chrom, first, last)
      end

      def {{ method.id }}(chrom : String, range : Range(B, E), &) forall B, E
        first, last = Internal.range_bounds(range)
        {{ method.id }}(chrom, first, last) { |value| yield value }
      end
    {% end %}

    {% for method in %w[intervals entries values] %}
      def {{ method.id }}(chrom : String, range : Range(B, E)) forall B, E
        first, last = Internal.range_bounds(range)
        {{ method.id }}(chrom, first, last)
      end
    {% end %}

    private def read_file_magic(bytes : Bytes) : Nil
      magic = IO::ByteFormat::LittleEndian.decode(UInt32, bytes[0, 4])
      unless {Internal::BW_MAGIC, Internal::BB_MAGIC}.includes?(magic)
        magic = IO::ByteFormat::BigEndian.decode(UInt32, bytes[0, 4])
        unless {Internal::BW_MAGIC, Internal::BB_MAGIC}.includes?(magic)
          raise FormatError.new("not a bigWig or bigBed file")
        end
        @big_endian = true
      end
      @file_type = magic == Internal::BW_MAGIC ? FileType::BigWig : FileType::BigBed
    end

    private def read_zoom_levels(levels : Int32) : Nil
      return if levels == 0
      cursor = Internal::Cursor.new(@source.read(64_i64, levels * 24, @options.max_block_bytes), @big_endian)
      levels.times do
        reduction = cursor.u32.to_i64
        cursor.skip(4)
        data = cursor.offset
        index = cursor.offset
        raise FormatError.new("invalid zoom reduction") if reduction == 0
        @source.check_range(data, 4)
        @source.check_range(index, 48)
        @zoom_levels << Internal::ZoomLevel.new(reduction, data, index)
      end
      @zoom_levels.sort_by!(&.reduction)
    end

    private def read_summary(offset : Int64)
      return {0_u64, Float64::NAN, Float64::NAN, Float64::NAN, Float64::NAN} if offset == 0
      cursor = Internal::Cursor.new(@source.read(offset, 40, @options.max_block_bytes), @big_endian)
      {cursor.u64, cursor.f64, cursor.f64, cursor.f64, cursor.f64}
    end

    private def validate_index(data_offset : Int64, version : Int32) : Nil
      if @index_offset != 0
        @source.check_range(@index_offset, 48)
      elsif data_offset != 0 && version >= 2
        count = Internal::Cursor.new(@source.read(data_offset, 8, @options.max_block_bytes), @big_endian).u64
        raise FormatError.new("nonempty file has no data index") unless count == 0
      end
    end

    private def read_header : Nil
      bytes = @source.read(0_i64, 64, @options.max_block_bytes)
      read_file_magic(bytes)
      c = Internal::Cursor.new(bytes, @big_endian)
      c.skip(4)
      version = c.u16.to_i
      raise UnsupportedError.new("unsupported BBI version #{version}") unless 1 <= version <= 4
      levels = c.u16.to_i
      chrom_offset = c.offset
      data_offset = c.offset
      @index_offset = c.offset
      fields = c.u16.to_i
      defined_fields = c.u16.to_i
      @sql_offset = c.offset
      summary_offset = c.offset
      buffer_size = c.u32
      if buffer_size > @options.max_block_bytes
        raise LimitError.new("header decompression buffer exceeds max_block_bytes")
      end
      @buffer_size = buffer_size.to_i
      extension_offset = c.offset
      unless chrom_offset >= 64 && (data_offset == 0 || data_offset >= 64) && (@index_offset == 0 || @index_offset >= 64)
        raise FormatError.new("invalid header offsets")
      end
      read_zoom_levels(levels)
      bases, minimum, maximum, sum, sum_squared = read_summary(summary_offset)
      @header = Header.new(version, levels, bases, minimum, maximum, sum, sum_squared,
        fields, defined_fields, summary_offset != 0)
      @source.check_range(extension_offset, 2) if extension_offset != 0
      @source.check_range(@sql_offset, 1) if @sql_offset != 0
      validate_index(data_offset, version)
      read_chromosomes(chrom_offset)
    end

    private def read_chromosome_tree_header(offset : Int64) : Tuple(UInt32, UInt32, UInt64, Int64)
      header = Internal::Cursor.new(@source.read(offset, 32, @options.max_block_bytes), @big_endian)
      raise FormatError.new("invalid chromosome B+tree signature") unless header.u32 == Internal::BPT_MAGIC
      fanout = header.u32
      key_size = header.u32
      value_size = header.u32
      count = header.u64
      unless fanout > 0 && key_size > 0 && key_size <= @options.max_block_bytes && value_size == 8
        raise FormatError.new("invalid chromosome B+tree header")
      end
      {fanout, key_size, count, offset + 32}
    end

    private def read_chromosome_node(offset : Int64, depth : Int32, key_size : UInt32,
                                     fanout : UInt32) : Array(Int64)
      if depth >= Internal::MAX_DEPTH
        raise FormatError.new("cyclic or excessively deep chromosome B+tree")
      end
      node = Internal::Cursor.new(@source.read(offset, 4, @options.max_block_bytes), @big_endian)
      flag = node.u8
      raise FormatError.new("invalid chromosome B+tree node flag") unless flag <= 1
      node.skip(1)
      items = node.u16.to_i
      raise FormatError.new("chromosome B+tree node exceeds fanout") if items > fanout
      size = items.to_i64 * (key_size.to_i64 + 8)
      data = Internal::Cursor.new(@source.read(offset + 4, size, @options.max_block_bytes), @big_endian)
      children = Array(Int64).new
      items.times do
        key = data.take(key_size)
        if flag == 1
          name = String.new(key[0, key.index(0_u8) || key.size])
          id = data.u32
          length = data.u32.to_i64
          if name.empty? || @chroms.has_key?(name) || @chrom_sizes.has_key?(id)
            raise FormatError.new("empty or duplicate chromosome name/ID")
          end
          @chroms[name] = Internal::Chromosome.new(name, id, length)
          @chrom_sizes[id] = length
        else
          children << data.offset
        end
      end
      children
    end

    private def read_chromosomes(offset : Int64) : Nil
      fanout, key_size, count, root = read_chromosome_tree_header(offset)
      stack = [{root, 0}]
      visited = Set(Int64).new
      until stack.empty?
        node_offset, depth = stack.pop
        if !visited.add?(node_offset)
          raise FormatError.new("cyclic or excessively deep chromosome B+tree")
        end
        children = read_chromosome_node(node_offset, depth, key_size, fanout)
        if @chroms.size.to_u64 > count
          raise FormatError.new("chromosome count exceeds header count")
        end
        children.reverse_each { |child| stack << {child, depth + 1} }
      end
      raise FormatError.new("chromosome count does not match header") unless @chroms.size.to_u64 == count
    end
  end
end

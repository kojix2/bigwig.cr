module BigWig::Internal
  def self.write_data_block(io : IO, data : Bytes, options : BigWig::WriteOptions) : Tuple(Int64, Int64)
    offset = io.pos.to_i64
    if options.compression
      Compress::Zlib::Writer.open(io, level: options.compression_level) { |zlib_writer| zlib_writer.write(data) }
    else
      io.write(data)
    end
    {offset, io.pos.to_i64 - offset}
  end

  # :nodoc:
  # One aggregate per reduction level; records spill directly to disk. Gaps
  # are skipped rather than represented by millions of empty zoom windows.
  class ZoomBuilder
    getter reduction : Int64
    getter count : Int64 = 0_i64
    getter? valid = true
    getter spool : TempFile
    @tid : UInt32? = nil
    @first = 0_i64
    @last = 0_i64
    @aggregate = Aggregate.new
    @finished = false

    def initialize(@reduction : Int64)
      @spool = TempFile.new
    end

    def add(tid : UInt32, chrom_size : Int64, interval : BigWig::Interval) : Nil
      return unless @valid
      position = interval.start
      while position < interval.stop
        first = position // @reduction * @reduction
        last = Math.min(first + @reduction, chrom_size)
        if @tid != tid || @first != first
          flush_record
          return unless @valid
          @tid, @first, @last = tid, first, last
        end
        stop = Math.min(interval.stop, last)
        @aggregate.add(interval.value.to_f64, (stop - position).to_f64)
        position = stop
      end
    end

    private def flush_record : Nil
      tid = @tid
      return unless tid && @aggregate.count > 0
      begin
        minimum = @aggregate.min.to_f32
        maximum = @aggregate.max.to_f32
        sum = @aggregate.sum.to_f32
        squares = @aggregate.sum_squared.to_f32
      rescue OverflowError
        @valid = false
        return
      end
      # Zoom fields are float32 even though global summaries are float64.
      # Omit an unrepresentable level instead of emitting infinities.
      unless minimum.finite? && maximum.finite? && sum.finite? && squares.finite?
        @valid = false
        return
      end
      Internal.put(@spool.io, tid)
      Internal.put(@spool.io, @first.to_u32)
      Internal.put(@spool.io, @last.to_u32)
      Internal.put(@spool.io, @aggregate.count.to_u32)
      Internal.put(@spool.io, minimum)
      Internal.put(@spool.io, maximum)
      Internal.put(@spool.io, sum)
      Internal.put(@spool.io, squares)
      @count += 1
      @valid = false if @count > UInt32::MAX
      @aggregate = Aggregate.new
    end

    def finish : Nil
      return if @finished
      flush_record if @valid
      @finished = true
      @spool.io.flush
      @spool.io.rewind
    end

    def close : Nil
      @spool.close
    end
  end

  # :nodoc:
  module ZoomWriter
    def self.append(io : IO, zoom : ZoomBuilder, options : BigWig::WriteOptions) : ZoomLevel
      index = TempFile.new
      begin
        data_offset = io.pos.to_i64
        Internal.put(io, zoom.count.to_u32)
        buffer = Bytes.new(options.block_items * 32)
        remaining = zoom.count
        block_count = 0_i64
        zoom.spool.io.rewind
        while remaining > 0
          count = Math.min(remaining, options.block_items).to_i
          bytes = buffer[0, count * 32]
          zoom.spool.io.read_fully(bytes)
          first = Cursor.new(bytes)
          start_chrom, start = first.u32, first.u32.to_i64
          last = Cursor.new(bytes[(count - 1) * 32, 32])
          end_chrom = last.u32
          last.skip(4)
          stop = last.u32.to_i64
          offset, size = Internal.write_data_block(io, bytes, options)
          BlockRef.new(start_chrom, start, end_chrom, stop, offset, size).write(index.io)
          remaining -= count
          block_count += 1
        end
        index_offset = RTreeWriter.write(io, index.io, block_count,
          options.index_block_size, options.block_items)
        ZoomLevel.new(zoom.reduction, data_offset, index_offset)
      ensure
        index.close
      end
    end
  end
end

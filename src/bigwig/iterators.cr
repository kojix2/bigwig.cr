module BigWig::Internal
  # :nodoc:
  class IntervalIterator
    include Iterator(BigWig::Interval)

    @cursor = Cursor.new
    @items = 0
    @kind = 0_u8
    @tid = 0_u32
    @section_start = 0_i64
    @section_stop = 0_i64
    @span = 0_i64
    @step = 0_i64
    @position = 0_i64
    @previous_tid : UInt32? = nil
    @previous_stop = 0_i64
    @finished = false

    def initialize(@source : Source, @blocks : BlockIterator, @loader : BlockLoader,
                   @query : Region, @big_endian : Bool, @chrom_sizes : Hash(UInt32, Int64))
    end

    private def load_section : Bool
      if @cursor.remaining == 0
        ref = @blocks.next
        if ref.is_a?(Iterator::Stop)
          @finished = true
          return false
        end
        @cursor = Cursor.new(@loader.load(ref), @big_endian)
      end
      read_section
      true
    end

    private def read_interval : BigWig::Interval
      case @kind
      when 1
        first = @cursor.u32.to_i64
        last = @cursor.u32.to_i64
        value = @cursor.f32
      when 2
        first = @cursor.u32.to_i64
        last = first + @span
        value = @cursor.f32
      when 3
        first = @position
        last = first + @span
        value = @cursor.f32
        @position += @step
      else
        raise BigWig::FormatError.new("unknown bigWig section type")
      end
      @items -= 1
      unless @section_start <= first < last <= @section_stop
        raise BigWig::FormatError.new("interval outside its section bounds")
      end
      if (previous = @previous_tid) && (previous > @tid || (previous == @tid && first < @previous_stop))
        raise BigWig::FormatError.new("bigWig intervals are unsorted or overlapping")
      end
      @previous_tid = @tid
      @previous_stop = last
      BigWig::Interval.new(first, last, value)
    end

    def next : BigWig::Interval | Iterator::Stop
      raise BigWig::ClosedError.new("reader is closed") if @source.closed?
      return stop if @finished
      loop do
        if @items == 0
          return stop unless load_section
          next if @items == 0
        end
        interval = read_interval
        if @tid == @query.tid && interval.start < @query.stop && interval.stop > @query.start
          return interval
        end
      end
    end

    private def read_section : Nil
      @tid = @cursor.u32
      @section_start = @cursor.u32.to_i64
      @section_stop = @cursor.u32.to_i64
      @step = @cursor.u32.to_i64
      @span = @cursor.u32.to_i64
      @kind = @cursor.u8
      @cursor.skip(1)
      @items = @cursor.u16.to_i
      size = @chrom_sizes[@tid]? || raise BigWig::FormatError.new("section has an unknown chromosome ID")
      unless @section_start <= @section_stop <= size
        raise BigWig::FormatError.new("section exceeds chromosome bounds")
      end
      item_size = case @kind
                  when 1 then 12
                  when 2 then 8
                  when 3 then 4
                  else        raise BigWig::FormatError.new("unknown bigWig section type #{@kind}")
                  end
      if @kind != 1 && @span == 0
        raise BigWig::FormatError.new("zero-width fixed-span section")
      end
      if @kind == 3 && (@step == 0 || @step < @span)
        raise BigWig::FormatError.new("invalid fixedStep section")
      end
      if @items * item_size > @cursor.remaining
        raise BigWig::FormatError.new("truncated bigWig section")
      end
      @position = @section_start
    end
  end

  # :nodoc:
  class EntryIterator
    include Iterator(BigWig::BedEntry)
    @cursor = Cursor.new
    @finished = false

    def initialize(@source : Source, @blocks : BlockIterator, @loader : BlockLoader,
                   @query : Region, @big_endian : Bool, @chrom_sizes : Hash(UInt32, Int64))
    end

    def next : BigWig::BedEntry | Iterator::Stop
      raise BigWig::ClosedError.new("reader is closed") if @source.closed?
      return stop if @finished
      loop do
        if @cursor.remaining == 0
          ref = @blocks.next
          if ref.is_a?(Iterator::Stop)
            @finished = true
            return stop
          end
          @cursor = Cursor.new(@loader.load(ref), @big_endian)
        end
        tid = @cursor.u32
        first = @cursor.u32.to_i64
        last = @cursor.u32.to_i64
        rest = @cursor.cstring_bytes
        size = @chrom_sizes[tid]? || raise BigWig::FormatError.new("bigBed entry has an unknown chromosome ID")
        unless first <= last <= size
          raise BigWig::FormatError.new("bigBed entry exceeds chromosome bounds")
        end
        if tid == @query.tid && first < @query.stop && last > @query.start
          return BigWig::BedEntry.new(first, last, String.new(rest))
        end
      end
    end
  end

  # :nodoc:
  class ZoomIterator
    include Iterator(ZoomRecord)
    @cursor = Cursor.new
    @finished = false
    @record_tid = 0_u32
    @previous_tid : UInt32? = nil
    @previous_stop = 0_i64

    def initialize(@source : Source, @blocks : BlockIterator, @loader : BlockLoader,
                   @query : Region, @big_endian : Bool, @chrom_sizes : Hash(UInt32, Int64))
    end

    private def load_zoom_block : Bool
      return true unless @cursor.remaining == 0
      ref = @blocks.next
      if ref.is_a?(Iterator::Stop)
        @finished = true
        return false
      end
      bytes = @loader.load(ref)
      raise BigWig::FormatError.new("truncated zoom block") unless bytes.size % 32 == 0
      @cursor = Cursor.new(bytes, @big_endian)
      true
    end

    private def read_zoom_record : ZoomRecord
      tid = @cursor.u32
      @record_tid = tid
      first = @cursor.u32.to_i64
      last = @cursor.u32.to_i64
      count = @cursor.u32
      record = ZoomRecord.new(first, last, count, @cursor.f32, @cursor.f32, @cursor.f32, @cursor.f32)
      size = @chrom_sizes[tid]? || raise BigWig::FormatError.new("zoom record has an unknown chromosome ID")
      unless first < last <= size && count <= last - first
        raise BigWig::FormatError.new("invalid zoom record bounds/count")
      end
      if (previous = @previous_tid) && (previous > tid || (previous == tid && first < @previous_stop))
        raise BigWig::FormatError.new("zoom records are unsorted or overlapping")
      end
      @previous_tid = tid
      @previous_stop = last
      record
    end

    def next : ZoomRecord | Iterator::Stop
      raise BigWig::ClosedError.new("reader is closed") if @source.closed?
      return stop if @finished
      loop do
        return stop unless load_zoom_block
        record = read_zoom_record
        if @record_tid == @query.tid && record.start < @query.stop && record.stop > @query.start
          return record
        end
      end
    end
  end

  # :nodoc:
  class ValueIterator
    include Iterator(Float32)
    @current : BigWig::Interval? = nil
    @input_finished = false

    def initialize(@source : Source, @intervals : IntervalIterator,
                   @position : Int64, @stop : Int64)
    end

    def next : Float32 | Iterator::Stop
      raise BigWig::ClosedError.new("reader is closed") if @source.closed?
      return stop if @position >= @stop
      loop do
        if current = @current
          break if current.stop > @position
          @current = nil
        end
        break if @input_finished
        value = @intervals.next
        if value.is_a?(Iterator::Stop)
          @input_finished = true
        else
          @current = value
        end
      end
      result = Float32::NAN
      if current = @current
        result = current.value if current.start <= @position < current.stop
      end
      @position += 1
      result
    end
  end
end

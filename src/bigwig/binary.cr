module BigWig::Internal
  BW_MAGIC    = 0x888FFC26_u32
  BB_MAGIC    = 0x8789F2EB_u32
  BPT_MAGIC   = 0x78CA8C91_u32
  RTREE_MAGIC = 0x2468ACE0_u32
  MAX_COORD   = UInt32::MAX.to_i64
  MAX_DEPTH   = 64
  LE          = IO::ByteFormat::LittleEndian

  # A checked, allocation-free cursor over a block. Never casts an unaligned
  # file pointer to a native struct; byte order is explicit for every field.
  struct Cursor
    getter position : Int32 = 0

    def initialize(@bytes : Bytes = Bytes.empty, @big_endian : Bool = false)
    end

    def remaining : Int32
      @bytes.size - @position
    end

    def take(count : Int) : Bytes
      unless 0 <= count <= remaining
        raise BigWig::FormatError.new("truncated binary record at byte #{@position}")
      end
      result = @bytes[@position, count.to_i]
      @position += count.to_i
      result
    end

    def skip(count : Int) : Nil
      take(count)
    end

    {% for name, type in {u8: UInt8, u16: UInt16, u32: UInt32, u64: UInt64, f32: Float32, f64: Float64} %}
      def {{ name.id }} : {{ type }}
        bytes = take(sizeof({{ type }}))
        if @big_endian
          IO::ByteFormat::BigEndian.decode({{ type }}, bytes)
        else
          IO::ByteFormat::LittleEndian.decode({{ type }}, bytes)
        end
      end
    {% end %}

    def offset : Int64
      value = u64
      raise BigWig::FormatError.new("file offset exceeds Int64") if value > Int64::MAX
      value.to_i64
    end

    def cstring_bytes : Bytes
      rest = @bytes[@position, remaining]
      length = rest.index(0_u8) || raise BigWig::FormatError.new("unterminated bigBed string")
      value = take(length)
      skip(1)
      value
    end
  end

  def self.coordinate(value : Int, name : String) : Int64
    unless 0 <= value <= MAX_COORD
      raise ArgumentError.new("#{name} must be in 0..#{MAX_COORD}")
    end
    value.to_i64
  end

  def self.range_bounds(range : Range(B, E)) : Tuple(Int64, Int64) forall B, E
    {% unless B < Int && E < Int %}
      {% raise "BigWig queries require a bounded integer Range" %}
    {% end %}
    start = coordinate(range.begin, "start")
    stop = coordinate(range.end, "stop")
    unless range.excludes_end?
      raise ArgumentError.new("inclusive range end exceeds the coordinate limit") if stop == MAX_COORD
      stop += 1
    end
    {start, stop}
  end

  def self.zeros(io : IO, count : Int64) : Nil
    buffer = uninitialized UInt8[4096]
    buffer.to_slice.fill(0_u8)
    left = count
    while left > 0
      size = Math.min(left, buffer.size).to_i
      io.write(buffer.to_slice[0, size])
      left -= size
    end
  end

  def self.put(io : IO, value) : Nil
    io.write_bytes(value, LE)
  end

  def self.patch(io : IO, offset : Int64, value) : Nil
    previous = io.pos
    begin
      io.seek(offset)
      put(io, value)
    ensure
      io.seek(previous)
    end
  end

  class Chromosome
    getter name : String
    getter id : UInt32
    getter size : Int64

    def initialize(@name, @id, @size)
    end
  end

  class Region
    getter tid : UInt32
    getter start : Int64
    getter stop : Int64

    def initialize(@tid, @start, @stop)
    end
  end

  class ZoomLevel
    getter reduction : Int64
    getter data_offset : Int64
    getter index_offset : Int64

    def initialize(@reduction, @data_offset, @index_offset)
    end
  end

  # Value type: R-tree traversal and construction create one descriptor per
  # visited block, so storing these inline avoids hot-path heap allocations.
  record BlockRef,
    start_chrom : UInt32, start : Int64,
    end_chrom : UInt32, stop : Int64,
    offset : Int64, size : Int64 do
    def overlaps?(region : Region) : Bool
      return false if region.start == region.stop
      return false if region.tid < start_chrom || region.tid > end_chrom
      return false if region.tid == start_chrom && start >= region.stop
      return false if region.tid == end_chrom && stop <= region.start
      true
    end

    def write(io : IO) : Nil
      Internal.put(io, start_chrom)
      Internal.put(io, start.to_u32)
      Internal.put(io, end_chrom)
      Internal.put(io, stop.to_u32)
      Internal.put(io, offset.to_u64)
      Internal.put(io, size.to_u64)
    end

    def self.read(io : IO) : self
      bytes = uninitialized UInt8[32]
      io.read_fully(bytes.to_slice)
      c = Cursor.new(bytes.to_slice)
      new(c.u32, c.u32.to_i64, c.u32, c.u32.to_i64, c.offset, c.offset)
    end
  end

  # Value type: approximate statistics stream one instance per zoom record.
  record ZoomRecord,
    start : Int64, stop : Int64, count : UInt32,
    min : Float32, max : Float32, sum : Float32, sum_squared : Float32
end

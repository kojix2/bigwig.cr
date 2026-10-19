require "option_parser"
require "compress/zlib"

# Build fixtures directly from the file format, without using bigwig.cr.
module Fixtures
  BW    = 0x888ffc26_u32
  BB    = 0x8789f2eb_u32
  BPT   = 0x78ca8c91_u32
  RTREE = 0x2468ace0_u32
  MAX   = UInt32::MAX.to_i64

  alias Chrom = Tuple(String, Int64)
  alias Block = Tuple(Int64, Int64, Int64, Int64, Int64, Int64)
  alias ChromEntry = Tuple(String, Int64, Int64)

  record Row, start : Int64, stop : Int64, value : Float64
  record Group, tid : Int64, kind : Int32, rows : Array(Row), span : Int64 = 0_i64, step : Int64 = 0_i64
  record BedRow, start : Int64, stop : Int64, rest : String
  record ZoomRow, tid : Int64, start : Int64, stop : Int64, valid : Int64,
    min : Float64, max : Float64, sum : Float64, sum_squares : Float64
  record Zoom, reduction : Int64, rows : Array(ZoomRow)
  record RawGroup, tid : Int64, start : Int64, end_tid : Int64, stop : Int64, bytes : Bytes

  class ChromNode
    getter entries : Array(ChromEntry)?
    getter children : Array(ChromNode)?

    def initialize(@entries : Array(ChromEntry)? = nil, @children : Array(ChromNode)? = nil)
    end

    def first_key : String
      if entries = @entries
        entries.first[0]
      else
        @children.not_nil!.first.first_key
      end
    end
  end

  class IndexNode
    getter blocks : Array(Block)?
    getter children : Array(IndexNode)?

    def initialize(@blocks : Array(Block)? = nil, @children : Array(IndexNode)? = nil)
    end
  end

  class Builder
    @io = IO::Memory.new
    @format : IO::ByteFormat

    def initialize(big_endian = false)
      @format = big_endian ? IO::ByteFormat::BigEndian : IO::ByteFormat::LittleEndian
    end

    def position : Int64
      @io.pos.to_i64
    end

    def bytes : Bytes
      @io.to_slice.dup
    end

    def raw(data : Bytes) : Int64
      start = position
      @io.write(data)
      start
    end

    def reserve(count : Int) : Int64
      raw(Bytes.new(count))
    end

    def patch(at : Int64, &)
      saved = position
      @io.seek(at)
      yield
      @io.seek(saved)
    end

    def u8(value : Int) : Nil
      @io.write_byte(value.to_u8)
    end

    def u16(value : Int) : Nil
      @io.write_bytes(value.to_u16, @format)
    end

    def u32(value : Int) : Nil
      @io.write_bytes(value.to_u32, @format)
    end

    def u64(value : Int) : Nil
      @io.write_bytes(value.to_u64, @format)
    end

    def f32(value : Float) : Nil
      @io.write_bytes(value.to_f32, @format)
    end

    def f64(value : Float) : Nil
      @io.write_bytes(value.to_f64, @format)
    end

    def string(value : String) : Nil
      @io.write(value.to_slice)
    end

    def chromosome_tree(chroms : Array(Chrom), fanout = 2) : Int64
      key_size = chroms.max_of { |name, _| name.bytesize } + 1
      start = position
      u32(BPT)
      u32(fanout)
      u32(key_size)
      u32(8)
      u64(chroms.size)
      u64(0)

      entries = chroms.map_with_index { |(name, size), id| {name, id.to_i64, size} }.sort_by(&.[0])
      nodes = entries.each_slice(fanout).map { |part| ChromNode.new(entries: part) }.to_a
      while nodes.size > 1
        nodes = nodes.each_slice(fanout).map { |part| ChromNode.new(children: part) }.to_a
      end
      emit_chrom(nodes.first, key_size)
      start
    end

    private def emit_chrom(node : ChromNode, key_size : Int32) : Int64
      leaf = !!node.entries
      count = leaf ? node.entries.not_nil!.size : node.children.not_nil!.size
      start = reserve(4 + count * (key_size + 8))
      patch(start) do
        u8(leaf ? 1 : 0)
        u8(0)
        u16(count)
      end
      count.times do |i|
        entry = node.entries.try(&.[i])
        child = node.children.try(&.[i])
        key = entry ? entry[0] : child.not_nil!.first_key
        offset = child ? emit_chrom(child, key_size) : 0_i64
        patch(start + 4 + i * (key_size + 8)) do
          string(key)
          reserve(key_size - key.bytesize)
          if entry
            u32(entry[1])
            u32(entry[2])
          else
            u64(offset)
          end
        end
      end
      start
    end

    def index(blocks : Array(Block), fanout = 2, slot = 1) : Int64
      start = reserve(48)
      nodes = blocks.each_slice(fanout).map { |part| IndexNode.new(blocks: part) }.to_a
      nodes << IndexNode.new(blocks: [] of Block) if nodes.empty?
      while nodes.size > 1
        nodes = nodes.each_slice(fanout).map { |part| IndexNode.new(children: part) }.to_a
      end
      bounds = emit_index(nodes.first)
      patch(start) do
        u32(RTREE)
        u32(fanout)
        u64(blocks.size)
        u32(bounds[0])
        u32(bounds[1])
        u32(bounds[2])
        u32(bounds[3])
        u64(start)
        u32(slot)
        u32(0)
      end
      start
    end

    private def emit_index(node : IndexNode) : Block
      leaf = !!node.blocks
      count = leaf ? node.blocks.not_nil!.size : node.children.not_nil!.size
      width = leaf ? 32 : 24
      start = reserve(4 + count * width)
      patch(start) do
        u8(leaf ? 1 : 0)
        u8(0)
        u16(count)
      end
      descriptors = node.blocks || node.children.not_nil!.map { |child| emit_index(child) }
      descriptors.each_with_index do |desc, i|
        patch(start + 4 + i * width) do
          4.times { |field| u32(desc[field]) }
          u64(desc[4])
          u64(desc[5]) if leaf
        end
      end
      if descriptors.empty?
        {0_i64, 0_i64, 0_i64, 0_i64, start, 0_i64}
      else
        {descriptors.first[0], descriptors.first[1], descriptors.last[2], descriptors.last[3], start, 0_i64}
      end
    end
  end

  def self.compressed(data : Bytes) : Bytes
    output = IO::Memory.new
    Compress::Zlib::Writer.open(output) { |writer| writer.write(data) }
    output.to_slice.dup
  end

  def self.section(group : Group, big_endian : Bool) : Bytes
    b = Builder.new(big_endian)
    b.u32(group.tid)
    b.u32(group.rows.first.start)
    b.u32(group.rows.last.stop)
    b.u32(group.step)
    b.u32(group.span)
    b.u8(group.kind)
    b.u8(0)
    b.u16(group.rows.size)
    group.rows.each do |row|
      b.u32(row.start) unless group.kind == 3
      b.u32(row.stop) if group.kind == 1
      b.f32(row.value)
    end
    b.bytes
  end

  def self.make_file(chroms : Array(Chrom), groups : Array(Group), *,
                     big_endian = false, compress = true, zooms = [] of Zoom,
                     bed_rows : Array(BedRow)? = nil, sql : String? = nil,
                     multi_section = false) : Bytes
    b = Builder.new(big_endian)
    b.reserve(64 + 24 * zooms.size)
    sql_pos = 0_i64
    if text = sql
      sql_pos = b.position
      b.string(text)
      b.u8(0)
    end
    summary_pos = b.reserve(40)
    chrom_pos = b.chromosome_tree(chroms)
    data_pos = b.position
    b.u64(bed_rows ? 1 : groups.size)

    raw_groups = [] of RawGroup
    if rows = bed_rows
      data = Builder.new(big_endian)
      rows.each do |row|
        data.u32(0)
        data.u32(row.start)
        data.u32(row.stop)
        data.string(row.rest)
        data.u8(0)
      end
      raw_groups << RawGroup.new(0_i64, rows.first.start, 0_i64, rows.last.stop, data.bytes)
    else
      groups.each do |group|
        raw_groups << RawGroup.new(group.tid, group.rows.first.start, group.tid,
          group.rows.last.stop, section(group, big_endian))
      end
    end
    if multi_section
      merged = IO::Memory.new
      raw_groups.each { |item| merged.write(item.bytes) }
      raw_groups = [RawGroup.new(raw_groups.first.tid, raw_groups.first.start,
        raw_groups.last.end_tid, raw_groups.last.stop, merged.to_slice.dup)]
    end

    blocks = [] of Block
    largest = 0
    raw_groups.each do |item|
      encoded = compress ? compressed(item.bytes) : item.bytes
      offset = b.raw(encoded)
      blocks << {item.tid, item.start, item.end_tid, item.stop, offset, encoded.size.to_i64}
      largest = Math.max(largest, item.bytes.size)
    end
    index_pos = b.index(blocks)

    zoom_headers = [] of Tuple(Int64, Int64, Int64)
    zooms.each do |zoom|
      data_pos_zoom = b.position
      b.u32(zoom.rows.size)
      zoom_blocks = [] of Block
      zoom.rows.each_slice(2) do |chunk|
        data = Builder.new(big_endian)
        chunk.each do |row|
          data.u32(row.tid)
          data.u32(row.start)
          data.u32(row.stop)
          data.u32(row.valid)
          data.f32(row.min)
          data.f32(row.max)
          data.f32(row.sum)
          data.f32(row.sum_squares)
        end
        raw = data.bytes
        encoded = compress ? compressed(raw) : raw
        offset = b.raw(encoded)
        zoom_blocks << {chunk.first.tid, chunk.first.start, chunk.last.tid,
                        chunk.last.stop, offset, encoded.size.to_i64}
        largest = Math.max(largest, raw.size)
      end
      zoom_headers << {zoom.reduction, data_pos_zoom, b.index(zoom_blocks, slot: 2)}
    end
    b.u32(bed_rows ? BB : BW)

    covered = 0_i64
    minimum = Float64::INFINITY
    maximum = -Float64::INFINITY
    sum = 0.0
    squares = 0.0
    if rows = bed_rows
      rows.each { |row| covered += row.stop - row.start }
      minimum = maximum = 1.0
    else
      groups.each do |group|
        group.rows.each do |row|
          width = row.stop - row.start
          covered += width
          minimum = Math.min(minimum, row.value)
          maximum = Math.max(maximum, row.value)
          sum += width * row.value
          squares += width * row.value * row.value
        end
      end
    end
    b.patch(summary_pos) do
      b.u64(covered)
      b.f64(minimum)
      b.f64(maximum)
      b.f64(sum)
      b.f64(squares)
    end
    b.patch(0_i64) do
      b.u32(bed_rows ? BB : BW)
      b.u16(4)
      b.u16(zooms.size)
      b.u64(chrom_pos)
      b.u64(data_pos)
      b.u64(index_pos)
      b.u16(bed_rows ? 6 : 0)
      b.u16(bed_rows ? 3 : 0)
      b.u64(sql_pos)
      b.u64(summary_pos)
      b.u32(compress ? largest : 0)
      b.u64(0)
    end
    zoom_headers.each_with_index do |(reduction, data, index), i|
      b.patch(64_i64 + 24_i64 * i) do
        b.u32(reduction)
        b.u32(0)
        b.u64(data)
        b.u64(index)
      end
    end
    b.bytes
  end

  def self.fixtures : Hash(String, Bytes)
    chroms = [{"zchr", 1000_i64}, {"achr", MAX}, {"empty", 10_i64}]
    groups = [
      Group.new(0, 1, [Row.new(0, 4, 1.0), Row.new(8, 10, 2.0)]),
      Group.new(0, 2, [Row.new(20, 22, 3.0), Row.new(25, 27, -1.0)], 2),
      Group.new(0, 3, [Row.new(40, 42, 0.5), Row.new(43, 45, 1.5),
                       Row.new(46, 48, 2.5)], 2, 3),
      Group.new(1, 1, [Row.new(3_000_000_000_i64, 3_000_000_003_i64, 4.0)]),
      Group.new(1, 1, [Row.new(MAX - 2, MAX, 5.0)]),
    ]
    result = Hash(String, Bytes).new
    {false, true}.each do |big_endian|
      endian_name = big_endian ? "big" : "little"
      {false, true}.each do |compress|
        compression_name = compress ? "zlib" : "raw"
        result["mixed-#{endian_name}-#{compression_name}.bw"] =
          make_file(chroms, groups, big_endian: big_endian, compress: compress)
      end
    end
    result["multi-section.bw"] = make_file(chroms, groups.first(3), multi_section: true)
    zoom = Zoom.new(16, [
      ZoomRow.new(0, 0, 16, 16, 2.0, 2.0, 32.0, 64.0),
      ZoomRow.new(0, 32, 48, 16, 4.0, 4.0, 64.0, 256.0),
    ])
    result["zoom-positive.bw"] = make_file([{"chr", 128_i64}],
      [Group.new(0, 1, [Row.new(0, 16, 2.0), Row.new(32, 48, 4.0)])], zooms: [zoom])
    sql = "table genes\n\"Test annotations\"\n(\nstring chrom; \"Reference sequence\"\nuint chromStart; \"Start, zero-based\"\nuint chromEnd; \"End, excluded\"\nstring name; \"Name\"\nuint score; \"Score\"\nchar[1] strand; \"Strand\"\n)\n"
    bed_rows = [
      BedRow.new(1, 4, "geneA\t10\t+"),
      BedRow.new(4, 4, "insertion\t0\t."),
      BedRow.new(5, 9, "geneB\t20\t-"),
    ]
    {false, true}.each do |big_endian|
      name = big_endian ? "big" : "little"
      result["minimal-#{name}.bb"] =
        make_file([{"chr", 100_i64}], [] of Group, big_endian: big_endian,
          bed_rows: bed_rows, sql: sql)
    end
    result
  end
end

write = false
check = false
target = File.expand_path("../../spec/fixtures", __DIR__)
parser = OptionParser.new do |options|
  options.banner = "Usage: crystal run scripts/fixtures/generate.cr -- [--write] [--check] [--output DIR]"
  options.on("--write", "Write generated fixtures") { write = true }
  options.on("--check", "Compare generated bytes with existing fixtures") { check = true }
  options.on("--output DIR", "Fixture directory (default: spec/fixtures)") { |path| target = path }
  options.on("-h", "--help", "Show this help") { puts options; exit 0 }
end

begin
  args = ARGV.dup
  parser.parse(args)
  raise ArgumentError.new("unexpected arguments: #{args.join(" ")}") unless args.empty?
  fixtures = Fixtures.fixtures
  Dir.mkdir_p(target) if write
  fixtures.each do |name, data|
    path = File.join(target, name)
    File.write(path, data) if write
    if check
      actual = File.read(path).to_slice
      raise "fixture differs: #{path}" unless actual == data
    end
  end
  puts "#{fixtures.size} deterministic fixtures #{check ? "verified" : write ? "generated" : "available"}"
rescue ex : OptionParser::Exception | ArgumentError
  STDERR.puts "generate: #{ex.message}"
  exit 2
rescue ex
  STDERR.puts "generate: #{ex.message}"
  exit 1
end

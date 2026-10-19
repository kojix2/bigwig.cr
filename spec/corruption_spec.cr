require "./spec_helper"

describe "malformed input and resource limits" do
  it "rejects short headers and non-BBI magic" do
    {0, 3, 63}.each do |size|
      expect_raises(BigWig::FormatError) { BigWig.open(IO::Memory.new(Bytes.new(size))) }
    end
    expect_raises(BigWig::FormatError) { BigWig.open(IO::Memory.new(Bytes.new(64))) }
  end

  it "rejects invalid header offsets and versions" do
    io = fixture_io("mixed-little-raw.bw")
    patch(io, 8, UInt64::MAX)
    expect_raises(BigWig::FormatError) { BigWig.open(io) }
    io.closed?.should be_false
    io = fixture_io("mixed-little-raw.bw")
    patch(io, 4, 99_u16)
    expect_raises(BigWig::UnsupportedError) { BigWig.open(io) }
  end

  it "rejects excessive decompression reservations before allocating" do
    io = fixture_io("mixed-little-zlib.bw")
    patch(io, 52, UInt32::MAX)
    expect_raises(BigWig::LimitError) { BigWig.open(io) }
  end

  it "rejects a corrupted compressed stream rather than returning no data" do
    io = fixture_io("mixed-little-zlib.bw")
    data = get_u64(io, 16).to_i + 8
    patch(io, data, 0_u32)
    BigWig.open(io) do |reader|
      expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
    end
  end

  it "rejects a block that expands past the declared header buffer" do
    io = fixture_io("mixed-little-zlib.bw")
    patch(io, 52, 25_u32)
    BigWig.open(io) do |reader|
      expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
    end
  end

  it "rejects malformed and cyclic index nodes" do
    io = fixture_io("mixed-little-raw.bw")
    index = get_u64(io, 24).to_i
    patch(io, index + 48, 2_u8)
    BigWig.open(io) do |reader|
      expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
    end
    io = fixture_io("mixed-little-raw.bw")
    index = get_u64(io, 24).to_i
    root = index + 48
    patch(io, root + 4 + 16, root.to_u64)
    BigWig.open(io) do |reader|
      expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
    end
  end

  it "rejects a cyclic chromosome tree" do
    io = fixture_io("mixed-little-raw.bw")
    tree = get_u64(io, 8).to_i
    key_size = IO::ByteFormat::LittleEndian.decode(UInt32, io.to_slice[tree + 8, 4]).to_i
    root = tree + 32
    patch(io, root + 4 + key_size, root.to_u64)
    expect_raises(BigWig::FormatError) { BigWig.open(io) }
  end

  it "rejects invalid section types, item counts, and overlapping data" do
    {1, 2, 3}.each do |mutation|
      io = fixture_io("mixed-little-raw.bw")
      data = get_u64(io, 16).to_i + 8
      case mutation
      when 1
        patch(io, data + 20, 9_u8)
      when 2
        patch(io, data + 22, 1000_u16)
      when 3
        patch(io, data + 24 + 12, 3_u32)
      end
      BigWig.open(io) do |reader|
        expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
      end
    end
  end

  it "detects a file truncated after its header was read" do
    with_output_path do |path|
      File.copy(fixture("mixed-little-raw.bw"), path)
      File.open(path, "r+b") do |file|
        reader = BigWig.open(file)
        begin
          file.truncate(100)
          expect_raises(BigWig::FormatError) { reader.intervals("zchr") }
        ensure
          reader.close
        end
      end
    end
  end

  it "checks limits and closes explicitly owned IO on initialization failure" do
    expect_raises(ArgumentError) { BigWig::ReadOptions.new(cache_bytes: 0).validate }
    expect_raises(ArgumentError) { BigWig::WriteOptions.new(block_items: 0).validate }
    expect_raises(ArgumentError) { BigWig::WriteOptions.new(index_block_size: 1).validate }
    expect_raises(ArgumentError) { BigWig::WriteOptions.new(max_zooms: -1).validate }
    io = IO::Memory.new(Bytes.new(8))
    expect_raises(BigWig::FormatError) { BigWig.open(io, sync_close: true) }
    io.closed?.should be_true
  end
end

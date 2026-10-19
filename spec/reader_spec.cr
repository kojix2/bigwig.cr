require "./spec_helper"

describe BigWig::Reader do
  it "opens by signature and returns the block value" do
    result = BigWig.open(fixture("test.bw")) do |reader|
      reader.bigwig?.should be_true
      reader.bigbed?.should be_false
      reader.header.version.should eq(4)
      reader.header.bases_covered.should eq(154_u64)
      reader.header.has_summary?.should be_true
      reader.chromosomes.should eq({"1" => 195471971_i64, "10" => 130694993_i64})
      reader.chromosome_size("missing").should be_nil
      17
    end
    result.should eq(17)
  end

  it "does not expose mutable chromosome metadata" do
    BigWig.open(fixture("test.bw")) do |reader|
      reader.chromosomes.clear
      reader.chromosome_size("1").should eq(195471971_i64)
    end
  end

  it "preserves original overlapping interval boundaries" do
    BigWig.open(fixture("test.bw")) do |reader|
      reader.intervals("1", 101, 102).should eq([BigWig::Interval.new(100, 150, 1.4)])
      reader.intervals("1", 3, 100).should be_empty
      reader.intervals("1", 150, 151).should eq([BigWig::Interval.new(150, 151, 1.5)])
    end
  end

  it "supports bounded inclusive and exclusive integer ranges" do
    BigWig.open(fixture("test.bw")) do |reader|
      reader.values("1", 0...3).should eq([0.1_f32, 0.2_f32, 0.3_f32])
      reader.values("1", 0..2).should eq(reader.values("1", 0, 3))
      reader.each_interval("1", 0...3).to_a.should eq(reader.intervals("1", 0, 3))
      seen = [] of BigWig::Interval
      reader.each_interval("1", 0..2) { |item| seen << item }
      seen.size.should eq(3)
      expect_raises(ArgumentError) { reader.values("1", 3..1) }
    end
  end

  it "fills dense values and reusable buffers with NaN gaps" do
    BigWig.open(fixture("test.bw")) do |reader|
      expected = [Float32::NAN, 1.4_f32, 1.4_f32, 1.4_f32]
      same_values(reader.values("1", 99, 103), expected)
      buffer = Slice(Float32).new(4, 99_f32)
      reader.read_values("1", 99, buffer).should eq(4)
      same_values(buffer, expected)
      same_values(reader.each_value("1", 99...103), expected)
      seen = [] of Float32
      reader.each_value("1", 99, 103) { |value| seen << value }
      same_values(seen, expected)
      reader.value("1", 3).nan?.should be_true
      reader.value("1", 0).should eq(0.1_f32)
    end
  end

  it "supports empty ranges, but not missing chromosomes or out-of-bounds ranges" do
    BigWig.open(fixture("test.bw")) do |reader|
      reader.values("1", 10, 10).should be_empty
      reader.each_interval("1", 10, 10).next.should be_a(Iterator::Stop)
      reader.each_value("1", 10, 10).next.should be_a(Iterator::Stop)
      expect_raises(ArgumentError) { reader.intervals("missing", 0, 1) }
      expect_raises(ArgumentError) { reader.values("1", -1, 1) }
      expect_raises(ArgumentError) { reader.intervals("1", 195471970, 195471972) }
      expect_raises(ArgumentError) { reader.value("1", 195471971) }
      reader.read_values("1", 195471971, Slice(Float32).empty).should eq(0)
    end
  end

  it "keeps separate cursors independent and returned records valid" do
    BigWig.open(fixture("mixed-little-zlib.bw")) do |reader|
      a = reader.each_interval("zchr")
      b = reader.each_interval("achr")
      saved = a.next.as(BigWig::Interval)
      b.next.as(BigWig::Interval).start.should eq(3_000_000_000_i64)
      a.next.as(BigWig::Interval).start.should eq(8)
      b.next.as(BigWig::Interval).stop.should eq(UInt32::MAX.to_i64)
      a.to_a.size.should eq(5)
      saved.should eq(BigWig::Interval.new(0, 4, 1))
    end
  end

  it "invalidates unfinished iterators when explicitly closed" do
    reader = BigWig::Reader.new(fixture("test.bw"))
    iterator = reader.each_interval("1")
    saved = iterator.next.as(BigWig::Interval)
    reader.close
    reader.close
    reader.closed?.should be_true
    saved.value.should eq(0.1_f32)
    expect_raises(BigWig::ClosedError) { iterator.next }
    expect_raises(BigWig::ClosedError) { reader.header }
    expect_raises(BigWig::ClosedError) { reader.intervals("1") }
  end

  it "borrows IO by default and preserves its position" do
    io = fixture_io("test.bw")
    io.seek(7)
    BigWig.open(io) do |reader|
      reader.value("1", 1).should eq(0.2_f32)
      io.pos.should eq(7)
    end
    io.closed?.should be_false
    io.pos.should eq(7)
    BigWig.open(io, sync_close: true, &.value("1", 0))
    io.closed?.should be_true
  end

  it "preserves a borrowed File position through positioned reads" do
    File.open(fixture("test.bw")) do |file|
      file.seek(9)
      BigWig.open(file, &.intervals("1").size.should(eq(5)))
      file.pos.should eq(9)
      file.closed?.should be_false
    end
  end

  it "limits eager allocations without restricting lazy queries" do
    options = BigWig::ReadOptions.new(max_array_bytes: 8_i64)
    BigWig.open(fixture("mixed-little-raw.bw"), options: options) do |reader|
      expect_raises(BigWig::LimitError) { reader.values("zchr", 0, 3) }
      expect_raises(BigWig::LimitError) { reader.intervals("zchr") }
      expect_raises(BigWig::LimitError) { reader.stats("zchr", bins: 2) }
      reader.each_interval("zchr").to_a.size.should eq(7)
    end
  end

  %w[mixed-little-raw.bw mixed-little-zlib.bw mixed-big-raw.bw mixed-big-zlib.bw].each do |name|
    it "decodes all section encodings and multilevel trees: #{name}" do
      BigWig.open(fixture(name)) do |reader|
        reader.intervals("zchr").should eq(mixed_intervals)
        reader.chromosomes.keys.should eq(["zchr", "achr", "empty"])
        reader.chromosome_size("achr").should eq(UInt32::MAX.to_i64)
        reader.value("achr", UInt32::MAX.to_i64 - 1).should eq(5_f32)
        reader.intervals("empty").should be_empty
        expect_raises(ArgumentError) { reader.intervals("achr", 0_i64, 4_294_967_296_i64) }
        expect_raises(ArgumentError) { reader.intervals("achr", 0_i64..UInt32::MAX.to_i64) }
      end
    end
  end

  it "decodes more than one section inside a compressed block" do
    BigWig.open(fixture("multi-section.bw")) do |reader|
      reader.intervals("zchr").should eq(mixed_intervals)
    end
  end
end

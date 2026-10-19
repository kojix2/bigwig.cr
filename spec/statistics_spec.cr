require "./spec_helper"

private def dense_stat(values : Array(Float32), type : Symbol) : Float64
  present = values.reject(&.nan?).map(&.to_f64)
  return Float64::NAN if present.empty?
  case type
  when :mean
    present.sum / present.size
  when :sum
    present.sum
  when :min
    present.min
  when :max
    present.max
  when :coverage
    present.size.to_f64 / values.size
  when :std
    return 0.0 if present.size == 1
    mean = present.sum / present.size
    Math.sqrt(present.sum { |value| (value - mean) ** 2 } / (present.size - 1))
  else
    raise "unknown test statistic"
  end
end

describe "statistics" do
  it "matches original fixture expectations" do
    BigWig.open(fixture("test.bw")) do |reader|
      reader.mean("1", 0, 3).should be_close(0.2000000054637591, 1e-12)
      reader.max("1", 0, 3).should eq(0.3_f32.to_f64)
      reader.min("1", 0, 3).should eq(0.1_f32.to_f64)
      reader.sum("1", 0, 3, exact: true).should be_close(0.6000000163912773, 1e-12)
      reader.coverage("1", 0, 100, exact: true).should eq(0.03)
      reader.mean("1").should be_close(1.3351851569281683, 1e-7)
      reader.std("1", 0, 1, exact: true).should eq(0.0)
    end
  end

  it "uses covered bases, not zeros for gaps" do
    BigWig.open(fixture("mixed-little-raw.bw")) do |reader|
      reader.mean("zchr", 0...10, exact: true).should be_close(8.0 / 6, 1e-12)
      reader.sum("zchr", 0...10, exact: true).should eq(8)
      reader.coverage("zchr", 0...10, exact: true).should eq(0.6)
      reader.min("zchr", 0...30).should eq(-1)
      reader.max("zchr", 0...30).should eq(3)
      reader.std("zchr", 0...10, exact: true).should be_close(Math.sqrt(4.0 / 15), 1e-12)
      reader.stats("empty", type: :coverage).first.nan?.should be_true
      reader.stats("empty", type: :sum).first.nan?.should be_true
    end
  end

  it "gives zero-width bins NaN with deterministic integer boundaries" do
    BigWig.open(fixture("test.bw")) do |reader|
      values = reader.stats("1", 0, 2, bins: 4, exact: true)
      values.size.should eq(4)
      values[0].nan?.should be_true
      values[1].should eq(0.1_f32.to_f64)
      values[2].nan?.should be_true
      values[3].should eq(0.2_f32.to_f64)
      reader.stats("1", 10, 10, bins: 3).all?(&.nan?).should be_true
      expect_raises(ArgumentError) { reader.stats("1", bins: 0) }
      expect_raises(ArgumentError) { reader.stats("1", type: :median) }
    end
  end

  it "computes zoom sum as an integral, not the C reference's sum of means" do
    BigWig.open(fixture("zoom-positive.bw")) do |reader|
      reader.zoom_levels.should eq([16_i64])
      reader.sum("chr").should eq(96)
      reader.sum("chr", exact: true).should eq(96)
      reader.mean("chr").should eq(3)
      reader.coverage("chr").should eq(0.25)
      reader.std("chr").should be_close(Math.sqrt(32.0 / 31), 1e-12)
      reader.sum("chr", 8, 72).should eq(80)
      reader.mean("chr", 8, 72).should be_close(80.0 / 24, 1e-12)
      bins = reader.stats("chr", 0...128, type: :sum, bins: 2)
      bins[0].should eq(96)
      bins[1].nan?.should be_true
    end
  end

  it "matches a dense oracle across randomized ranges and bins" do
    random = Random.new(9042)
    BigWig.open(fixture("mixed-little-zlib.bw")) do |reader|
      40.times do
        first = random.rand(0..80)
        last = random.rand((first + 1)..100)
        {1, 3, 11}.each do |bins|
          {:mean, :min, :max, :sum, :coverage, :std}.each do |type|
            actual = reader.stats("zchr", first, last, bins: bins, type: type, exact: true)
            bins.times do |i|
              lo = first + (last - first) * i // bins
              hi = first + (last - first) * (i + 1) // bins
              expected = dense_stat(reader.values("zchr", lo, hi), type)
              if expected.nan?
                actual[i].nan?.should be_true
              else
                actual[i].should be_close(expected, 1e-10)
              end
            end
          end
        end
      end
    end
  end

  it "keeps exact variance stable for a large baseline" do
    io = IO::Memory.new
    BigWig.create(io, {"chr" => 2}, options: BigWig::WriteOptions.new(max_zooms: 0)) do |writer|
      writer.add("chr", 0, 1, 100_000_000.0)
      writer.add("chr", 1, 2, 100_000_016.0)
    end
    BigWig.open(io) do |reader|
      reader.std("chr", exact: true).should be_close(Math.sqrt(128.0), 1e-9)
    end
  end
end

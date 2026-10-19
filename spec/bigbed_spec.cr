require "./spec_helper"

describe "bigBed" do
  it "reads the original Ruby binding fixture and AutoSQL" do
    BigWig.open(fixture("test.bigBed")) do |reader|
      reader.bigbed?.should be_true
      reader.bigwig?.should be_false
      entries = reader.entries("chr1", 10_000_000, 10_020_000)
      entries.size.should eq(3)
      entries.first.should eq(BigWig::BedEntry.new(10009333_i64, 10009640_i64, "61035\t130\t-\t0.026\t0.42\t404"))
      entries.last.stop.should eq(10024307)
      reader.sql.try(&.includes?("RnaElements")).should be_true
      reader.each_entry("chr1", 10_000_000...10_020_000).to_a.should eq(entries)
      expect_raises(BigWig::UnsupportedError) { reader.values("chr1", 0, 10) }
      expect_raises(BigWig::UnsupportedError) { reader.stats("chr1", 0, 10) }
    end
  end

  %w[minimal-little.bb minimal-big.bb].each do |name|
    it "reads #{name}, including interior zero-width entries" do
      reader = BigWig.open(fixture(name))
      begin
        iterator = reader.each_entry("chr", 0, 10)
        first = iterator.next.as(BigWig::BedEntry)
        iterator.to_a.size.should eq(2)
        reader.close
        first.rest.should eq("geneA\t10\t+")
      ensure
        reader.close
      end
      BigWig.open(fixture(name)) do |bigbed_reader|
        bigbed_reader.entries("chr", 4, 5).should be_empty
        bigbed_reader.entries("chr", 3, 5).map(&.rest).should eq(["geneA\t10\t+", "insertion\t0\t."])
        bigbed_reader.entries("chr", 4, 4).should be_empty
        count = 0
        bigbed_reader.each_entry("chr", 0..9) { |_| count += 1 }
        count.should eq(3)
      end
    end
  end

  it "returns nil for absent AutoSQL" do
    io = fixture_io("minimal-little.bb")
    patch(io, 36, 0_u64)
    BigWig.open(io, &.sql.should(be_nil))
  end

  it "rejects bigBed operations on bigWig" do
    BigWig.open(fixture("test.bw")) do |reader|
      expect_raises(BigWig::UnsupportedError) { reader.entries("1") }
      expect_raises(BigWig::UnsupportedError) { reader.sql }
    end
  end
end

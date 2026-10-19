require "./spec_helper"

describe BigWig::Writer do
  {false, true}.each do |compressed|
    {0, 4}.each do |zooms|
      it "round-trips all encodings (compression=#{compressed}, max_zooms=#{zooms})" do
        io = IO::Memory.new
        options = BigWig::WriteOptions.new(compression: compressed, max_zooms: zooms, block_items: 2, index_block_size: 2)
        BigWig.create(io, {"chr" => 1000}, options: options) { |writer| write_mixed(writer) }
        io.closed?.should be_false
        BigWig.open(io) do |reader|
          reader.intervals("chr").should eq(mixed_intervals)
          reader.header.bases_covered.should eq(16_u64)
          reader.header.sum_data.should eq(21)
          reader.sum("chr", exact: true).should eq(21)
          if zooms == 0
            reader.zoom_levels.should be_empty
          else
            reader.zoom_levels.should_not be_empty
            reader.sum("chr").should eq(21)
          end
        end
      end
    end
  end

  it "accepts lazy input without consuming it up front" do
    io = IO::Memory.new
    consumed = 0
    input = (0...1000).each.map do |i|
      consumed += 1
      BigWig::Interval.new(i * 2, i * 2 + 1, i % 7)
    end
    options = BigWig::WriteOptions.new(block_items: 3, index_block_size: 2, max_zooms: 2)
    BigWig.create(io, {"chr" => 2000}, options: options) do |writer|
      consumed.should eq(0)
      writer.add("chr", input)
      consumed.should eq(1000)
    end
    BigWig.open(io) do |reader|
      reader.each_interval("chr").to_a.size.should eq(1000)
      reader.value("chr", 1998).should eq((999 % 7).to_f32)
      reader.sum("chr", exact: true).should eq((0...1000).sum { |i| i % 7 })
    end
  end

  it "writes multilevel chromosome trees with declaration-order IDs" do
    declarations = [{"z", 100_i64}, {"b", 100_i64}, {"a", 100_i64}, {"c", 100_i64}, {"huge", UInt32::MAX.to_i64}]
    io = IO::Memory.new
    options = BigWig::WriteOptions.new(block_items: 1, index_block_size: 2, max_zooms: 0)
    BigWig.create(io, declarations.each, options: options) do |writer|
      declarations.first(4).each_with_index { |pair, index| writer.add(pair[0], 1, 3, index) }
      writer.add("huge", 4_294_967_293_i64, UInt32::MAX, 7)
    end
    BigWig.open(io) do |reader|
      reader.chromosomes.keys.should eq(["z", "b", "a", "c", "huge"])
      reader.value("c", 1).should eq(3_f32)
      reader.value("huge", 4_294_967_294_i64).should eq(7_f32)
    end
  end

  it "writes valid empty files without zoom records" do
    io = IO::Memory.new
    BigWig.create(io, {"chr" => 100, "empty" => 0}) { |_| }
    BigWig.open(io) do |reader|
      reader.header.bases_covered.should eq(0)
      reader.header.zoom_levels.should eq(0)
      reader.intervals("chr").should be_empty
      reader.values("chr", 0, 3).all?(&.nan?).should be_true
      reader.mean("chr").nan?.should be_true
    end
  end

  it "distinguishes stored zeros from missing data in exact and zoom statistics" do
    io = IO::Memory.new
    BigWig.create(io, {"chr" => 1000}) do |writer|
      writer.add("chr", 0, 10, 0)
      writer.add("chr", 20, 30, -2)
      writer.add("chr", 30, 40, 2)
    end
    BigWig.open(io) do |reader|
      reader.value("chr", 0).should eq(0_f32)
      reader.value("chr", 10).nan?.should be_true
      reader.sum("chr").should eq(0)
      reader.sum("chr", exact: true).should eq(0)
      reader.coverage("chr").should be_close(0.03, 1e-12)
    end
  end

  it "omits zoom levels whose float32 summaries would overflow" do
    io = IO::Memory.new
    BigWig.create(io, {"chr" => 1000}) { |writer| writer.add("chr", 0, 2, 1e30) }
    BigWig.open(io) do |reader|
      reader.zoom_levels.should be_empty
      reader.header.sum_squared.finite?.should be_true
      reader.value("chr", 0).should eq(1e30_f32)
    end
  end

  it "rejects invalid declarations before creating output" do
    with_output_path do |path|
      expect_raises(ArgumentError) { BigWig.create(path, [] of Tuple(String, Int32)) }
      expect_raises(ArgumentError) { BigWig.create(path, [{"chr", 10}, {"chr", 20}]) }
      expect_raises(ArgumentError) { BigWig.create(path, {"chr" => -1}) }
      expect_raises(ArgumentError) { BigWig.create(path, {"a\0b" => 10}) }
      File.exists?(path).should be_false
    end
  end

  it "validates ordering and values before modifying the pending section" do
    io = IO::Memory.new
    BigWig.create(io, {"a" => 100, "b" => 100}) do |writer|
      writer.add("a", 10, 20, 1)
      expect_raises(ArgumentError) { writer.add("a", 0, 5, 2) }
      expect_raises(ArgumentError) { writer.add("a", 15, 25, 2) }
      expect_raises(ArgumentError) { writer.add("a", 20, 20, 2) }
      expect_raises(ArgumentError) { writer.add("a", 20, 101, 2) }
      expect_raises(ArgumentError) { writer.add("a", 20, 30, Float64::NAN) }
      expect_raises(ArgumentError) { writer.add("a", 20, 30, Float64::INFINITY) }
      expect_raises(ArgumentError) { writer.add("a", 20, 30, 1e100) }
      expect_raises(ArgumentError) { writer.add_span("a", start: 20, span: 2, step: 1, values: [1]) }
      expect_raises(ArgumentError) { writer.add_spans("a", [{20, 1}], span: 0) }
      writer.add("a", 20, 30, 2)
      writer.add("b", 0, 1, 3)
      expect_raises(ArgumentError) { writer.add("a", 30, 40, 4) }
    end
    BigWig.open(io, &.intervals("a").size.should(eq(2)))
  end

  it "does not publish a path until close succeeds" do
    with_output_path do |path|
      result = BigWig.create(path, {"chr" => 100}) do |writer|
        writer.add("chr", BigWig::Interval.new(1, 3, 2))
        writer.flush
        File.exists?(path).should be_false
        :done
      end
      result.should eq(:done)
      BigWig.open(path, &.value("chr", 1).should(eq(2_f32)))
    end
  end

  it "preserves existing paths unless overwrite is explicitly enabled" do
    with_output_path do |path|
      File.write(path, "original")
      expect_raises(File::AlreadyExistsError) { BigWig.create(path, {"chr" => 100}) }
      File.read(path).should eq("original")
      options = BigWig::WriteOptions.new(overwrite: true)
      expect_raises(Exception, "stop") do
        BigWig.create(path, {"chr" => 100}, options: options) do |writer|
          writer.add("chr", 0, 1, 1)
          raise "stop"
        end
      end
      File.read(path).should eq("original")
      BigWig.create(path, {"chr" => 100}, options: options) { |writer| writer.add("chr", 0, 1, 1) }
      BigWig.open(path, &.value("chr", 0).should(eq(1_f32)))
    end
  end

  it "aborts on block failure and cleans the staging file" do
    with_output_path do |path|
      expect_raises(Exception, "stop") do
        BigWig.create(path, {"chr" => 100}) do |writer|
          writer.add("chr", 0, 10, 1)
          raise "stop"
        end
      end
      File.exists?(path).should be_false
    end
  end

  it "handles explicit close and abort idempotently" do
    io = IO::Memory.new
    writer = BigWig.create(io, {"chr" => 100}, sync_close: true)
    writer.add("chr", 0, 1, 1)
    writer.abort
    writer.abort
    writer.close
    writer.aborted?.should be_true
    writer.closed?.should be_true
    io.closed?.should be_true
    expect_raises(BigWig::ClosedError) { writer.add("chr", 1, 2, 2) }
  end

  it "rejects nonempty borrowed IO without truncating it" do
    io = IO::Memory.new("existing".to_slice.dup)
    expect_raises(ArgumentError) { BigWig.create(io, {"chr" => 100}) }
    io.to_s.should eq("existing")
    io.closed?.should be_false
  end
  it "rejects write-only IO for zooms but accepts it with max_zooms zero" do
    with_output_path do |path|
      File.open(path, "wb") do |file|
        expect_raises(ArgumentError, "readable") { BigWig.create(file, {"chr" => 100}) }
        file.closed?.should be_false
        file.size.should eq(0)
      end
      File.open(path, "wb") do |file|
        options = BigWig::WriteOptions.new(max_zooms: 0)
        BigWig.create(file, {"chr" => 100}, options: options) { |writer| writer.add("chr", 0, 1, 2) }
        file.closed?.should be_false
      end
      BigWig.open(path, &.value("chr", 0).should(eq(2.0_f32)))
    end
  end

  it "preserves a destination created by another writer before publication" do
    with_output_path do |path|
      writer = BigWig.create(path, {"chr" => 100})
      begin
        writer.add("chr", 0, 1, 1)
        File.write(path, "other writer")
        expect_raises(File::AlreadyExistsError) { writer.close }
        File.read(path).should eq("other writer")
        writer.aborted?.should be_true
      ensure
        writer.abort
      end
    end
  end

  it "closes explicitly owned IO when chromosome validation fails" do
    io = IO::Memory.new
    expect_raises(ArgumentError) { BigWig.create(io, {"chr" => -1}, sync_close: true) }
    io.closed?.should be_true
  end
end

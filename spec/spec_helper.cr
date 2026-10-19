require "spec"
require "file_utils"
require "../src/bigwig"

def fixture(name : String) : String
  File.join(__DIR__, "fixtures", name)
end

def fixture_io(name : String) : IO::Memory
  IO::Memory.new(File.read(fixture(name)).to_slice.dup)
end

def mixed_intervals : Array(BigWig::Interval)
  [
    BigWig::Interval.new(0, 4, 1),
    BigWig::Interval.new(8, 10, 2),
    BigWig::Interval.new(20, 22, 3),
    BigWig::Interval.new(25, 27, -1),
    BigWig::Interval.new(40, 42, 0.5),
    BigWig::Interval.new(43, 45, 1.5),
    BigWig::Interval.new(46, 48, 2.5),
  ]
end

def write_mixed(writer : BigWig::Writer, chrom : String = "chr") : Nil
  writer.add(chrom, 0, 4, 1)
  writer.add(chrom, 8, 10, 2)
  writer.add_spans(chrom, [{20, 3.0}, {25, -1.0}].each, span: 2)
  writer.add_span(chrom, start: 40, span: 2, step: 3, values: [0.5, 1.5, 2.5].each)
end

def with_output_path(&)
  file = File.tempfile("bigwig-spec", ".bw")
  path = file.path
  file.close
  File.delete(path)
  begin
    yield path
  ensure
    File.delete(path) if File.exists?(path)
  end
end

def get_u64(io : IO::Memory, offset : Int) : UInt64
  IO::ByteFormat::LittleEndian.decode(UInt64, io.to_slice[offset, 8])
end

def patch(io : IO::Memory, offset : Int, value) : Nil
  io.seek(offset)
  io.write_bytes(value, IO::ByteFormat::LittleEndian)
  io.rewind
end

def same_values(actual : Enumerable(Float32), expected : Enumerable(Float32)) : Nil
  a, b = actual.to_a, expected.to_a
  a.size.should eq(b.size)
  a.zip(b).each do |x, y|
    if y.nan?
      x.nan?.should be_true
    else
      x.should eq(y)
    end
  end
end

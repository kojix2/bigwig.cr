require "option_parser"
require "../src/bigwig"

chrom : String? = nil
first = 0_i64
last : Int64? = nil
chunk_size = 65_536
omit_gaps = false

parser = OptionParser.new do |options|
  options.banner = "Usage: dump_values [options] FILE_OR_URL"
  options.on("-c CHROM", "--chrom=CHROM", "Chromosome to query (required)") { |value| chrom = value }
  options.on("-s START", "--start=START", "Zero-based start (default: 0)") do |value|
    first = value.to_i64? || raise ArgumentError.new("invalid start: #{value}")
  end
  options.on("-e END", "--end=END", "Exclusive end (default: chromosome length)") do |value|
    last = value.to_i64? || raise ArgumentError.new("invalid end: #{value}")
  end
  options.on("--chunk-size=N", "Bases read per chunk (default: 65536)") do |value|
    chunk_size = value.to_i? || raise ArgumentError.new("invalid chunk size: #{value}")
  end
  options.on("--omit-gaps", "Skip positions without data") { omit_gaps = true }
  options.on("-h", "--help", "Show this help") { puts options; exit 0 }
end

begin
  args = ARGV.dup
  parser.parse(args)
  unless args.size == 1 && chrom && first >= 0 && chunk_size > 0 && (!last || last.not_nil! >= first)
    STDERR.puts parser
    exit 2
  end
  name = chrom.not_nil!
  BigWig.open(args[0]) do |reader|
    raise ArgumentError.new("input must be a bigWig file") unless reader.bigwig?
    size = reader.chromosome_size(name) || raise ArgumentError.new("unknown chromosome: #{name}")
    stop = (last || size).not_nil!
    raise ArgumentError.new("region exceeds chromosome length #{size}") if stop > size

    buffer = Slice(Float32).new(chunk_size)
    position = first
    while position < stop
      length = Math.min(buffer.size.to_i64, stop - position).to_i
      part = buffer[0, length]
      reader.read_values(name, position, part)
      part.each_with_index do |value, index|
        next if omit_gaps && value.nan?
        STDOUT << name << '\t' << (position + index) << '\t' << value << '\n'
      end
      position += length
    end
  end
rescue ex : OptionParser::Exception | ArgumentError
  STDERR.puts "dump_values: #{ex.message}"
  exit 2
rescue ex
  STDERR.puts "dump_values: #{ex.message}"
  exit 1
end

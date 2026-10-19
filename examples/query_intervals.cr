require "option_parser"
require "../src/bigwig"

chrom : String? = nil
first = 0_i64
last : Int64? = nil
list_chromosomes = false

parser = OptionParser.new do |options|
  options.banner = "Usage: query_intervals [options] FILE_OR_URL"
  options.on("-c CHROM", "--chrom=CHROM", "Chromosome to query") { |value| chrom = value }
  options.on("-s START", "--start=START", "Zero-based start (default: 0)") do |value|
    first = value.to_i64? || raise ArgumentError.new("invalid start: #{value}")
  end
  options.on("-e END", "--end=END", "Exclusive end (default: chromosome length)") do |value|
    last = value.to_i64? || raise ArgumentError.new("invalid end: #{value}")
  end
  options.on("-l", "--list-chromosomes", "Print chromosome names and sizes") { list_chromosomes = true }
  options.on("-h", "--help", "Show this help") { puts options; exit 0 }
end

begin
  args = ARGV.dup
  parser.parse(args)
  unless args.size == 1 && (list_chromosomes || chrom) && first >= 0 && (!last || last.not_nil! >= first)
    STDERR.puts parser
    exit 2
  end
  if list_chromosomes && (chrom || first != 0 || last)
    STDERR.puts "--list-chromosomes cannot be combined with region options"
    exit 2
  end

  BigWig.open(args[0]) do |reader|
    if list_chromosomes
      reader.chromosomes.each { |name, size| STDOUT << name << '\t' << size << '\n' }
    elsif name = chrom
      if reader.bigwig?
        reader.each_interval(name, first, last) do |record|
          STDOUT << name << '\t' << record.start << '\t' << record.stop << '\t' << record.value << '\n'
        end
      else
        reader.each_entry(name, first, last) do |record|
          STDOUT << name << '\t' << record.start << '\t' << record.stop
          STDOUT << '\t' << record.rest unless record.rest.empty?
          STDOUT << '\n'
        end
      end
    end
  end
rescue ex : OptionParser::Exception | ArgumentError
  STDERR.puts "query_intervals: #{ex.message}"
  exit 2
rescue ex
  STDERR.puts "query_intervals: #{ex.message}"
  exit 1
end

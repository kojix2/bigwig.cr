require "option_parser"
require "../src/bigwig"

overwrite = false
compression = true
max_zooms = 10

parser = OptionParser.new do |options|
  options.banner = "Usage: bedgraph_to_bigwig [options] CHROM_SIZES BEDGRAPH OUTPUT_BW"
  options.on("--overwrite", "Replace an existing output file") { overwrite = true }
  options.on("--no-compression", "Write uncompressed data blocks") { compression = false }
  options.on("--no-zoom", "Omit zoom levels") { max_zooms = 0 }
  options.on("-h", "--help", "Show this help") { puts options; exit 0 }
end

begin
  args = ARGV.dup
  parser.parse(args)
  unless args.size == 3 && args[0] != "-" && args[2] != "-"
    STDERR.puts parser
    exit 2
  end
rescue ex : OptionParser::Exception
  STDERR.puts "bedgraph_to_bigwig: #{ex.message}"
  exit 2
end

sizes_path, bedgraph_path, output_path = args

begin
  chromosomes = Array(Tuple(String, Int64)).new
  File.open(sizes_path) do |input|
    input.each_line.with_index(1) do |line, number|
      line = line.strip
      next if line.empty? || line.starts_with?('#')
      fields = line.split
      unless fields.size == 2 && (size = fields[1].to_i64?) && size > 0
        raise ArgumentError.new("#{sizes_path}:#{number}: expected CHROM and positive SIZE")
      end
      chromosomes << {fields[0], size}
    end
  end

  write_options = BigWig::WriteOptions.new(
    overwrite: overwrite, compression: compression, max_zooms: max_zooms
  )
  BigWig.create(output_path, chromosomes, options: write_options) do |writer|
    input = bedgraph_path == "-" ? STDIN : File.open(bedgraph_path)
    begin
      input.each_line.with_index(1) do |line, number|
        line = line.strip
        next if line.empty? || line.starts_with?('#') ||
                line.starts_with?("track ") || line.starts_with?("browser ")
        fields = line.split
        unless fields.size == 4 && (start = fields[1].to_i64?) &&
               (stop = fields[2].to_i64?) && (value = fields[3].to_f64?)
          raise ArgumentError.new("#{bedgraph_path}:#{number}: expected CHROM START END VALUE")
        end
        begin
          writer.add(fields[0], start, stop, value)
        rescue ex : ArgumentError
          raise ArgumentError.new("#{bedgraph_path}:#{number}: #{ex.message}")
        end
      end
    ensure
      input.close unless input == STDIN
    end
  end
rescue ex
  STDERR.puts "bedgraph_to_bigwig: #{ex.message}"
  exit 1
end

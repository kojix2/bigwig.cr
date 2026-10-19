# Produce actual Crystal-written files for the optional external C reader check.
require "../../src/bigwig"

unless ARGV.size == 1
  STDERR.puts "Usage: crystal run scripts/reference/generate.cr -- OUTPUT_DIRECTORY"
  exit 2
end

begin
  directory = ARGV[0]
  Dir.mkdir_p(directory)
  {
    {"crystal-raw.bw", false, 0},
    {"crystal-zlib.bw", true, 0},
    {"crystal-zoom.bw", true, 4},
  }.each do |name, compressed, zooms|
    options = BigWig::WriteOptions.new(
      compression: compressed, max_zooms: zooms,
      block_items: 2, index_block_size: 2, overwrite: true
    )
    BigWig.create(File.join(directory, name), {"chr" => 1_000}, options: options) do |writer|
      writer.add("chr", 0, 4, 1)
      writer.add("chr", 8, 10, 2)
      writer.add_spans("chr", [{20, 3.0}, {25, -1.0}].each, span: 2)
      writer.add_span("chr", start: 40, span: 2, step: 3, values: [0.5, 1.5, 2.5].each)
    end
    puts File.join(directory, name)
  end
rescue ex
  STDERR.puts "#{ex.class}: #{ex.message}"
  exit 1
end

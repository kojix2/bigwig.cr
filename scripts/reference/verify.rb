#!/usr/bin/env ruby
# Validate fixtures and Crystal-written files through bio-bigwig (libBigWig).
require "json"
require "optparse"
require "fileutils"

begin
  require "bio/bigwig"
rescue LoadError => e
  warn "bio-bigwig is required: gem install bio-bigwig (#{e.message})"
  exit 2
end

ROOT = File.expand_path("../..", __dir__)
KINDS = %i[mean std max min coverage sum].freeze

def check(condition, message)
  raise message unless condition
end

def close_enough?(actual, expected)
  (actual.nan? && expected.nan?) ||
    (actual - expected).abs <= [1e-7, 2e-6 * [actual.abs, expected.abs].max].max
end

def exact_stat(rows, start, stop, kind)
  hits = rows.filter_map do |left, right, value|
    width = [right, stop].min - [left, start].max
    [width, value] if width.positive?
  end
  covered = hits.sum { |width, _| width }
  return Float::NAN if covered.zero?

  total = hits.sum { |width, value| width * value }
  mean = total / covered
  case kind
  when :sum then total
  when :mean then mean
  when :coverage then covered.to_f / (stop - start)
  when :min then hits.map(&:last).min
  when :max then hits.map(&:last).max
  when :std
    covered > 1 ? Math.sqrt(hits.sum { |width, value| width * (value - mean)**2 } / (covered - 1)) : 0.0
  end
end

def check_bigwig(path, chromosomes)
  checks = 0
  Bio::BigWig.open(path) do |bw|
    check(bw.is_bigwig?, "#{path}: expected bigWig")
    random = Random.new(1234)
    chromosomes.each do |chrom, (size, rows)|
      check(bw.chroms(chrom) == size, "#{path}: incorrect size for #{chrom}")
      actual = bw.intervals(chrom, 0, size)
      check(actual == rows, "#{path}: incorrect intervals for #{chrom}: #{actual.inspect}")
      checks += 2

      limit = [size, 200].min
      queries = [[0, limit]]
      queries.concat(rows.map { |start, stop, _| [[0, start - 2].max, [size, stop + 2].min] })
      30.times do
        start = random.rand(limit)
        queries << [start, random.rand((start + 1)..limit)]
      end
      queries.each do |start, stop|
        KINDS.each do |kind|
          bins = [5, stop - start].min
          values = bw.stats(chrom, start, stop, nbins: bins, type: kind, exact: true)
          check(values.size == bins, "#{path}: #{chrom}: incorrect #{kind} bin count")
          values.each_with_index do |value, index|
            low = start + (stop - start) * index / bins
            high = start + (stop - start) * (index + 1) / bins
            expected = exact_stat(rows, low, high, kind)
            check(close_enough?(value, expected),
                  "#{path}: #{chrom}: #{low}-#{high} #{kind}: #{value.inspect} != #{expected.inspect}")
            checks += 1
          end
        end
      end
    end
  end
  checks
end

def check_bigbed(path, chrom, start, stop, expected)
  Bio::BigWig.open(path) do |bb|
    check(bb.is_bigbed?, "#{path}: expected bigBed")
    actual = bb.entries(chrom, start, stop)
    check(actual == expected, "#{path}: incorrect entries: #{actual.inspect}")
    check(bb.sql&.include?("table"), "#{path}: missing AutoSQL")
  end
  2
end

def f32(value)
  [value].pack("f").unpack1("f")
end

options = {}
parser = OptionParser.new do |opts|
  opts.banner = "Usage: ruby scripts/reference/verify.rb [--interop-dir DIR] [--report FILE]"
  opts.on("--interop-dir DIR", "Also check Crystal-written files in DIR") { |dir| options[:interop_dir] = dir }
  opts.on("--report FILE", "Write the JSON report to FILE") { |file| options[:report] = file }
  opts.on("-h", "--help", "Show this help") { puts opts; exit 0 }
end

begin
  parser.parse!
  check(ARGV.empty?, "unexpected arguments: #{ARGV.join(" ")}")
  fixtures = File.join(ROOT, "spec/fixtures")
  old = {
    "1" => [195_471_971, [[0, 1, f32(0.1)], [1, 2, f32(0.2)], [2, 3, f32(0.3)],
                           [100, 150, f32(1.4)], [150, 151, f32(1.5)]]],
    "10" => [130_694_993, [[200, 300, f32(2.0)]]]
  }
  mixed = {
    "zchr" => [1_000, [[0, 4, 1.0], [8, 10, 2.0], [20, 22, 3.0], [25, 27, -1.0],
                       [40, 42, 0.5], [43, 45, 1.5], [46, 48, 2.5]]],
    "achr" => [2**32 - 1, [[3_000_000_000, 3_000_000_003, 4.0], [2**32 - 3, 2**32 - 1, 5.0]]],
    "empty" => [10, []]
  }
  report = {
    kind: "bio-bigwig fixture/reference validation",
    crystal_compiler_executed: false,
    checks: 0,
    files: []
  }
  {
    "test.bw" => old,
    "mixed-little-raw.bw" => mixed,
    "mixed-little-zlib.bw" => mixed,
    "zoom-positive.bw" => { "chr" => [128, [[0, 16, 2.0], [32, 48, 4.0]]] }
  }.each do |name, expected|
    report[:checks] += check_bigwig(File.join(fixtures, name), expected)
    report[:files] << name
  end
  report[:checks] += check_bigbed(File.join(fixtures, "test.bigBed"), "chr1",
                                10_000_000, 10_020_000,
                                [[10_009_333, 10_009_640, "61035\t130\t-\t0.026\t0.42\t404"],
                                 [10_014_007, 10_014_289, "61047\t136\t-\t0.029\t0.42\t404"],
                                 [10_014_373, 10_024_307, "61048\t630\t-\t5.420\t0.00\t2672399"]])
  report[:checks] += check_bigbed(File.join(fixtures, "minimal-little.bb"), "chr", 0, 10,
                                [[1, 4, "geneA\t10\t+"], [4, 4, "insertion\t0\t."],
                                 [5, 9, "geneB\t20\t-"]])
  report[:files].concat(%w[test.bigBed minimal-little.bb])

  Bio::BigWig.open(File.join(fixtures, "zoom-positive.bw")) do |bw|
    raw = bw.stats("chr", 0, 128, type: :sum, exact: true).first
    zoom = bw.stats("chr", 0, 128, type: :sum, exact: false).first
    check(raw == 96.0 && zoom == 6.0, "unexpected zoom sum: exact=#{raw}, zoom=#{zoom}")
    report[:confirmed_reference_bug] = {
      file: "zoom-positive.bw", exact_sum: raw, zoom_sum: zoom, crystal_contract_sum: 96.0
    }
    report[:checks] += 1
  end

  report[:crystal_output_files_checked] = []
  if dir = options[:interop_dir]
    expected = { "chr" => mixed.fetch("zchr") }
    %w[crystal-raw.bw crystal-zlib.bw crystal-zoom.bw].each do |name|
      report[:checks] += check_bigwig(File.join(dir, name), expected)
      report[:crystal_output_files_checked] << name
    end
    report[:note] = "Crystal-written files were checked by bio-bigwig; this script does not run Crystal."
  else
    report[:note] = "No Crystal code was run; these checks validate fixtures and numerical expectations only."
  end
  report[:status] = "passed"
  output = JSON.pretty_generate(report) + "\n"
  if path = options[:report]
    FileUtils.mkdir_p(File.dirname(path))
    File.write(path, output)
  end
  puts output
rescue OptionParser::ParseError => e
  warn "verify: #{e.message}"
  exit 2
rescue StandardError => e
  warn "verify: #{e.message}"
  exit 1
end

# bigwig.cr

[![test](https://github.com/kojix2/bigwig.cr/actions/workflows/test.yml/badge.svg)](https://github.com/kojix2/bigwig.cr/actions/workflows/test.yml)
[![Lines of Code](https://img.shields.io/endpoint?url=https%3A%2F%2Ftokei.kojix2.net%2Fbadge%2Fgithub%2Fkojix2%2Fbigwig.cr%2Flines)](https://tokei.kojix2.net/github/kojix2/bigwig.cr)
[![DOI](https://zenodo.org/badge/1425265010.svg)](https://doi.org/10.5281/zenodo.23470329)

A native Crystal implementation based on [libBigWig](https://github.com/dpryan79/libBigWig). 

The code was generated with LLM coding agent.

## Install

Add to `shard.yml`:

```yaml
dependencies:
  bigwig:
    github: kojix2/bigwig
```

## Read

```crystal
require "bigwig"

BigWig.open("signal.bw") do |bw|
  bw.each_interval("chr1", 1_000...2_000) do |interval|
    puts "#{interval.start}\t#{interval.stop}\t#{interval.value}"
  end

  values = bw.values("chr1", 1_000, 2_000)
  means = bw.stats("chr1", 0, 1_000_000, type: :mean, bins: 100)
end
```

Coordinates are zero-based and half-open. `values` returns one `Float32` per base, with `NaN` for gaps. `stats` returns `Array(Float64)`; it may use zoom data for approximate results. Pass `exact: true` for exact results.

Use `each_entry` for bigBed:

```crystal
BigWig.open("annotations.bb") do |bb|
  bb.each_entry("chr1", 1_000, 2_000) do |entry|
    puts "#{entry.start}\t#{entry.stop}\t#{entry.rest}"
  end
end
```

## Write

```crystal
chromosomes = {"chr1" => 248_956_422_i64}

BigWig.create("output.bw", chromosomes) do |bw|
  bw.add("chr1", 100, 200, 1.5)
  bw.add("chr1", 200, 300, 2.0)
end
```

Add intervals in chromosome declaration order, sorted and non-overlapping within each chromosome. Existing output files are preserved by default. Pass `options: BigWig::WriteOptions.new(overwrite: true)` to replace one. bigBed writing is not supported.

## Input and output

`BigWig.open` accepts a local path, an HTTP(S) URL, or seekable `IO`. `BigWig.create` accepts a local path or `IO`. Caller-provided `IO` stays open unless `sync_close: true` is passed.

```crystal
BigWig.open("https://example.org/signal.bw") do |bw|
  puts bw.max("chr1", 0, 100_000)
end
```

HTTP(S) reading requires a server that supports Range requests; there is no full-download fallback. It uses Crystal's standard `HTTP::Client`, which does not use `HTTP_PROXY` or `HTTPS_PROXY` automatically.

## Command-line examples

The files in `examples/` are standalone CLI tools. Run any of them with `--help` for all options.

```sh
crystal run examples/query_intervals.cr -- --list-chromosomes signal.bw
crystal run examples/query_intervals.cr -- --chrom chr1 --start 1000 --end 2000 signal.bw
crystal run examples/dump_values.cr -- --chrom chr1 --start 1000 --end 2000 --omit-gaps signal.bw
crystal run examples/bedgraph_to_bigwig.cr -- chrom.sizes input.bedgraph output.bw
```

`query_intervals` prints overlapping bigWig intervals or bigBed entries as tab-separated rows. `dump_values` prints chromosome, zero-based position, and value; it includes `NaN` gaps unless `--omit-gaps` is set. The converter accepts `-` instead of a bedGraph path to read standard input. Its input must be sorted by chromosome declaration order and contain no overlapping intervals.

## Development

```sh
gem install bio-bigwig
ruby scripts/check.rb
```

`scripts/check.rb` is also the CI entry point. It runs specs, checks fixture bytes, compiles examples, writes three test files with this library, and verifies them with bio-bigwig. The reference report is saved to `build/reference-validation.json`. Regenerate fixtures with `crystal run scripts/fixtures/generate.cr -- --write`.

## License

[MIT](LICENSE).

require "set"
require "compress/zlib"
require "http/client"

require "./bigwig/version"
require "./bigwig/errors"
require "./bigwig/types"
require "./bigwig/binary"
require "./bigwig/source"
require "./bigwig/http"
require "./bigwig/index"
require "./bigwig/iterators"
require "./bigwig/reader"
require "./bigwig/statistics"
require "./bigwig/write_index"
require "./bigwig/zoom_writer"
require "./bigwig/writer"

# Native Crystal bigWig/bigBed I/O, independent of libBigWig.
module BigWig
  # The library owns files and HTTP connections opened from a path/URL.
  def self.open(path : String | Path, *, options : ReadOptions = ReadOptions.new,
                headers : HTTP::Headers = HTTP::Headers.new) : Reader
    options.validate
    text = path.to_s
    if path.is_a?(String) && (text.starts_with?("http://") || text.starts_with?("https://"))
      source = Internal::HTTPSource.new(text, options, headers)
    else
      file = ::File.open(path, "rb")
      begin
        source = Internal::FileSource.new(file)
      rescue ex
        file.close
        raise ex
      end
    end
    Reader.new(source, options)
  end

  # Only seekable IO is accepted. It is not silently copied into memory.
  # Caller-provided IO, including File, is read under a lock with its position restored.
  def self.open(io : IO, *, sync_close : Bool = false,
                options : ReadOptions = ReadOptions.new) : Reader
    options.validate
    begin
      source = Internal::IOSource.new(io, sync_close)
      Reader.new(source, options)
    rescue ex
      io.close if sync_close && !io.closed?
      raise ex
    end
  end

  # Like File.open, the block overload returns the value of the block.
  def self.open(input, **options, &)
    reader = open(input, **options)
    begin
      yield reader
    ensure
      reader.close
    end
  end

  # Path outputs are staged in the destination directory and published only
  # after successful finalization. Existing outputs are preserved by default.
  def self.create(path : String | Path, chromosomes, *, options : WriteOptions = WriteOptions.new) : Writer
    Writer.to_path(path, chromosomes, options)
  end

  def self.create(io : IO, chromosomes, *, sync_close : Bool = false,
                  options : WriteOptions = WriteOptions.new) : Writer
    Writer.new(io, chromosomes, sync_close: sync_close, options: options)
  rescue ex
    io.close if sync_close && !io.closed?
    raise ex
  end

  def self.create(output, chromosomes, **options, &)
    writer = create(output, chromosomes, **options)
    begin
      value = yield writer
      writer.close
      value
    rescue ex
      writer.abort
      raise ex
    ensure
      writer.abort unless writer.closed?
    end
  end
end

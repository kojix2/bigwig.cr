module BigWig::Internal
  # :nodoc:
  abstract class Source
    abstract def size : Int64
    abstract def read_exact_at(offset : Int64, buffer : Bytes) : Nil
    abstract def close : Nil
    abstract def closed? : Bool

    def check_range(offset : Int64, length : Int) : Nil
      raise BigWig::ClosedError.new("source is closed") if closed?
      unless offset >= 0 && length >= 0 && offset <= size && length.to_i64 <= size - offset
        raise BigWig::FormatError.new("read outside file: offset=#{offset}, length=#{length}, size=#{size}")
      end
    end

    def read(offset : Int64, length : Int, limit : Int32) : Bytes
      check_range(offset, length)
      raise BigWig::LimitError.new("block exceeds max_block_bytes (#{limit})") if length > limit
      bytes = Bytes.new(length.to_i)
      read_exact_at(offset, bytes)
      bytes
    end
  end

  # IO ownership is explicit; reads preserve the wrapped stream's position.
  # A mutex makes interleaved cursors independent, including over IO::Memory.
  class IOSource < Source
    getter size : Int64
    @closed = false
    @mutex = Mutex.new

    def initialize(@io : IO, @sync_close : Bool = false)
      original = @io.pos
      begin
        @io.seek(0, IO::Seek::End)
        @size = @io.pos.to_i64
      ensure
        @io.seek(original)
      end
    end

    def closed? : Bool
      @closed || @io.closed?
    end

    def read_exact_at(offset : Int64, buffer : Bytes) : Nil
      @mutex.synchronize do
        check_range(offset, buffer.size)
        previous = @io.pos
        begin
          @io.seek(offset)
          @io.read_fully(buffer)
        rescue IO::EOFError
          raise BigWig::FormatError.new("unexpected EOF at byte #{offset}")
        ensure
          @io.seek(previous) unless @io.closed?
        end
      end
    end

    def close : Nil
      @mutex.synchronize do
        return if @closed
        @closed = true
        @io.close if @sync_close
      end
    end
  end

  # Owned files have no caller-visible cursor to preserve, so absolute reads
  # can use File#read_at directly.
  class FileSource < Source
    getter size : Int64
    @closed = false
    @mutex = Mutex.new

    def initialize(@file : ::File)
      @size = @file.size.to_i64
    end

    def closed? : Bool
      @closed || @file.closed?
    end

    def read_exact_at(offset : Int64, buffer : Bytes) : Nil
      @mutex.synchronize do
        check_range(offset, buffer.size)
        return if buffer.empty?
        @file.read_at(offset, buffer.size, &.read_fully(buffer))
      end
    rescue IO::EOFError
      raise BigWig::FormatError.new("unexpected EOF at byte #{offset}")
    rescue ex : ArgumentError
      # File#read_at checks the current filesystem size, whereas this source
      # validates against its opening size. A mismatch indicates truncation.
      raise BigWig::FormatError.new("file changed during positioned read: #{ex.message}")
    end

    def close : Nil
      @mutex.synchronize do
        return if @closed
        @closed = true
        @file.close
      end
    end
  end
end

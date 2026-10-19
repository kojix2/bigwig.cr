module BigWig
  # Streaming bigWig writer. Input must follow chromosome declaration order
  # and have increasing, non-overlapping coordinates. close finalizes indexes
  # and zoom levels; abort discards a path-backed temporary output.
  class Writer
    @chroms : Array(Internal::Chromosome)
    @lookup : Hash(String, Internal::Chromosome)
    @index_spool : Internal::TempFile
    @chrom_offset = 0_i64
    @data_offset = 0_i64
    @index_offset = 0_i64
    @summary_offset = 0_i64
    @section_count = 0_i64
    @entry_count = 0_i64
    @total_width = 0_u64
    @aggregate = Internal::Aggregate.new
    @section = IO::Memory.new
    @items = 0
    @kind = 0_u8
    @tid = 0_u32
    @first = 0_i64
    @last = 0_i64
    @span = 0_i64
    @step = 0_i64
    @next_start = 0_i64
    @previous_tid : UInt32? = nil
    @previous_stop = 0_i64
    @closed = false
    @aborted = false
    @temporary_path : String? = nil

    # For existing, empty, seekable IO. The IO remains open by default.
    # destination is used by the path factory, not by normal IO callers.
    def initialize(@io : IO, chromosomes, *, @options : WriteOptions = WriteOptions.new,
                   @sync_close : Bool = false, @destination : String? = nil)
      @options.validate
      @chroms = self.class.normalize_chromosomes(chromosomes)
      @lookup = Hash(String, Internal::Chromosome).new
      @chroms.each { |chrom| @lookup[chrom.name] = chrom }
      @index_spool = Internal::TempFile.new
      if @destination && @io.is_a?(::File)
        @temporary_path = @io.as(::File).path
      end
      begin
        unless @io.pos == 0
          raise ArgumentError.new("output IO must start at position zero")
        end
        @io.seek(0, IO::Seek::End)
        unless @io.pos == 0
          raise ArgumentError.new("output IO must be empty; existing IO is not truncated implicitly")
        end
        @io.seek(0)
        if @options.max_zooms > 0
          begin
            @io.read_byte
            @io.rewind
          rescue IO::Error
            raise ArgumentError.new("zoom generation requires readable output IO; open read/write or set max_zooms: 0")
          end
        end
        Internal.zeros(@io, 64_i64 + @options.max_zooms.to_i64 * 24)
        @summary_offset = @io.pos.to_i64
        Internal.zeros(@io, 40_i64)
        @chrom_offset = Internal::ChromTreeWriter.write(@io, @chroms, @options.index_block_size)
        @data_offset = @io.pos.to_i64
        Internal.put(@io, 0_u64)
      rescue ex
        abort_after_failure
        raise ex
      end
    end

    # :nodoc:
    def self.normalize_chromosomes(chromosomes) : Array(Internal::Chromosome)
      output = Array(Internal::Chromosome).new
      seen = Set(String).new
      chromosomes.each do |name, length|
        raise ArgumentError.new("chromosome name cannot be empty or contain NUL") if name.empty? || name.includes?('\0')
        raise ArgumentError.new("duplicate chromosome: #{name}") unless seen.add?(name)
        size = Internal.coordinate(length, "chromosome length")
        output << Internal::Chromosome.new(name, output.size.to_u32, size)
      end
      raise ArgumentError.new("at least one chromosome is required") if output.empty?
      output
    end

    # :nodoc:
    def self.to_path(path : String | Path, chromosomes, options : WriteOptions) : Writer
      options.validate
      target = ::File.expand_path(path.to_s)
      if !options.overwrite && (::File.exists?(target) || ::File.symlink?(target))
        raise ::File::AlreadyExistsError.new("output already exists", file: target)
      end
      # Validate before creating temporary files. This also supports a one-shot
      # enumerable of chromosome declarations without consuming it twice.
      normalized = normalize_chromosomes(chromosomes)
      declarations = normalized.map { |chrom| {chrom.name, chrom.size} }
      io = ::File.tempfile(".bigwig", ".tmp", dir: ::File.dirname(target))
      begin
        new(io, declarations, options: options, sync_close: true, destination: target)
      rescue ex
        io.close unless io.closed?
        ::File.delete(io.path) if ::File.exists?(io.path)
        raise ex
      end
    end

    def closed? : Bool
      @closed
    end

    def aborted? : Bool
      @aborted
    end

    private def check_open : Nil
      raise ClosedError.new("writer is closed") if @closed
    end

    private def abort_after_failure : Nil
      abort
    end

    private def float32(value : Int | Float) : Float32
      value.to_f32
    rescue OverflowError
      raise ArgumentError.new("value must be representable as Float32")
    end

    # Add a bedGraph-style interval, using the file's float32 value precision.
    def add(chrom : String, start : Int, stop : Int, value : Int | Float) : Nil
      chromosome = resolve_chromosome(chrom)
      first = Internal.coordinate(start, "start")
      last = Internal.coordinate(stop, "stop")
      append(chromosome, first, last, float32(value), 1_u8, 0_i64, 0_i64)
    end

    def add(chrom : String, interval : Interval) : Nil
      add(chrom, interval.start, interval.stop, interval.value)
    end

    # Accepts a lazy Enumerable; it is never converted to an array.
    def add(chrom : String, intervals : Enumerable(Interval)) : Nil
      intervals.each { |interval| add(chrom, interval) }
    end

    # fixedStep encoding. Coordinates remain zero-based (unlike WIG text).
    def add_span(chrom : String, *, start : Int, span : Int, values,
                 step : Int = span) : Nil
      chromosome = resolve_chromosome(chrom)
      position = Internal.coordinate(start, "start")
      span64 = Internal.coordinate(span, "span")
      step64 = Internal.coordinate(step, "step")
      unless span64 > 0 && step64 >= span64 && position <= chromosome.size
        raise ArgumentError.new("require span > 0, step >= span, and start within the chromosome")
      end
      values.each do |value|
        append(chromosome, position, position + span64, float32(value), 3_u8, span64, step64)
        position += step64
      end
    end

    # variableStep encoding. entries lazily yields {start, value} tuples.
    def add_spans(chrom : String, entries, *, span : Int) : Nil
      chromosome = resolve_chromosome(chrom)
      span64 = Internal.coordinate(span, "span")
      raise ArgumentError.new("span must be positive") if span64 == 0
      entries.each do |entry|
        first = Internal.coordinate(entry[0], "start")
        append(chromosome, first, first + span64, float32(entry[1]), 2_u8, span64, 0_i64)
      end
    end

    private def resolve_chromosome(chrom : String) : Internal::Chromosome
      check_open
      @lookup[chrom]? || raise ArgumentError.new("unknown chromosome: #{chrom}")
    end

    private def validate_append(chrom : Internal::Chromosome, first : Int64,
                                last : Int64, value : Float32) : Nil
      unless 0 <= first < last <= chrom.size && value.finite?
        raise ArgumentError.new("require a nonempty interval inside the chromosome and a finite float32 value")
      end
      if (previous = @previous_tid) && (chrom.id < previous || (chrom.id == previous && first < @previous_stop))
        raise ArgumentError.new("write in chromosome declaration order with sorted, non-overlapping intervals")
      end
    end

    private def prepare_section(chrom : Internal::Chromosome, first : Int64,
                                kind : UInt8, span : Int64, step : Int64) : Nil
      if @items > 0 && (@tid != chrom.id || @kind != kind || @span != span || @step != step ||
         @items == @options.block_items || (kind == 3 && first != @next_start))
        flush_section
      end
      return unless @items == 0
      @section.clear
      Internal.zeros(@section, 24_i64)
      @tid, @first, @kind, @span, @step = chrom.id, first, kind, span, step
    end

    private def write_interval_start(first : Int64, last : Int64, kind : UInt8) : Nil
      case kind
      when 1
        Internal.put(@section, first.to_u32)
        Internal.put(@section, last.to_u32)
      when 2
        Internal.put(@section, first.to_u32)
      end
    end

    private def append(chrom : Internal::Chromosome, first : Int64, last : Int64,
                       value : Float32, kind : UInt8, span : Int64, step : Int64) : Nil
      check_open
      validate_append(chrom, first, last, value)
      begin
        prepare_section(chrom, first, kind, span, step)
        write_interval_start(first, last, kind)
        Internal.put(@section, value)
        @items += 1
        @last = last
        @next_start = first + step
        @previous_tid, @previous_stop = chrom.id, last
        @entry_count += 1
        @total_width += (last - first).to_u64
        @aggregate.add(value.to_f64, (last - first).to_f64)
      rescue ex
        abort_after_failure
        raise ex
      end
    end

    private def flush_section : Nil
      return if @items == 0
      @section.rewind
      Internal.put(@section, @tid)
      Internal.put(@section, @first.to_u32)
      Internal.put(@section, @last.to_u32)
      Internal.put(@section, @step.to_u32)
      Internal.put(@section, @span.to_u32)
      Internal.put(@section, @kind)
      Internal.put(@section, 0_u8)
      Internal.put(@section, @items.to_u16)
      offset, size = Internal.write_data_block(@io, @section.to_slice, @options)
      Internal::BlockRef.new(@tid, @first, @tid, @last, offset, size).write(@index_spool.io)
      @section_count += 1
      @items = 0
      @section.clear
    end

    # Flush pending data. This does not finalize the file or publish a path.
    def flush : Nil
      check_open
      begin
        flush_section
        @io.flush
      rescue ex
        abort_after_failure
        raise ex
      end
    end

    private def write_header(zooms : Array(Internal::ZoomLevel)) : Nil
      previous = @io.pos
      begin
        @io.rewind
        Internal.put(@io, Internal::BW_MAGIC)
        Internal.put(@io, 4_u16)
        Internal.put(@io, zooms.size.to_u16)
        Internal.put(@io, @chrom_offset.to_u64)
        Internal.put(@io, @data_offset.to_u64)
        Internal.put(@io, @index_offset.to_u64)
        Internal.put(@io, 0_u16)
        Internal.put(@io, 0_u16)
        Internal.put(@io, 0_u64)
        Internal.put(@io, @summary_offset.to_u64)
        maximum_block = Math.max(24 + @options.block_items * 12, @options.block_items * 32)
        Internal.put(@io, @options.compression ? maximum_block.to_u32 : 0_u32)
        Internal.put(@io, 0_u64)
        zooms.each do |zoom|
          Internal.put(@io, zoom.reduction.to_u32)
          Internal.put(@io, 0_u32)
          Internal.put(@io, zoom.data_offset.to_u64)
          Internal.put(@io, zoom.index_offset.to_u64)
        end
        @io.seek(@summary_offset)
        Internal.put(@io, @total_width)
        Internal.put(@io, @entry_count == 0 ? 0.0 : @aggregate.min)
        Internal.put(@io, @entry_count == 0 ? 0.0 : @aggregate.max)
        Internal.put(@io, @aggregate.sum)
        Internal.put(@io, @aggregate.sum_squared)
        @io.seek(@data_offset)
        Internal.put(@io, @section_count.to_u64)
      ensure
        @io.seek(previous)
      end
    end

    private def build_zooms : Array(Internal::ZoomLevel)
      output = Array(Internal::ZoomLevel).new
      return output if @options.max_zooms == 0 || @entry_count == 0
      max_size = @chroms.max_of(&.size)
      return output if max_size == 0
      reduction = Math.max(10_u64, (@total_width // @entry_count.to_u64) * 16)
      reduction = Math.min(reduction, max_size.to_u64)
      builders = Array(Internal::ZoomBuilder).new
      begin
        @options.max_zooms.times do
          break if reduction > max_size.to_u64 || reduction > UInt32::MAX
          builders << Internal::ZoomBuilder.new(reduction.to_i64)
          reduction *= 4
        end
        # The primary index and header are already valid at this point.
        # File/IO sources preserve the writer's current position.
        BigWig.open(@io) do |reader|
          @chroms.each do |chrom|
            reader.each_interval(chrom.name) do |interval|
              builders.each(&.add(chrom.id, chrom.size, interval))
            end
          end
        end
        builders.each do |builder|
          builder.finish
          if builder.valid? && builder.count > 0
            output << Internal::ZoomWriter.append(@io, builder, @options)
          end
        end
        output
      ensure
        builders.each(&.close)
      end
    end

    # Finalize and publish. Explicit close inside a block commits immediately;
    # a later exception in that user block cannot undo an already committed file.
    def close : Nil
      return if @closed
      begin
        flush_section
        @index_offset = Internal::RTreeWriter.write(@io, @index_spool.io,
          @section_count, @options.index_block_size)
        write_header(Array(Internal::ZoomLevel).new)
        @io.flush
        zooms = build_zooms
        Internal.put(@io, Internal::BW_MAGIC)
        write_header(zooms)
        @io.flush
        @io.close if @sync_close
        publish
        @closed = true
      rescue ex
        abort_after_failure
        raise ex
      ensure
        @index_spool.close
      end
    end

    private def publish : Nil
      target = @destination
      temporary = @temporary_path
      return unless target && temporary
      if @options.overwrite
        ::File.rename(temporary, target)
      else
        # Unlike exists? + rename, link is an atomic no-clobber operation.
        # A filesystem without hard links fails rather than overwriting.
        ::File.link(temporary, target)
        ::File.delete(temporary)
      end
      @temporary_path = nil
    end

    # Does not finalize. Existing caller-owned IO remains open and may contain
    # incomplete bytes; a path-backed destination is not published.
    def abort : Nil
      return if @closed
      @aborted = true
      @closed = true
      begin
        @io.close if @sync_close && !@io.closed?
      ensure
        @index_spool.close
        if temporary = @temporary_path
          ::File.delete(temporary) if ::File.exists?(temporary)
          @temporary_path = nil
        end
      end
    end
  end
end

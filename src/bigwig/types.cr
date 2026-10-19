module BigWig
  # Coordinates are zero-based; stop is excluded. Returned records retain the
  # original interval boundaries, not boundaries clipped to the query.
  # This remains a value type because readers can emit millions of intervals.
  record Interval, start : Int64, stop : Int64, value : Float32 do
    def initialize(start : Int, stop : Int, value : Int | Float)
      @start = start.to_i64
      @stop = stop.to_i64
      @value = value.to_f32
    end
  end

  # `rest` owns the tab-separated fields after chrom/start/stop. It remains
  # valid after advancing the iterator or closing the reader.
  # This remains a value type because readers create one for every matching entry.
  record BedEntry, start : Int64, stop : Int64, rest : String

  enum FileType
    BigWig
    BigBed
  end

  # The stored file-wide summary is absent in some older files. In that case
  # has_summary? is false, bases_covered is zero, and numeric summaries are NaN.
  class Header
    getter version : Int32
    getter zoom_levels : Int32
    getter bases_covered : UInt64
    getter min_value : Float64
    getter max_value : Float64
    getter sum_data : Float64
    getter sum_squared : Float64
    getter field_count : Int32
    getter defined_field_count : Int32
    # ameba:disable Naming/QueryBoolMethods
    getter has_summary : Bool

    # ameba:enable Naming/QueryBoolMethods

    def initialize(@version, @zoom_levels, @bases_covered, @min_value, @max_value,
                   @sum_data, @sum_squared, @field_count, @defined_field_count,
                   @has_summary : Bool)
    end

    def has_summary? : Bool
      @has_summary
    end

    def copy_with(version _version : Int32 = @version,
                  zoom_levels _zoom_levels : Int32 = @zoom_levels,
                  bases_covered _bases_covered : UInt64 = @bases_covered,
                  min_value _min_value : Float64 = @min_value,
                  max_value _max_value : Float64 = @max_value,
                  sum_data _sum_data : Float64 = @sum_data,
                  sum_squared _sum_squared : Float64 = @sum_squared,
                  field_count _field_count : Int32 = @field_count,
                  defined_field_count _defined_field_count : Int32 = @defined_field_count,
                  has_summary _has_summary : Bool = @has_summary) : self
      self.class.new(_version, _zoom_levels, _bases_covered, _min_value, _max_value,
        _sum_data, _sum_squared, _field_count, _defined_field_count, _has_summary)
    end

    def_equals_and_hash @version, @zoom_levels, @bases_covered, @min_value, @max_value,
      @sum_data, @sum_squared, @field_count, @defined_field_count, @has_summary
  end

  # Allocation limits are explicit and can be raised for unusual files.
  # cache_bytes bounds the HTTP LRU cache, not decompressed query buffers.
  class ReadOptions
    getter max_block_bytes : Int32
    getter max_array_bytes : Int64
    getter http_block_bytes : Int32
    getter cache_bytes : Int32
    getter connect_timeout : Time::Span
    getter read_timeout : Time::Span

    def initialize(@max_block_bytes = 64 * 1024 * 1024,
                   @max_array_bytes = 256_i64 * 1024 * 1024,
                   @http_block_bytes = 64 * 1024,
                   @cache_bytes = 1024 * 1024,
                   @connect_timeout = 10.seconds,
                   @read_timeout = 30.seconds)
    end

    def validate : Nil
      raise ArgumentError.new("max_block_bytes must be positive") unless max_block_bytes > 0
      raise ArgumentError.new("max_array_bytes must be positive") unless max_array_bytes > 0
      unless 1 <= http_block_bytes <= max_block_bytes && cache_bytes >= http_block_bytes
        raise ArgumentError.new("require 1 <= http_block_bytes <= max_block_bytes and cache_bytes >= http_block_bytes")
      end
      unless connect_timeout > 0.seconds && read_timeout > 0.seconds
        raise ArgumentError.new("HTTP timeouts must be positive")
      end
    end

    def copy_with(max_block_bytes _max_block_bytes : Int32 = @max_block_bytes,
                  max_array_bytes _max_array_bytes : Int64 = @max_array_bytes,
                  http_block_bytes _http_block_bytes : Int32 = @http_block_bytes,
                  cache_bytes _cache_bytes : Int32 = @cache_bytes,
                  connect_timeout _connect_timeout : Time::Span = @connect_timeout,
                  read_timeout _read_timeout : Time::Span = @read_timeout) : self
      self.class.new(_max_block_bytes, _max_array_bytes, _http_block_bytes, _cache_bytes,
        _connect_timeout, _read_timeout)
    end

    def_equals_and_hash @max_block_bytes, @max_array_bytes, @http_block_bytes, @cache_bytes,
      @connect_timeout, @read_timeout
  end

  # max_zooms is a maximum: empty files and exceptionally large numeric values
  # can produce fewer levels. Temporary files keep index/zoom memory bounded.
  class WriteOptions
    # ameba:disable Naming/QueryBoolMethods
    getter compression : Bool
    # ameba:enable Naming/QueryBoolMethods
    getter compression_level : Int32
    getter block_items : Int32
    getter index_block_size : Int32
    getter max_zooms : Int32
    # ameba:disable Naming/QueryBoolMethods
    getter overwrite : Bool

    # ameba:enable Naming/QueryBoolMethods

    def initialize(@compression = true,
                   @compression_level = 6,
                   @block_items = 1024,
                   @index_block_size = 64,
                   @max_zooms = 10,
                   @overwrite = false)
    end

    def validate : Nil
      raise ArgumentError.new("compression_level must be in 0..9") unless 0 <= compression_level <= 9
      raise ArgumentError.new("block_items must be in 1..65535") unless 1 <= block_items <= 65_535
      raise ArgumentError.new("index_block_size must be in 2..65535") unless 2 <= index_block_size <= 65_535
      raise ArgumentError.new("max_zooms must be in 0..16") unless 0 <= max_zooms <= 16
    end

    def copy_with(compression _compression : Bool = @compression,
                  compression_level _compression_level : Int32 = @compression_level,
                  block_items _block_items : Int32 = @block_items,
                  index_block_size _index_block_size : Int32 = @index_block_size,
                  max_zooms _max_zooms : Int32 = @max_zooms,
                  overwrite _overwrite : Bool = @overwrite) : self
      self.class.new(_compression, _compression_level, _block_items, _index_block_size,
        _max_zooms, _overwrite)
    end

    def_equals_and_hash @compression, @compression_level, @block_items, @index_block_size,
      @max_zooms, @overwrite
  end
end

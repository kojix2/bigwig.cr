module BigWig::Internal
  # :nodoc:
  # Weighted Welford variance and compensated sums. A constant interval counts
  # as (stop-start) observations without visiting each base.
  struct Aggregate
    getter count = 0.0
    getter sum = 0.0
    getter sum_squared = 0.0
    getter min = Float64::INFINITY
    getter max = -Float64::INFINITY
    @mean = 0.0
    @m2 = 0.0
    @sum_error = 0.0
    @squares_error = 0.0

    def add(value : Float64, weight : Float64) : Nil
      return if weight <= 0
      merge(weight, value * weight, value * value * weight, value, value, 0.0)
    end

    def add(record : BigWig::Interval, first : Int64, last : Int64) : Nil
      add(record.value.to_f64, (last - first).to_f64)
    end

    def add(record : ZoomRecord, first : Int64, last : Int64) : Nil
      return if first >= last || record.count == 0
      scalar = (last - first).to_f64 / (record.stop - record.start)
      count = record.count.to_f64 * scalar
      sum = record.sum.to_f64 * scalar
      squares = record.sum_squared.to_f64 * scalar
      m2 = Math.max(0.0, squares - sum * sum / count)
      merge(count, sum, squares, record.min.to_f64, record.max.to_f64, m2)
    end

    private def merge(count : Float64, sum : Float64, squares : Float64,
                      minimum : Float64, maximum : Float64, m2 : Float64) : Nil
      if @count == 0
        @mean = sum / count
        @m2 = m2
        @min = minimum
        @max = maximum
      else
        delta = sum / count - @mean
        total = @count + count
        @m2 += m2 + delta * delta * (@count / total) * count
        @mean += delta * (count / total)
        @min = Math.min(@min, minimum)
        @max = Math.max(@max, maximum)
      end
      @count += count
      y = sum - @sum_error
      t = @sum + y
      @sum_error = (t - @sum) - y
      @sum = t
      y2 = squares - @squares_error
      t2 = @sum_squared + y2
      @squares_error = (t2 - @sum_squared) - y2
      @sum_squared = t2
    end

    def result(type : Symbol, width : Int64, approximate : Bool) : Float64
      return Float64::NAN if width == 0 || @count == 0
      case type
      when :mean     then @sum / @count
      when :min      then @min
      when :max      then @max
      when :sum      then @sum
      when :coverage then @count / width
      when :std
        return approximate ? Float64::NAN : 0.0 if @count <= 1
        Math.sqrt(Math.max(0.0, @m2 / (@count - 1)))
      else
        raise ArgumentError.new("unknown statistic: #{type}")
      end
    end
  end
end

module BigWig
  class Reader
    # Always returns bins Float64s. Empty bins (including coverage and sum)
    # are NaN. Mean excludes uncovered bases; coverage divides by bin width.
    # std is the sample standard deviation (n-1), with one exact base => 0.
    #
    # With exact: false, use the coarsest zoom not exceeding half a bin's
    # width. Partly intersecting zoom records are weighted by overlap/width;
    # minima/maxima may include data just outside the requested bin.
    def stats(chrom : String, start : Int = 0, stop : Int? = nil, *,
              type : Symbol = :mean, bins : Int = 1, exact : Bool = false) : Array(Float64)
      require_bigwig
      unless {:mean, :min, :max, :sum, :coverage, :std}.includes?(type)
        raise ArgumentError.new("type must be :mean, :min, :max, :sum, :coverage, or :std")
      end
      raise ArgumentError.new("bins must be in 1..Int32::MAX") unless 1 <= bins <= Int32::MAX
      check_allocation(bins.to_i64 * sizeof(Float64))
      query = region(chrom, start, stop)
      width = query.stop - query.start
      zoom = exact ? nil : @zoom_levels.reverse.find { |level| level.reduction <= width // bins // 2 }
      if zoom
        iterator = Internal::ZoomIterator.new(@source, blocks(query, zoom.index_offset),
          loader, query, @big_endian, @chrom_sizes)
        summarize(iterator, query, bins.to_i, type, true)
      else
        summarize(interval_iterator(query), query, bins.to_i, type, false)
      end
    end

    def stats(chrom : String, range : Range(B, E), **options) forall B, E
      first, last = Internal.range_bounds(range)
      stats(chrom, first, last, **options)
    end

    # The iterator is consumed once for all bins. It retains at most one
    # interval straddling the next boundary, not an array of all overlaps.
    private def summarize(iterator, query : Internal::Region, bins : Int32,
                          type : Symbol, approximate : Bool) : Array(Float64)
      current = iterator.next
      width = query.stop - query.start
      Array(Float64).new(bins) do |i|
        first = query.start + width * i.to_i64 // bins
        last = query.start + width * (i.to_i64 + 1) // bins
        acc = Internal::Aggregate.new
        unless first == last
          loop do
            item = current
            break if item.is_a?(Iterator::Stop)
            break if item.stop > first
            current = iterator.next
          end
          loop do
            item = current
            break if item.is_a?(Iterator::Stop)
            break if item.start >= last
            acc.add(item, Math.max(first, item.start), Math.min(last, item.stop))
            break if item.stop > last
            current = iterator.next
          end
        end
        acc.result(type, last - first, approximate)
      end
    end

    def mean(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :mean, exact: exact)[0]
    end

    def mean(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      mean(chrom, first, last, exact: exact)
    end

    def min(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :min, exact: exact)[0]
    end

    def min(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      min(chrom, first, last, exact: exact)
    end

    def max(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :max, exact: exact)[0]
    end

    def max(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      max(chrom, first, last, exact: exact)
    end

    def sum(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :sum, exact: exact)[0]
    end

    def sum(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      sum(chrom, first, last, exact: exact)
    end

    def coverage(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :coverage, exact: exact)[0]
    end

    def coverage(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      coverage(chrom, first, last, exact: exact)
    end

    def std(chrom : String, start : Int = 0, stop : Int? = nil, *, exact : Bool = false) : Float64
      stats(chrom, start, stop, type: :std, exact: exact)[0]
    end

    def std(chrom : String, range : Range(B, E), *, exact : Bool = false) : Float64 forall B, E
      first, last = Internal.range_bounds(range)
      std(chrom, first, last, exact: exact)
    end
  end
end

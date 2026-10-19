module BigWig
  class Error < Exception
  end

  # Malformed, truncated, or inconsistent file content.
  class FormatError < Error
  end

  class ClosedError < Error
  end

  class UnsupportedError < Error
  end

  class LimitError < Error
  end

  class HTTPError < Error
  end

  class RangeNotSupportedError < HTTPError
  end

  class SourceChangedError < HTTPError
  end
end

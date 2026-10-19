require "http/client"
require "uri"

module BigWig::Internal
  # :nodoc:
  # Bounded range cache. No whole-file download fallback, and TLS verification
  # is left enabled. Close connections before rejecting malformed responses so
  # HTTP::Client cannot drain an unexpectedly large response body.
  class HTTPSource < Source
    private class Redirect < Exception
      getter location : String

      def initialize(@location : String)
        super("HTTP redirect")
      end
    end

    getter size : Int64 = 0_i64
    @client : HTTP::Client
    @uri : URI
    @closed = false
    @mutex = Mutex.new
    @blocks = Hash(Int64, Bytes).new
    @etag : String? = nil
    @modified : String? = nil
    @reuse = Bytes.empty

    def initialize(url : String, @options : BigWig::ReadOptions,
                   headers : HTTP::Headers = HTTP::Headers.new)
      @uri = URI.parse(url)
      validate_uri(@uri)
      @headers = headers.dup
      @client = make_client(@uri)
      begin
        fetch(0_i64, 1, discover: true)
      rescue ex
        @client.close
        @closed = true
        raise ex
      end
    end

    def closed? : Bool
      @closed
    end

    private def validate_uri(uri : URI) : Nil
      unless {"http", "https"}.includes?(uri.scheme) && uri.host
        raise ArgumentError.new("expected an absolute HTTP(S) URL")
      end
      if uri.user || uri.password
        raise ArgumentError.new("use HTTP headers instead of credentials in URLs")
      end
    end

    private def make_client(uri : URI) : HTTP::Client
      client = HTTP::Client.new(uri)
      client.connect_timeout = @options.connect_timeout
      client.read_timeout = @options.read_timeout
      client.write_timeout = @options.read_timeout
      client
    end

    def read_exact_at(offset : Int64, buffer : Bytes) : Nil
      @mutex.synchronize do
        check_range(offset, buffer.size)
        copied = 0
        while copied < buffer.size
          position = offset + copied
          first = position // @options.http_block_bytes * @options.http_block_bytes
          data = @blocks.delete(first)
          unless data
            if @blocks.size >= @options.cache_bytes // @options.http_block_bytes
              _, @reuse = @blocks.shift
            end
            length = Math.min(@options.http_block_bytes.to_i64, @size - first).to_i
            data = fetch(first, length)
            @reuse = Bytes.empty
          end
          @blocks[first] = data
          from = (position - first).to_i
          count = Math.min(buffer.size - copied, data.size - from)
          buffer[copied, count].copy_from(data[from, count])
          copied += count
        end
      end
    end

    private def read_response(response : HTTP::Client::Response, first : Int64,
                              length : Int32, discover : Bool) : Bytes
      validate_status(response)
      encoding = response.headers["Content-Encoding"]?
      if encoding && encoding.downcase != "identity"
        raise BigWig::HTTPError.new("encoded HTTP ranges are not supported")
      end
      start, stop, total = parse_content_range(response, first, length)
      unless start == first && stop == first + length - 1 && total > first + length - 1
        raise BigWig::HTTPError.new("incorrect byte-range bounds")
      end
      validate_remote_file(response, total, discover)
      if count = response.headers["Content-Length"]?
        raise BigWig::HTTPError.new("incorrect Content-Length") unless count.to_i64? == length
      end
      read_response_body(response, length)
    end

    private def validate_status(response : HTTP::Client::Response) : Nil
      status = response.status_code
      if {301, 302, 303, 307, 308}.includes?(status)
        location = response.headers["Location"]? || raise BigWig::HTTPError.new("redirect has no Location")
        raise Redirect.new(location)
      end
      raise BigWig::RangeNotSupportedError.new("server ignored HTTP Range") if status == 200
      raise BigWig::SourceChangedError.new("remote file changed (HTTP #{status})") if status == 412
      raise BigWig::HTTPError.new("HTTP #{status} for range request") unless status == 206
    end

    private def parse_content_range(response : HTTP::Client::Response, first : Int64,
                                    length : Int32) : Tuple(Int64, Int64, Int64)
      content_range = response.headers["Content-Range"]? || raise BigWig::HTTPError.new("missing Content-Range")
      match = /\Abytes (\d+)-(\d+)\/(\d+)\z/.match(content_range) || raise BigWig::HTTPError.new("invalid Content-Range")
      start = match[1].to_i64? || raise BigWig::HTTPError.new("invalid Content-Range start")
      stop = match[2].to_i64? || raise BigWig::HTTPError.new("invalid Content-Range stop")
      total = match[3].to_i64? || raise BigWig::HTTPError.new("invalid remote file size")
      {start, stop, total}
    end

    private def validate_remote_file(response : HTTP::Client::Response, total : Int64,
                                     discover : Bool) : Nil
      if !discover && total != @size
        raise BigWig::SourceChangedError.new("remote file size changed")
      end
      if (etag = @etag) && response.headers["ETag"]? != etag
        raise BigWig::SourceChangedError.new("remote ETag changed or disappeared")
      end
      if (modified = @modified) && response.headers["Last-Modified"]? != modified
        raise BigWig::SourceChangedError.new("remote Last-Modified changed or disappeared")
      end
      if discover
        @size = total
        @etag = response.headers["ETag"]?
        @modified = response.headers["Last-Modified"]?
      end
    end

    private def read_response_body(response : HTTP::Client::Response, length : Int32) : Bytes
      bytes = @reuse.size >= length ? @reuse[0, length] : Bytes.new(length)
      response.body_io.read_fully(bytes)
      unless response.body_io.read_byte.nil?
        raise BigWig::HTTPError.new("HTTP range body is longer than declared")
      end
      bytes
    end

    private def request_once(first : Int64, length : Int32, discover : Bool,
                             headers : HTTP::Headers) : Bytes
      rejected : Exception? = nil
      request = HTTP::Request.new("GET", @uri.request_target, headers)
      begin
        @client.exec(request) do |response|
          read_response(response, first, length, discover)
        rescue ex
          rejected = ex
          @client.close
          raise ex
        end
      rescue ex
        # Closing a malformed response can itself raise. Keep the original
        # protocol error (or redirect), not a secondary cleanup exception.
        raise(rejected || ex)
      end
    end

    private def fetch(first : Int64, length : Int32, discover : Bool = false) : Bytes
      redirects = 0
      loop do
        headers = @headers.dup
        headers["Range"] = "bytes=#{first}-#{first + length - 1}"
        headers["Accept-Encoding"] = "identity"
        if (etag = @etag) && !etag.starts_with?("W/")
          headers["If-Match"] = etag
        elsif modified = @modified
          headers["If-Unmodified-Since"] = modified
        end
        begin
          return request_once(first, length, discover, headers)
        rescue ex : Redirect
          redirects += 1
          raise BigWig::HTTPError.new("too many HTTP redirects") if redirects > 5
          next_uri = @uri.resolve(ex.location)
          validate_uri(next_uri)
          if @uri.scheme == "https" && next_uri.scheme != "https"
            raise BigWig::HTTPError.new("refusing HTTPS-to-HTTP redirect")
          end
          if {@uri.scheme, @uri.host, @uri.port} != {next_uri.scheme, next_uri.host, next_uri.port}
            # Do not leak Authorization, cookies, or custom token headers.
            @headers = HTTP::Headers.new
          end
          @uri = next_uri
          @client = make_client(@uri)
        end
      end
    rescue ex
      @client.close
      @closed = true
      @blocks.clear
      @reuse = Bytes.empty
      raise ex
    end

    def close : Nil
      @mutex.synchronize do
        return if @closed
        @closed = true
        @client.close
        @blocks.clear
        @reuse = Bytes.empty
      end
    end
  end
end

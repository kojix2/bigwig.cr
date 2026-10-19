require "./spec_helper"
require "http/server"

private def with_http_server(handler : HTTP::Server::Context ->, &)
  server = HTTP::Server.new(&handler)
  address = server.bind_tcp("127.0.0.1", 0)
  spawn { server.listen }
  begin
    yield "http://127.0.0.1:#{address.port}"
  ensure
    server.close
  end
end

private def serve_range(context : HTTP::Server::Context, data : Bytes, etag : String? = nil) : Nil
  match = /\Abytes=(\d+)-(\d+)\z/.match(context.request.headers["Range"]) || raise "invalid test request"
  first, last = match[1].to_i, match[2].to_i
  context.response.status_code = 206
  context.response.headers["Content-Range"] = "bytes #{first}-#{last}/#{data.size}"
  context.response.headers["ETag"] = etag if etag
  context.response.content_length = last - first + 1
  context.response.write(data[first, last - first + 1])
end

describe "HTTP range sources" do
  it "uses bounded ranges, caches hits, and works through the public API" do
    data = File.read(fixture("test.bw")).to_slice
    requests = 0
    handler = ->(context : HTTP::Server::Context) do
      requests += 1
      context.request.headers["Accept-Encoding"].should eq("identity")
      serve_range(context, data, "\"stable\"")
    end
    with_http_server(handler) do |base|
      options = BigWig::ReadOptions.new(http_block_bytes: 1024, cache_bytes: 65536)
      BigWig.open("#{base}/test.bw", options: options) do |reader|
        reader.value("1", 0).should eq(0.1_f32)
        previous = requests
        reader.value("1", 0).should eq(0.1_f32)
        requests.should eq(previous)
        reader.intervals("1").size.should eq(5)
      end
      requests.should be < 30
    end
  end

  it "rejects a server that ignores Range instead of downloading the whole file" do
    handler = ->(context : HTTP::Server::Context) do
      context.response.status_code = 200
      context.response.print "not a range"
      nil
    end
    with_http_server(handler) do |base|
      expect_raises(BigWig::RangeNotSupportedError) { BigWig.open("#{base}/file.bw") }
    end
  end

  it "rejects wrong Content-Range and encoded response bodies" do
    {false, true}.each do |encoded|
      handler = ->(context : HTTP::Server::Context) do
        context.response.status_code = 206
        context.response.headers["Content-Range"] = encoded ? "bytes 0-0/10" : "bytes 1-1/10"
        context.response.headers["Content-Encoding"] = "gzip" if encoded
        context.response.content_length = 1
        context.response.write_byte(0_u8)
        nil
      end
      with_http_server(handler) do |base|
        expect_raises(BigWig::HTTPError) { BigWig.open("#{base}/file.bw") }
      end
    end
  end

  it "rejects missing range metadata and unexpected statuses" do
    {206, 404, 416}.each do |status|
      handler = ->(context : HTTP::Server::Context) do
        context.response.status_code = status
        context.response.content_length = 0
        nil
      end
      with_http_server(handler) do |base|
        expect_raises(BigWig::HTTPError) { BigWig.open("#{base}/file.bw") }
      end
    end
  end

  it "detects ETag changes on cache misses and sends If-Match" do
    data = File.read(fixture("test.bw")).to_slice
    etag = "\"first\""
    conditional = false
    handler = ->(context : HTTP::Server::Context) do
      conditional ||= context.request.headers["If-Match"]? == "\"first\""
      serve_range(context, data, etag)
    end
    with_http_server(handler) do |base|
      options = BigWig::ReadOptions.new(http_block_bytes: 64, cache_bytes: 64)
      source = BigWig::Internal::HTTPSource.new("#{base}/file.bw", options)
      begin
        bytes = Bytes.new(4)
        source.read_exact_at(0_i64, bytes)
        conditional.should be_true
        etag = "\"second\""
        expect_raises(BigWig::SourceChangedError) { source.read_exact_at(128_i64, bytes) }
        source.closed?.should be_true
      ensure
        source.close
      end
    end
  end

  it "handles relative redirects" do
    data = File.read(fixture("test.bw")).to_slice
    handler = ->(context : HTTP::Server::Context) do
      if context.request.path == "/start"
        context.response.status_code = 302
        context.response.headers["Location"] = "/actual"
        context.response.content_length = 0
      else
        serve_range(context, data)
      end
      nil
    end
    with_http_server(handler) do |base|
      BigWig.open("#{base}/start") { |reader| reader.value("1", 1).should eq(0.2_f32) }
    end
  end

  it "does not forward custom secrets to a different redirect origin" do
    data = File.read(fixture("test.bw")).to_slice
    destination = ->(context : HTTP::Server::Context) do
      context.request.headers["X-Secret"]?.should be_nil
      serve_range(context, data)
    end
    with_http_server(destination) do |target|
      redirect = ->(context : HTTP::Server::Context) do
        context.response.status_code = 302
        context.response.headers["Location"] = "#{target}/data"
        context.response.content_length = 0
        nil
      end
      with_http_server(redirect) do |base|
        BigWig.open("#{base}/data", headers: HTTP::Headers{"X-Secret" => "private"}) do |reader|
          reader.value("1", 0).should eq(0.1_f32)
        end
      end
    end
  end

  it "bounds redirect loops" do
    handler = ->(context : HTTP::Server::Context) do
      context.response.status_code = 302
      context.response.headers["Location"] = "/again"
      context.response.content_length = 0
      nil
    end
    with_http_server(handler) do |base|
      expect_raises(BigWig::HTTPError, "too many") { BigWig.open("#{base}/again") }
    end
  end
end

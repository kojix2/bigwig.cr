module BigWig::Internal
  # :nodoc:
  # Temporary spools are scoped explicitly. No finalizer performs I/O.
  class TempFile
    getter io : ::File

    def initialize
      @io = ::File.tempfile("bigwig", ".tmp")
    end

    def close : Nil
      @io.close unless @io.closed?
    ensure
      ::File.delete(@io.path) if ::File.exists?(@io.path)
    end
  end

  # :nodoc:
  # Only chromosome metadata is resident; chromosome IDs retain declaration
  # order while B+tree keys are lexically sorted as required by the format.
  module ChromTreeWriter
    class Node
      getter first : Int32
      getter count : Int32
      getter children : Array(Node)
      property offset : Int64 = 0_i64

      def initialize(@first, @count, @children = Array(Node).new)
      end

      def leaf? : Bool
        @children.empty?
      end
    end

    def self.write(io : IO, chromosomes : Array(Chromosome), fanout : Int32) : Int64
      offset = io.pos.to_i64
      sorted = chromosomes.sort_by(&.name)
      key_size = sorted.max_of(&.name.bytesize) + 1
      put(io, BPT_MAGIC)
      put(io, fanout.to_u32)
      put(io, key_size.to_u32)
      put(io, 8_u32)
      put(io, sorted.size.to_u64)
      put(io, 0_u64)
      level = Array(Node).new
      first = 0
      while first < sorted.size
        count = Math.min(fanout, sorted.size - first)
        level << Node.new(first, count)
        first += count
      end
      while level.size > 1
        parents = Array(Node).new
        level.each_slice(fanout) do |children|
          parents << Node.new(children.first.first, children.size, children)
        end
        level = parents
      end
      root = level.first
      queue = [root]
      position = io.pos.to_i64
      index = 0
      while index < queue.size
        node = queue[index]
        node.offset = position
        position += 4_i64 + node.count.to_i64 * (key_size + 8)
        queue.concat(node.children)
        index += 1
      end
      key = Bytes.new(key_size)
      queue.each do |queued_node|
        put(io, queued_node.leaf? ? 1_u8 : 0_u8)
        put(io, 0_u8)
        put(io, queued_node.count.to_u16)
        queued_node.count.times do |i|
          chrom = queued_node.leaf? ? sorted[queued_node.first + i] : sorted[queued_node.children[i].first]
          key.fill(0_u8)
          key[0, chrom.name.bytesize].copy_from(chrom.name.to_slice)
          io.write(key)
          if queued_node.leaf?
            put(io, chrom.id)
            put(io, chrom.size.to_u32)
          else
            put(io, queued_node.children[i].offset.to_u64)
          end
        end
      end
      offset
    end

    private def self.put(io : IO, value) : Nil
      Internal.put(io, value)
    end
  end

  # :nodoc:
  # R-tree construction uses a disk spool of 32-byte descriptors. Leaf nodes
  # are written first, then parent levels. A reserved root at index+48 points
  # to them using absolute offsets. Resident descriptors are O(fanout).
  module RTreeWriter
    def self.write(io : IO, input : IO, count : Int64, fanout : Int32,
                   items_per_slot : Int32 = 1) : Int64
      offset = io.pos.to_i64
      Internal.zeros(io, 48_i64 + 4 + fanout.to_i64 * 32)
      owned = Array(TempFile).new
      begin
        current = input
        current_count = count
        leaf = true
        while current_count > fanout
          next_level = TempFile.new
          owned << next_level
          current.rewind
          groups = 0_i64
          remaining = current_count
          while remaining > 0
            n = Math.min(fanout.to_i64, remaining).to_i
            descriptor = write_node(io, current, n, leaf)
            descriptor.write(next_level.io)
            groups += 1
            remaining -= n
          end
          if previous = owned[-2]?
            previous.close
          end
          current = next_level.io
          current_count = groups
          leaf = false
        end
        end_position = io.pos
        current.rewind
        io.seek(offset + 48)
        root = write_node(io, current, current_count.to_i, leaf)
        io.seek(offset)
        Internal.put(io, RTREE_MAGIC)
        Internal.put(io, fanout.to_u32)
        Internal.put(io, count.to_u64)
        Internal.put(io, root.start_chrom)
        Internal.put(io, root.start.to_u32)
        Internal.put(io, root.end_chrom)
        Internal.put(io, root.stop.to_u32)
        Internal.put(io, offset.to_u64) # end of indexed data
        Internal.put(io, items_per_slot.to_u32)
        Internal.put(io, 0_u32)
        io.seek(end_position)
        offset
      ensure
        owned.each(&.close)
      end
    end

    private def self.write_node(io : IO, descriptors : IO, count : Int32, leaf : Bool) : BlockRef
      offset = io.pos.to_i64
      Internal.put(io, leaf ? 1_u8 : 0_u8)
      Internal.put(io, 0_u8)
      Internal.put(io, count.to_u16)
      first : BlockRef? = nil
      last : BlockRef? = nil
      count.times do
        item = BlockRef.read(descriptors)
        first ||= item
        last = item
        Internal.put(io, item.start_chrom)
        Internal.put(io, item.start.to_u32)
        Internal.put(io, item.end_chrom)
        Internal.put(io, item.stop.to_u32)
        Internal.put(io, item.offset.to_u64)
        Internal.put(io, item.size.to_u64) if leaf
      end
      if a = first
        b = last || raise BigWig::FormatError.new("missing final R-tree block")
        BlockRef.new(a.start_chrom, a.start, b.end_chrom, b.stop, offset, 0_i64)
      else
        BlockRef.new(0_u32, 0_i64, 0_u32, 0_i64, offset, 0_i64)
      end
    end
  end
end

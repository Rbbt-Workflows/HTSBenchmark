require 'HTSBenchmark/SV'
module HTSBenchmark
  def self.slice_fasta(input, output, ranges = {})

    fragments = HTSBenchmark.collect_fragments(input, ranges)

    Open.write(output) do |f|
      fragments.each do |chr, fragments|
        f.puts ">#{chr}"
        f.puts fragments.sort_by{|range,s| range.first}.collect{|range,s| s } * ""
      end
    end
  end

  def self.slice_vcf(input, output, ranges = {})

    sorted_range_with_offset = {}
    sizes = {}
    ranges.each do |chr, list|
      sorted_list = list.sort_by{|s,e| s }
      list_with_offsets = []
      offset = 0
      sorted_list.each do |s,e|
        list_with_offsets << [s, e, offset]
        offset += e - s + 1
      end
      chr = chr.sub(/^chr/,'')
      sizes[chr] = offset
      sorted_range_with_offset[chr] = list_with_offsets
    end

    input = StringIO.new input if String  === input && ! Misc.is_filename?(input)
    Open.open(output, :mode => 'w') do |sout|
      chr = nil
      size = nil
      pointer = nil
      TSV.traverse input, :type => :array, :bar => true do |line|
        if line[0] == '#'
          if m = line.match(/##contig=<ID=(.*),length=.*/)
            chr = m.captures.first
            size = sizes[chr.sub(/^chr/, '')]
            #sout.puts "##contig=<ID=#{chr},length=#{size}>" if size
          else
            sout.puts line
          end
        else
          orig_chr, _sep, position_and_rest = line.partition("\t")
          position, _sep, rest = position_and_rest.partition("\t")

          position = position.to_i
          chr = orig_chr.sub(/^chr/, '')
          sorted_range_with_offset[chr].each do |start,eend,offset|
            if start <= position && position <= eend
              position = position - start + offset + 1
              sout.puts [orig_chr, position.to_s, rest] * "\t"
              break
            end
          end
        end
      end
    end
  end

end

module HTSBenchmark
  def self.minify_fasta(input, output, sizes = {})
    padding = sizes[:padding]

    Open.open(output, :mode => 'w') do |sout|
      chr = nil
      size = nil
      pointer = nil
      TSV.traverse input, :type => :array, :bar => true do |line|
        if line[0] == '>'
          chr = line.split(" ").first[1..-1]
          chr = chr.sub(/^chr/, '')
          size = sizes[chr]
          size = size + padding if size && padding
          sout.puts line.split(" ").first if size
          pointer = 0
        else
          sout.puts line[0..size - pointer - 1] if size && size > pointer
          pointer = pointer + line.size 
        end
      end
    end
  end

  def self.minify_vcf(input, output, sizes = {})

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
          chr, _sep, position = line.partition("\t")
          chr = chr.sub(/^chr/, '')
          size = sizes[chr]
          next if size.nil?
          position = position.to_i
          sout.puts line if size > position
        end
      end
    end
  end

  def self.minify_mutations(input, output, sizes = {})

    input = StringIO.new input if String  === input && ! Misc.is_filename?(input)
    Open.open(output, :mode => 'w') do |sout|
      chr = nil
      size = nil
      pointer = nil
      TSV.traverse input, :type => :array, :bar => true do |line|
        chr, _sep, position = line.partition(":")
        chr = chr.sub(/^chr/, '')
        size = sizes[chr]
        next if size.nil?
        position = position.to_i
        sout.puts line if size > position
      end
    end
  end


  def self.calculate_sizes(positions, min = 100, padding = 1_000)
    sizes = {}
    positions = StringIO.new(positions) if String === positions
    TSV.traverse positions, :type => :array, :bar => true do |position|
      chr, _sep, position = position.partition(":")
      chr = chr.sub(/^chr/, '')
      position = position.to_i
      sizes[chr] ||= [[], nil]

      list, max = sizes[chr]

      next if max && max < position && list.length == min

      list.push(position)
      list.sort!
      list.pop if list.length > min
      max = list.last
      sizes[chr] = [list, max]
    end
    res = {}
    sizes.collect{|k,v| res[k] = v.last + padding}
    res
  end
end

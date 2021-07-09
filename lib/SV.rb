module HTSBenchmark
  def self.collect_fragments(input, ranges)
    chr = nil
    size = nil
    pointer = nil
    fragments = {}
    chr_txt = nil
    TSV.traverse input, :type => :array, :bar => true do |line|
      line = line.chomp
      if line[0] == '>'
        if chr
          ranges[chr].each do |start,eend,id|
            fragments[chr] ||= {}
            fragments[chr][[start, eend,id]] = chr_txt[start-1..eend-1]
          end if ranges[chr]
        end
        chr = line.split(" ").first[1..-1]
        chr_txt = ""
      else
        chr_txt << line
      end
    end

    if chr
      ranges[chr].each do |start,eend,id|
        fragments[chr] ||= {}
        fragments[chr][[start, eend,id]] = chr_txt[start..eend]
      end if ranges[chr]
    end

    fragments
  end

  def self.remove_fragments(input, inserts = {})
    chr = nil
    size = nil
    pointer = nil
    chr_txt = nil
    Misc.open_pipe do |output|
      TSV.traverse input, :type => :array, :bar => true do |line|
        line = line.chomp
        if line[0] == '>'
          if chr
            offset = 0
            inserts[chr].sort_by{|s,e| s }.reverse.each do |start,eend|
              chr_txt[start-1..eend-1] = "" 
            end if inserts[chr]
            output << chr_txt << "\n"
          end
          output << line << "\n"
          chr = line.split(" ").first[1..-1]
          chr_txt = ""
        else
          chr_txt << line
        end
      end

      if chr
        inserts[chr].sort_by{|s,e| s }.reverse.each do |start,eend|
          chr_txt[start-1..eend-1] = "" 
        end if inserts[chr]
        output << chr_txt << "\n"
      end

    end
  end

  def self.insert_fragments(input, inserts = {})
    chr = nil
    size = nil
    pointer = nil
    chr_txt = nil
    Misc.open_pipe do |output|
      TSV.traverse input, :type => :array, :bar => true do |line|
        line = line.chomp
        if line[0] == '>'
          if chr
            offset = 0
            inserts[chr].sort_by{|s,i| s }.reverse.each do |start,insert|
              chr_txt[start-2..start-2] = chr_txt[start-2] + insert 
            end if inserts[chr]
            output << chr_txt << "\n"
          end
          output << line << "\n"
          chr = line.split(" ").first[1..-1]
          chr_txt = ""
        else
          chr_txt << line #<< "\n"
        end
      end

      if chr
        inserts[chr].sort_by{|s,i| s }.reverse.each do |start,insert|
          chr_txt[start-2..start-2] = chr_txt[start-2] + insert 
        end if inserts[chr]
        output << chr_txt << "\n"
      end
    end
  end

  def self.shift_ranges_by_deletes(ranges, deletes)
    new = {}
    ranges.keys.each do |chr|
      if deletes[chr].nil?
        new[chr] = ranges[chr]
      else
        shifts = deletes[chr].sort_by{|s,e| s }.collect do |s,e|
          [s, e - s + 1]
        end

        new[chr] = ranges[chr].collect{|s,e,id|
          acc = 0
          shifts.each do |ss,sl|
            break if ss > s
            if ss <= s
              if ss + sl >= s
                acc = nil
                break
              else
                acc += sl
              end
            end
          end
          acc.nil? ? [nil, nil, id] : [s - acc, e - acc, id]
        }
      end
    end

    new
  end

  def self.shift_ranges_by_inserts(ranges, inserts)
    new = {}

    ranges.keys.each do |chr|
      if inserts[chr].nil?
        new[chr] = ranges[chr]
      else
        shifts = inserts[chr].sort_by{|s,e| s }

        new[chr] = ranges[chr].collect{|s,e,id|
          e = e.length if String === e
          acc = 0
          shifts.each do |ss,sl|
            sl = sl.length if String === sl
            break if ss > e
            if ss <= s
              acc += sl
            else
              acc = nil
              break
            end
          end
          acc.nil? ? [nil, nil, id] : [s + acc, e + acc, id]
        }
      end
    end

    new
  end

  def self.duplicate_positions(positions, duplications)
    new = {}

    chr_offsets = {}
    duplications.each do |chr,list|
      list.each do|start,eend,tchr,tpos,teend,invert|
        chr_offsets[tchr] ||= []
        chr_offsets[tchr] << [tpos, eend - start + 1]
      end
    end

    new_positions = []
    positions.keys.each do |chr|
      new[chr] = positions[chr]

      if ! duplications[chr].nil?
        fragments = duplications[chr].sort_by{|s| s }
        positions[chr].each do |pos,eend,id|
          fragments.each do |start,eend,tchr,tpos,teend,invert|
            next unless start <= pos && eend >= pos
            chr_offset = Misc.sum(chr_offsets[tchr].sort_by{|t,l| t }.select{|t,l| t < tpos }.collect{|t,l| l }).to_i
            offset = invert ? eend - pos : pos - start
            new_pos = [tchr, chr_offset + tpos + offset, id]
            new_positions <<  new_pos
          end
        end
      end
    end

    new_positions
  end

  def self.apply_SVs(reference, mutations, svs)

    collect_fragments = {}
    remove_fragments = {}
    insert_fragments = []
    duplications = {}

    copies = {}
    TSV.traverse reference, :type => :array do |line|
      if m = line.match(/>((copy-\d+)_([^\s+]*))/)
        orig, copy, real = m.captures
        real = real.sub('chr', '')
        copies[real] ||= []
        copies[real] << orig
      end
    end

    TSV.traverse svs do |id,values|
      type, chr, start, eend, target_chr, target_start, target_end = values

      if !(chr.include?('copy-') || copies[chr.sub('chr', '')].nil?)
        chr_copies = copies[chr.sub('chr', '')]
        num = Misc.digest([chr, start, eend, type] * ":").chars.inject(0){|acc,e| acc += e.hex }
        chr = chr_copies[num % chr_copies.length]
      end

      target_chr = chr if target_chr == 'same' || target_chr == 'cis'

      if target_chr && !(target_chr.include?('copy-') || copies[target_chr.sub('chr', '')].nil?)
        chr_copies = copies[target_chr.sub('chr', '')]
        num = Misc.digest([chr, start, eend, target_chr, type] * ":").chars.inject(0){|acc,e| acc += e.hex }
        target_chr = chr_copies[num % chr_copies.length]
      end

      case type
      when "DEL"
        remove_fragments[chr] ||= [] 
        remove_fragments[chr] << [start.to_i, eend.to_i, id]
      when "INS"
        collect_fragments[id] = [chr, start.to_i, eend.to_i]
        duplications[chr] ||= []
        duplications[chr] << [start.to_i, eend.to_i, target_chr, target_start.to_i, target_end.to_i, false]
        insert_fragments << [id, target_chr, target_start.to_i, target_end.to_i, false]
      when "INV"
        collect_fragments[id] = [chr, start.to_i, eend.to_i]
        duplications[chr] ||= []
        duplications[chr] << [start.to_i, eend.to_i, target_chr, target_start.to_i, target_end.to_i, true]
        insert_fragments << [id, target_chr, target_start.to_i, target_end.to_i, true]
      end
    end

    collect_fragment_ranges = {}
    collect_fragments.each do |id,values|
      chr, start, eend = values
      collect_fragment_ranges[chr] ||= []
      collect_fragment_ranges[chr] << [start, eend, id]
    end

    fragments = {}

    HTSBenchmark.collect_fragments(reference, collect_fragment_ranges).each do |chr,sequences|
      collect_fragment_ranges[chr].zip(sequences).each do |range_info,sequence_info|
        start, eend, id = range_info
        range, sequence = sequence_info

        fragments[id] = sequence
      end
    end

    insertions = {}
    insert_fragments.each do |id,chr,s,e,inverse|
      insertions[chr] ||= []
      sequence = fragments[id]
      next if sequence.nil?
      sequence = sequence.reverse if inverse
      insertions[chr] << [s, sequence, id]
    end

    removals = {}
    remove_fragments.each do |chr,list|
      list.each do |s,e,id|
        removals[chr] ||= []
        removals[chr] << [s,e,id]
      end
    end

    mutations_chr = {}
    mutations.each do |mutation|
      chr, pos, alt = mutation.split(":")

      if !(chr.include?('copy-') || copies[chr.sub('chr', '')].nil?)
        chr_copies = copies[chr.sub('chr', '')]
        num = Misc.digest([chr, pos] * ":").chars.inject(0){|acc,e| acc += e.hex }
        chr = chr_copies[num % chr_copies.length]
      end

      mutations_chr[chr] ||= []
      mutations_chr[chr] << [pos.to_i, pos.to_i, mutation]
    end


    new_positions = HTSBenchmark.duplicate_positions(mutations_chr, duplications)

    mutations_chr = HTSBenchmark.shift_ranges_by_inserts(mutations_chr, insertions)

    new_positions.each do |chr,pos,id|
      mutations_chr[chr] ||= []
      mutations_chr[chr] << [pos, pos, id]
    end

    removals = HTSBenchmark.shift_ranges_by_inserts(removals, insertions)

    mutations_chr = HTSBenchmark.shift_ranges_by_deletes(mutations_chr, removals)

    mutation_translations = TSV.setup({}, :key_field => "Genomic Position", :fields => ["Translation (Genomic Position)"], :type => :flat)

    mutations_chr.each do |chr,list| 
      list.collect{|e| 
        start, eend, id = e
        iii [chr, start, eend, id] if id.include?("SV:")
        mutation_translations[id] ||= []
        next if start.nil?
        _chr, _pos, alt = id.split(":")
        mutation_translations[id] << [chr, start, alt] * ":"
      } 
    end

    reference = HTSBenchmark.insert_fragments(reference,  insertions)

    reference = HTSBenchmark.remove_fragments(reference,  removals)

    [reference, mutation_translations]
  end

end
require 'fc'
module HTSBenchmark

  def self.restore_VCF_positions(vcf, svs, bar = nil)
    collect_fragments, insert_fragments, remove_fragments, duplications = self.SV_regions(svs)

    insertions = {}
    insert_fragments.each do |id,chr,s,e,inverse|
      l = collect_fragments[id][2] - collect_fragments[id][1] + 1
      insertions[chr] ||= []
      insertions[chr] << [s, l]
    end

    last_chr = nil
    delete_shifts = insert_shifts = nil
    TSV.traverse vcf, :into => :stream, :type => :array, :bar => bar do |line|
      next line if line =~ /^#/
      chr, pos, *rest = line.split("\t")

      pos = pos.to_i

      if insertions && last_chr != chr
        delete_shifts = (remove_fragments[chr] || []).sort_by{|s,e| s }.collect{|s,e| [s, e - s + 1] }
        insert_shifts = (insertions[chr] || []).sort_by{|s,e| s }
        last_chr = chr
      end

      last_chr ||= chr

      if delete_shifts
        delete_shifts.each do |s,l|
          pos += l if pos > s
        end

        insert_shifts.each do |s,l|
          pos -= l if pos > s
        end
      end

      ([chr, pos] + rest) * "\t"
    end
  end

  def self.rename_FASTQ_reads(golden_bam, fastq1, fastq2, target_fastq1, target_fastq2, svs, bar = nil)
    if svs
      collect_fragments, insert_fragments, remove_fragments, duplications = self.SV_regions(svs)

      insertions = {}
      insert_fragments.each do |id,chr,s,e,inverse|
        l = collect_fragments[id][2] - collect_fragments[id][1] + 1
        insertions[chr] ||= []
        insertions[chr] << [s, l]
      end

    end

    last_chr = nil
    delete_shifts = insert_shifts = nil
    gb_io = TSV.traverse CMD.cmd(:samtools, "view #{golden_bam} | cut -f 1,2,3,4,8", :pipe => true),
      :bar => bar,
      :type => :array, :into => :stream do |line|
      name, flags, chr, pos, pair_pos = line.split("\t")
      pair = flags.to_i[7]
      next unless pair == 0

      pos = pos.to_i
      pair_pos = pair_pos.to_i

      parts = name.split("-")
      chr = name.scan(/((?:copy-\d+_)?chr(?:.*?))-/).flatten.first

      id = parts.pop

      if insertions && last_chr != chr
        delete_shifts = (remove_fragments[chr] || []).sort_by{|s,e| s }.collect{|s,e| [s, e - s + 1] }
        insert_shifts = (insertions[chr] || []).sort_by{|s,e| s }
        last_chr = chr
      end

      last_chr ||= chr

      if delete_shifts
        delete_shifts.each do |s,l|
          pos += l if pos > s
          pair_pos += l if pair_pos > s
        end

        insert_shifts.each do |s,l|
          pos -= l if pos > s
          pair_pos -= l if pair_pos > s
        end
      end

      loc = [pos.to_s, pair_pos.to_s] * "_"
      new_name = [parts, loc, id] * "-"

      [name, new_name] * "\t"
    end

    lines = FastContainers::PriorityQueue.new(:min)

    current = 1
    current_chr = nil
    sgb_io = Misc.open_pipe do |sin|
      while line = gb_io.gets
        m = line.match(/(.*)-(.*)-(.*)\t/)
        sample, chr, num = m.values_at 1,2,3
        chr = chr.sub('chr','').to_i
        current_chr ||= chr
        if current_chr != chr
          current_chr = chr
          current = 1
        end
        num = (num.to_f + 1.0) / 2 # read numbers are always odd
        lines.push line, num
        while lines.top_key == current
          line = lines.top
          lines.pop
          sin.puts line
          current += 1
        end
      end

      while lines.top_key == current
        line = lines.top
        lines.pop
        sin.puts line
        current += 1
      end
    end

    out1, in1 = Misc.pipe
    out2, in2 = Misc.pipe

    thr1 = Misc.consume_stream Open.gzip(out1), true, target_fastq1
    thr2 = Misc.consume_stream Open.gzip(out2), true, target_fastq2

    Open.open(fastq1) do |f1|
      Open.open(fastq2) do |f2|
        while true
          line1 = f1.gets.strip
          line2 = f2.gets.strip
          gb_line = sgb_io.gets.strip

          name = line1.scan(/^@(.*)\/1/).flatten.first
          gb_name, _sep, new_name = gb_line.partition("\t")

          raise "Read not found #{name}: #{line1}" unless line1.include? name

          line1.sub!(name, new_name) 
          line2.sub!(name, new_name) 

          in1.puts line1
          in1.puts f1.gets
          in1.puts f1.gets
          in1.puts f1.gets

          in2.puts line2
          in2.puts f2.gets
          in2.puts f2.gets
          in2.puts f2.gets

          break if f1.eof?
        end
      end
    end

    in1.close
    in2.close

    thr1.join
    thr2.join
  end

  dep :NEAT_simulate_DNA, :compute => :produce
  task :rename_FASTQ => :string do 
    neat = step(:NEAT_simulate_DNA)

    golden_bam = neat.file('output').glob("*.bam").last
    fastq1 = neat.file('output').glob("*read1.fq.gz").last
    fastq2 = neat.file('output').glob("*read2.fq.gz").last
    target_fastq1 = file('read1.fa.gz')
    target_fastq2 = file('read2.fa.gz')

    HTSBenchmark.rename_FASTQ_reads(golden_bam, fastq1, fastq2, target_fastq1, target_fastq2, nil, self.progress_bar("Annotating reads with position info"))
    
    "DONE"
  end

  input :golden_bam, :file, "", nil, :nofile => true
  input :fastq1, :file, "", nil, :nofile => true
  input :fastq2, :file, "", nil, :nofile => true
  task :rename_FASTQ_reads => :string do |golden_bam,fastq1,fastq2|
    target_fastq1 = file('read1.fa.gz')
    target_fastq2 = file('read2.fa.gz')

    HTSBenchmark.rename_FASTQ_reads(golden_bam, fastq1, fastq2, target_fastq1, target_fastq2, nil, self.progress_bar("Annotating reads with position info"))
    
    "DONE"
  end


end

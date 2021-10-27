require 'HTSBenchmark/miniref'
require 'HTSBenchmark/sliceref'

module HTSBenchmark
  def self.mutation_reference(mutations, reference_file)
    sizes = Samtools.contig_sizes(reference_file)

    contigs = sizes.keys
    copies = {}
    contigs.each do |contig|
      if m = contig.match(/^(copy-\d+)_(?:chr)?(.*)/)
        copy, orig = m.captures
      else
        orig = contig
      end
      copies[orig] ||= []
      copies[orig] << contig
    end

    mutation_ranges = mutations.collect do |m| 
      chr, pos, alt = m.split(":")
      chr = chr.sub(/^chr/,'')

      if alternatives = copies[chr]
        num = Misc.digest([chr, pos] * ":").chars.inject(0){|acc,e| acc += e.hex }
        chr = alternatives[num % alternatives.length]
      end

      chr = "chr" + chr if sizes[chr].nil? && chr !~ /^(copy|chr)/
      chr = chr.sub(/^chr/,'') if sizes[chr].nil? && chr =~ /^chr/

      next if sizes[chr].nil?

      pos = pos.to_i
      if alt["-"]
        pos = pos - 1
        eend = pos + alt.length
      else
        eend = pos
      end


      next if eend > sizes[chr]

      chr + ":" + pos.to_s + "-" + eend.to_s
    end.compact

    reference = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Reference"], :type => :single)

    TmpFile.with_file(mutation_ranges * "\n") do |ranges|
      pos_reference = TSV.setup({}, "Genomic Position~Reference#:type=:single")
      TSV.traverse CMD.cmd(:samtools, "faidx #{reference_file} -r #{ranges} 2> /dev/null | tr '\\n' '\\t' | tr '>' '\\n'", :pipe => true), :type => :array do |line|
        pos_info, ref = line.split("\t")
        next if ref.nil?
        chr, range = pos_info.split(":")
        pos = range.split("-").first
        pos_reference[[chr, pos]] = ref
      end

      mutations.each do |mutation|
        chr, pos, alt = mutation.split(":")

        if alternatives = copies[chr]
          num = Misc.digest([chr, pos] * ":").chars.inject(0){|acc,e| acc += e.hex }
          chr = alternatives[num % alternatives.length]
        end

        chr = "chr" + chr if sizes[chr].nil? && chr !~ /^(copy|chr)/
        chr = chr.sub(/^chr/,'') if sizes[chr].nil? && chr =~ /^chr/

        next if sizes[chr].nil?

        mutation = [chr, pos, alt] * ":"

        pos = pos.to_i
        if alt["-"]
          pos = pos - 1
          eend = pos + alt.length
        else
          eend = pos
        end

        next if eend > sizes[chr]
        
        ref = pos_reference[[chr, pos.to_s]]
        raise mutation if ref.nil?
        reference[mutation] = ref
      end
    end

    reference
  end

  def self.haploid_mutation(mutation)
    chr, pos, alt = mutation.split(":")
    if chr =~ /^copy-\d/
      [chr, pos, alt] * ":"
    else
      chr = chr.sub(/^chr/,'')
      num = Misc.digest([chr, pos] * ":").chars.inject(0){|acc,e| acc += e.hex }
      chr = %w(copy-1_chr copy-2_chr)[num % 2] + chr
      [chr, pos, alt] * ":"
    end
  end

  def self.haploid_SV(values)
    type, chr, start, eend, target_chr, target_start, target_end = values

    chr = chr.to_s.sub(/^chr/, '')
    target_chr = target_chr.to_s.sub(/^chr/, '') if target_chr && ! target_chr.empty?
    target_chr = target_chr.to_s

    if !chr.include?('copy-') 
      chr_copies = %w(copy-1_chr copy-2_chr)
      num = Misc.digest([chr, start, eend, type] * ":").chars.inject(0){|acc,e| acc += e.hex }
      chr = chr_copies[num % chr_copies.length] + chr
    end

    target_chr = chr if target_chr == 'same' || target_chr == 'cis'

    if target_chr && ! target_chr.empty? && ! target_chr.include?('copy-')
      chr_copies = %w(copy-1_chr copy-2_chr)
      num = Misc.digest([chr, start, eend, target_chr, type] * ":").chars.inject(0){|acc,e| acc += e.hex }
      target_chr = chr_copies[num % chr_copies.length] + target_chr
    end

    [type, chr, start, eend, target_chr, target_start, target_end]
  end

  input :mutations, :array, "List of mutations to partially cover"
  input :min, :integer, "Minimum number of variants per chromosome to cover", 100
  task :miniref_sizes => :tsv do |mutations,min|
    sizes = HTSBenchmark.calculate_sizes(mutations, min)
    TSV.setup(sizes, :key_field => "Chromosome", :fields => ["Size"], :type => :single, :cast => :to_i)
  end

  input :organism, :string, "Organism code, no build", "Hsa"
  input :reference, :string, "Reference code", "hg38", :jobname => true
  input :padding, :integer, "Extra bases to add to reference", 1_000
  input :do_vcf, :boolean, "Minimize also the vcfs", false
  input :sizes, :tsv, "Sizes of each chromosome's beggining to preserve"
  extension 'fa.gz'
  task :miniref => :binary do |organism,reference,padding,do_vcf,sizes|

    sizes[:padding] = padding

    output = file(reference)

    reference_path = Rbbt.share.organisms[organism][reference]
    files = reference_path.glob_all("**/*")

    files.each do |file|
      subpath = file.original.sub(reference_path, '')

      target = output[subpath].find.remove_extension('.gz')
      type = case file
             when /\.vcf(?:\.gz)?$/
               next unless do_vcf
               HTSBenchmark.minify_vcf(file, target, sizes)

             when /\.fa(?:sta)?(?:\.gz)?$/
               HTSBenchmark.minify_fasta(file, target, sizes)
             else
               next
             end
      CMD.cmd("bgzip #{target}")
    end

    Open.link output["hg38.fa.gz"], self.tmp_path
    nil
  end

  dep :miniref
  extension 'fa.gz'
  task :miniref_ploidy => :binary do
    sout = Misc.open_pipe do |sin|
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          sin.puts ">copy-1_" + line[1..-1]
        else
          sin.puts line
        end
      end
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          sin.puts ">copy-2_" + line[1..-1]
        else
          sin.puts line
        end
      end
    end

    CMD.cmd("bgzip -c > #{self.tmp_path}", :in => sout)
    nil
  end

  input :mutations, :array, "Mutations to make haploid"
  input :reference, :binary, "Reference file", nil, :nofile => true
  task :mutations_to_reference =>  :tsv do |mutations,reference|
    reference = Samtools.prepare_FASTA(reference, file('reference'))
    mutation_reference_tsv = HTSBenchmark.mutation_reference(mutations, reference).to_s
  end

  input :organism, :string, "Organism code, no build", "Hsa"
  input :reference, :string, "Reference code", nil, :jobname => true
  input :do_vcf, :boolean, "Minimize also the vcfs", false
  input :bed_file, :file, "Bed file of ranges"
  extension 'fa.gz'
  task :sliceref => :binary do |organism,reference,do_vcf,bed_file|

    output = file(reference)

    reference_path = Rbbt.share.organisms[organism][reference]
    files = reference_path.glob_all("**/*")

    files.each do |file|
      subpath = file.original.sub(reference_path, '')

      target = output[subpath].find.remove_extension('.gz')
      type = case file
             when /\.vcf(?:\.gz)?$/
               next unless do_vcf
               HTSBenchmark.slice_vcf(file, target, bed_file)

             when /\.fa(?:sta)?(?:\.gz)?$/
               HTSBenchmark.slice_fasta(file, target, bed_file)
             else
               next
             end

      CMD.cmd("bgzip #{target}")
    end

    Open.link output["hg38.fa.gz"], self.tmp_path
    nil
  end
end

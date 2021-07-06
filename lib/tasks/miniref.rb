require 'NEAT_genreads'
require 'miniref'

module HTSBenchmark

  input :vcf, :file, "VCF to partially cover"
  input :min, :integer, "Minimum number of variants per chromosome to cover", 100
  dep Sequence, :genomic_mutations, :vcf_file => :vcf, :filters => [], :quality => nil
  task :miniref_sizes => :json do |vcf,min|
    require 'miniref'
    HTSBenchmark.calculate_sizes(step(:genomic_mutations), min)
  end

  input :sizes, :text, "Chromosome sizes in JSON", nil
  input :vcf_file, :file, "VCF to minimize"
  extension :vcf
  task :minify_vcf => :text do |sizes,vcf_file|
    sizes = Open.read(sizes) if Misc.is_filename?(sizes)
    sizes = JSON.parse(sizes) unless Hash === sizes

    HTSBenchmark.minify_vcf(vcf_file, self.tmp_path, sizes)
    nil
  end

  input :organism, :string, "Organism code, no build", "Hsa"
  input :reference, :string, "Reference code", "hg38", :jobname => true
  input :padding, :integer, "Extra bases to add to reference", 1_000
  dep :miniref_sizes
  extension 'fa.gz'
  task :miniref => :binary do |organism,reference,padding|
    require 'miniref'

    sizes = step(:miniref_sizes).load
    sizes[:padding] = padding

    output = file(reference)


    reference_path = Rbbt.share.organisms[organism][reference]
    files = reference_path.glob_all("**/*")
    files.each do |file|
      subpath = file.original.sub(reference_path, '')

      target = output[subpath].find.remove_extension('.gz')
      type = case file
             when /\.vcf(?:\.gz)?$/
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

  input :vcf, :file, "VCF to change poidy to", nil, :nofile => true
  input :haploid_reference, :file, "FASTA file with haploid reference", nil, :nofile => true
  extension :vcf
  task :vcf_ploidy => :binary do |vcf,haploid_reference|
    copies = {}
    TSV.traverse haploid_reference, :type => :array do |line|
      if m = line.match(/>((copy-\d+)_([^\s+]*))/)
        orig, copy, real = m.captures
        real = real.sub('chr', '')
        copies[real] ||= []
        copies[real] << orig
      end
    end

    io = TSV.traverse vcf, :type => :array, :into => :stream do |line|
      next line if line[0] == "#"
      parts = line.split("\t")
      chr, pos, _id, ref, alt   = parts

      if chr.include? 'copy-'
        copy = chr
      else
        chr = chr.sub('chr', '')
        chr_copies = copies[chr]
        num = Misc.digest([chr, pos, ref, alt] * ":").chars.inject(0){|acc,e| acc += e.hex }
        copy = chr_copies[num % chr_copies.length]
      end

      ([copy] + parts[1..-1]) * "\t"
    end

    Misc.consume_stream(io, false, self.tmp_path)
    nil
  end
end

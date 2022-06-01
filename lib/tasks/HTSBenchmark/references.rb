require 'HTSBenchmark/miniref'
require 'HTSBenchmark/sliceref'
require 'HTSBenchmark/haploid'

module HTSBenchmark
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
    reference_path.produce
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
  input :haploid_chromosomes, :array, "Chromosomes to not duplicate", %w(M X Y)
  extension 'fa.gz'
  task :miniref_ploidy => :binary do |haploid_chromosomes|
    haploid_chromosomes = haploid_chromosomes.collect{|c| c.sub('chr', '') }
    sout = Misc.open_pipe do |sin|
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          sin.puts ">copy-1_" + line[1..-1]
        else
          sin.puts line
        end
      end

      skip = false
      TSV.traverse step(:miniref), :type => :array do |line|
        if line =~ />/
          chr = line[1..-1].split(/\s+/).first.sub('chr', '')
          if haploid_chromosomes.include?(chr)
            skip = true
          else
            skip = false
          end
          sin.puts ">copy-2_" + line[1..-1] unless skip
        else
          sin.puts line unless skip
        end
      end
    end

    CMD.cmd("bgzip -c > #{self.tmp_path}", :in => sout)
    nil
  end

  input :mutations, :array, "Mutations to make haploid"
  input :reference, :binary, "Reference file", nil, :nofile => true
  task :mutations_to_reference =>  :tsv do |mutations,reference|
    reference = reference.path if Step === reference
    reference = Samtools.prepare_FASTA(reference, file('reference'))
    mutation_reference_tsv = HTSBenchmark.mutation_reference(mutations, reference).to_s
  end

  input :organism, :string, "Organism code, no build", "Hsa"
  input :reference, :string, "Reference code", nil, :jobname => true
  input :do_vcf, :boolean, "Minimize also the vcfs", false
  input :bed_file, :file, "Bed file of ranges", nil, :nofile => true
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

  input :vcf_file, :file, "VCF file to restore"
  input :bed_file, :file, "Bed file of ranges"
  extension :vcf
  task :restore_sliced_vcf => :text do |vcf,bed|
    HTSBenchmark.restore_sliced_vcf(vcf, self.tmp_path, bed)
    nil
  end
end

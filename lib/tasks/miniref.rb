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
  dep :miniref_sizes, :compute => :produce
  task :miniref => :text do |organism,reference,padding|
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

end

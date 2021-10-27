module HTSBenchmark

  dep :simulate_tumor_normal_diploid, :do_vcf => true
  task :bundle_tumor_normal => :array do
    samples = step(:simulate_tumor_normal_diploid)

    neat = step(:NEAT_simulate_DNA)

    fastq_dir = file('FASTQ')
    Dir.glob(samples.files_dir + "/*.gz").each do |file|
      Open.link file, fastq_dir[File.basename(file)]
    end

    ref_dir = file('reference')
    miniref = step(:miniref)
    miniref.files.each do |file|
      Open.link File.join(miniref.files_dir, file), ref_dir[File.basename(file)] 
    end

    BWA.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)
    Samtools.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)
    GATK.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir)
    end

    Dir.glob(self.files_dir + "**/*")
  end

  dep :population, :bundle => true
  task :bundle_population => :array do
    samples = step(:population)

    fastq_dir = file('FASTQ')
    Dir.glob(samples.files_dir + "/*.gz").each do |file|
      Open.link file, fastq_dir[File.basename(file)]
    end

    ref_dir = file('reference')
    miniref = step(:miniref)
    miniref.files.each do |file|
      Open.link File.join(miniref.files_dir, file), ref_dir[File.basename(file)] 
    end

    BWA.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)
    Samtools.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)
    GATK.prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir)
    end

    Dir.glob(self.files_dir + "**/*")
  end
end

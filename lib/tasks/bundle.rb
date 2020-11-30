module HTSBenchmark

  dep :simulated_sample
  dep_task :aligned_BAM, HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder, :reference => 'hg38', :sequencing_center => "Simulation", :platform => "NEAT_GenReads" do |jobname,options,dependencies|
    neat = dependencies.flatten.select{|dep| dep.task_name.to_s == "simulated_sample"}.first
    options[:fastq1] = neat.file('output/tumor_read1.fq.gz')
    options[:fastq2] = neat.file('output/tumor_read2.fq.gz')
    {:inputs => options, :jobname => jobname}
  end

  dep :simulated_sample
  dep_task :aligned_BAM_normal, HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder, :reference => 'hg38', :sequencing_center => "Simulation", :platform => "NEAT_GenReads" do |jobname,options,dependencies|
    neat = dependencies.flatten.select{|dep| dep.task_name.to_s == "simulated_sample"}.first
    options[:fastq1] = neat.file('output/normal_read1.fq.gz')
    options[:fastq2] = neat.file('output/normal_read2.fq.gz')
    {:inputs => options, :jobname => jobname}
  end

  dep :aligned_BAM
  dep :aligned_BAM_normal
  task :bundle => :text do
    neat = step(:NEAT_genreads)
    miniref = step(:miniref)

    bam = step(:aligned_BAM).path
    bam_normal = step(:aligned_BAM_normal).path

    bundle = file('bundle')

    germline = neat.file('output/normal_golden.vcf')
    tumor = neat.file('output/tumor_golden.vcf')

    #somatic = self.rec_dependencies.select{|dep| dep.task_name.to_s == 'minify_vcf'}.collect{|dep| dep.path}

    reference = nil
    vcf_files = []

    orig_reference =  miniref.file('hg38/hg38.fa.gz')
    reference = BWA.prepare_FASTA orig_reference
    reference = GATK.prepare_FASTA orig_reference
    reference = Samtools.prepare_FASTA orig_reference

    miniref.file('hg38').known_sites.glob("*.vcf.gz").each do |file|
      vcf_files << GATK.prepare_VCF(file)
    end

    golden_bam_normal = neat.file('output/normal_golden.bam')
    golden_bam = neat.file('output/tumor_golden.bam')

    Open.ln_s File.dirname(reference), bundle.reference

    Open.mkdir bundle.known_sites
    vcf_files.each do |vcf|
      Dir.glob(File.dirname(vcf) + '/*').each do |file|
        Open.ln_s file, bundle.known_sites[File.basename(file)]
      end
    end

    neat.file("output").glob("tumor_read*gz").each do |file|
      Open.ln_s file, bundle.FASTQ.tumor[File.basename(file)]
    end

    neat.file("output").glob("normal_read*gz").each do |file|
      Open.ln_s file, bundle.FASTQ.normal[File.basename(file)]
    end

    Open.mkdir bundle.BAM
    Open.mkdir bundle.golden_BAM

    [golden_bam, golden_bam_normal].each do |file|
      indexed = Samtools.prepare_BAM(file)
      Dir.glob(indexed + "*").each do |file|
        Open.ln_s file, bundle.golden_BAM
      end
    end

    indexed = Samtools.prepare_BAM(bam)
    Dir.glob(indexed + "*").each do |file|
      name = File.basename(file).sub(/.*\./, 'tumor.')
      Open.ln_s file, bundle.BAM[name]
    end

    indexed = Samtools.prepare_BAM(bam_normal)
    Dir.glob(indexed + "*").each do |file|
      name = File.basename(file).sub(/.*\./, 'normal.')
      Open.ln_s file, bundle.BAM[name]
    end

    Open.mkdir bundle.truth
    Open.ln_s germline, bundle.truth["germline.vcf"]
    Open.ln_s tumor, bundle.truth["tumor.vcf"]
    CMD.cmd_log('vcftools', "--vcf  #{tumor} --recode --exclude-positions #{germline} --stdout > #{bundle.truth["somatic.vcf"]}")

    "DONE"
  end

  #dep :benchmark
  #task :bundle => :text do
  #  benchmark = step(:benchmark)
  #  neat = step(:NEAT_genreads)

  #  germline, somatic = self.rec_dependencies.select{|dep| dep.task_name.to_s == 'minify_vcf'}.collect{|dep| dep.path}

  #  bundle = file('bundle')
  #  stage = benchmark.file('stage')

  #  bam = bam_normal = reference = tumor = normal = nil
  #  vcf_files = []
  #  Misc.in_dir stage do

  #    bam = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 BAM -ck HTS_high -pf", :log => true, :pipe => true).read.strip
  #    bam_normal = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 BAM_normal -ck HTS_high -pf", :log => true, :pipe => true).read.strip

  #    orig_reference =  stage.share.organisms.Hsa.hg38["hg38.fa.gz"]
  #    reference = BWA.prepare_FASTA orig_reference
  #    reference = GATK.prepare_FASTA orig_reference
  #    reference = Samtools.prepare_FASTA orig_reference

  #    stage.share.organisms.Hsa.hg38.known_sites.glob("*.vcf.gz").each do |file|
  #      vcf_files << GATK.prepare_VCF(file)
  #    end
  #    vcf_files
  #  end

  #  golden_bam_normal = neat.file('output/normal_golden.bam')
  #  golden_bam = neat.file('output/tumor_golden.bam')

  #  Open.ln_s File.dirname(reference), bundle.reference

  #  Open.mkdir bundle.known_sites
  #  vcf_files.each do |vcf|
  #    Dir.glob(File.dirname(vcf) + '/*').each do |file|
  #      Open.ln_s file, bundle.known_sites[File.basename(file)]
  #    end
  #  end

  #  stage.share.data.studies.Miniref.glob("*/Miniref/*.gz").each do |file|
  #    Open.ln_s file, bundle.FASTQ.tumor[File.basename(file)]
  #  end

  #  stage.share.data.studies.Miniref.glob("*/Miniref_normal/*.gz").each do |file|
  #    Open.ln_s file, bundle.FASTQ.normal[File.basename(file)]
  #  end

  #  Open.mkdir bundle.BAM
  #  Open.mkdir bundle.golden_BAM

  #  [golden_bam, golden_bam_normal].each do |file|
  #    indexed = Samtools.prepare_BAM(file)
  #    Dir.glob(indexed + "*").each do |file|
  #      Open.ln_s file, bundle.golden_BAM
  #    end
  #  end

  #  indexed = Samtools.prepare_BAM(bam)
  #  Dir.glob(indexed + "*").each do |file|
  #    name = File.basename(file).sub('Miniref', 'tumor')
  #    Open.ln_s file, bundle.BAM[name]
  #  end

  #  indexed = Samtools.prepare_BAM(bam_normal)
  #  Dir.glob(indexed + "*").each do |file|
  #    name = File.basename(file).sub('Miniref', 'normal')
  #    Open.ln_s file, bundle.BAM[name]
  #  end

  #  Open.mkdir bundle.truth
  #  Open.ln_s germline, bundle.truth["germline.vcf"]
  #  Open.ln_s somatic, bundle.truth["somatic.vcf"]

  #  "DONE"
  #end
end

module HTSBenchmark

  input :somatic, :array, "Somatic variants"
  input :germline, :array, "Germline variants"
  task :merge_somatic_germline => :array do |somatic,germline|
    somatic = Open.read(somatic).split("\n") if Misc.is_filename?(somatic)
    germline = Open.read(germline).split("\n") if Misc.is_filename?(germline)
    somatic + germline
  end

  dep :genotype_germline_hg38, :compute => :produce
  dep :genotype_somatic_hg38, :compute => :produce
  dep :miniref_sizes, :mutations => :genotype_somatic_hg38
  dep :miniref_ploidy, :sizes => :miniref_sizes, :jobname => 'hg38', :organism => "Hsa"
  dep :merge_somatic_germline, :somatic => :genotype_somatic_hg38, :germline => :genotype_germline_hg38
  dep_task :simulate_tumor, HTSBenchmark, :NEAT_simulate_DNA, :reference => :miniref_ploidy, :mutations => :merge_somatic_germline, :haploid_reference => true, :depth => 60 do |jobname,options|
    {:inputs => options, :jobname => jobname + "_tumor"}
  end

  dep :genotype_germline_hg38, :compute => :produce
  dep :genotype_somatic_hg38, :compute => :produce
  dep :miniref_sizes, :mutations => :genotype_somatic_hg38
  dep :miniref_ploidy, :sizes => :miniref_sizes, :jobname => 'hg38', :organism => "Hsa"
  dep_task :simulate_normal, HTSBenchmark, :NEAT_simulate_DNA, :reference => :miniref_ploidy, :mutations => :genotype_germline_hg38, :haploid_reference => true, :depth => 30 do |jobname,options|
    {:inputs => options, :jobname => jobname + "_normal"}
  end

  dep :simulate_tumor
  dep :simulate_normal
  task :simulate_tumor_normal => :array do 

    normal, tumor = dependencies 

    nr1, nr2 = normal.load
    tr1, tr2 = tumor.load

    Open.link(nr1, file('normal_read1.fq.gz'))
    Open.link(nr2, file('normal_read2.fq.gz'))
    Open.link(tr1, file('tumor_read1.fq.gz'))
    Open.link(tr1, file('tumor_read2.fq.gz'))

    Dir.glob(files_dir + "/*.fq.gz")
  end

  dep :genotype_germline_hg38, :compute => :produce
  dep :genotype_somatic_hg38, :compute => :produce
  dep :miniref_sizes, :mutations => :genotype_somatic_hg38
  dep :miniref, :sizes => :miniref_sizes, :jobname => 'hg38', :organism => "Hsa"
  dep :merge_somatic_germline, :somatic => :genotype_somatic_hg38, :germline => :genotype_germline_hg38
  dep_task :simulate_tumor_diploid, HTSBenchmark, :NEAT_simulate_DNA, :reference => :miniref, :mutations => :merge_somatic_germline, :haploid_reference => false, :depth => 60 do |jobname,options|
    {:inputs => options, :jobname => jobname + "_tumor"}
  end

  dep :genotype_germline_hg38, :compute => :produce
  dep :genotype_somatic_hg38, :compute => :produce
  dep :miniref_sizes, :mutations => :genotype_somatic_hg38
  dep :miniref, :sizes => :miniref_sizes, :jobname => 'hg38', :organism => "Hsa"
  dep_task :simulate_normal_diploid, HTSBenchmark, :NEAT_simulate_DNA, :reference => :miniref, :mutations => :genotype_germline_hg38, :haploid_reference => false, :depth => 30 do |jobname,options|
    {:inputs => options, :jobname => jobname + "_normal"}
  end

  dep :simulate_tumor_diploid
  dep :simulate_normal_diploid
  task :simulate_tumor_normal_diploid => :array do 

    normal, tumor = dependencies 

    nr1, nr2 = normal.load
    tr1, tr2 = tumor.load

    Open.link(nr1, file('normal_read1.fq.gz'))
    Open.link(nr2, file('normal_read2.fq.gz'))
    Open.link(tr1, file('tumor_read1.fq.gz'))
    Open.link(tr1, file('tumor_read2.fq.gz'))

    normal_bam = normal.file('output').glob("*.bam").first
    tumor_bam = normal.file('output').glob("*.bam").first

    Open.link(normal_bam, file('normal.bam'))
    Open.link(tumor_bam, file('tumor.bam'))

    Dir.glob(files_dir + "/*.fq.gz")
  end
end

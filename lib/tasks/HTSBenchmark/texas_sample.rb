module HTSBenchmark

  TEXAS_SAMPLE = "DO262483"
  TEXAS_SAMPLE_PREFIX_TUMOR = "TCRB-CA.DO262483.SA622771.wgs"
  TEXAS_SAMPLE_PREFIX_NORMAL = "TCRB-CA.DO262483.SA622770.wgs"

  input :regions_to_slice, :file, "Regions to slice", nil, :required => true
  dep Sample, :BAM, :jobname => TEXAS_SAMPLE, :reference => "hg38"
  dep Sample, :BAM_normal, :jobname => TEXAS_SAMPLE, :reference => "hg38"
  dep HTS, :extract_BAM_region_with_mates_samtools, :bam => :BAM, :bed_file => :regions_to_slice
  dep HTS, :extract_BAM_region_with_mates_samtools, :bam => :BAM_normal, :bed_file => :regions_to_slice
  dep HTS, :revert_BAM, :bam_file => :placeholder, :by_group => false do |jobname,options,dependencies|
    slice = dependencies.flatten.select{|d| d.task_name.to_s == "extract_BAM_region_with_mates_samtools" }.first
    {:inputs => options.merge({:bam_file => slice })}
  end
  dep HTS, :revert_BAM, :bam_file => :placeholder, :by_group => false do |jobname,options,dependencies|
    slice = dependencies.flatten.select{|d| d.task_name.to_s == "extract_BAM_region_with_mates_samtools" }.last
    {:inputs => options.merge({:bam_file => slice })}
  end
  dep HTSBenchmark, :sliceref, :bed_file => :regions_to_slice, :reference => 'hg38', :do_vcf => true
  task :texas_benchmark => :text do
    bam, bam_normal, fastq, fastq_normal, ref = dependencies

    bam_dir = file('WGS').BAM
    Open.link bam.path, bam_dir.tumor["#{TEXAS_SAMPLE}.bam"]
    Open.link bam_normal.path, bam_dir.normal["#{TEXAS_SAMPLE}_normal.bam"]

    orig_reference = HTS.helpers[:reference_file].call('hg38')
    reference = GATK.prepare_FASTA orig_reference

    cram_dir = file('WGS').CRAM
    Open.mkdir cram_dir.tumor
    Open.mkdir cram_dir.normal
    CMD.cmd(:samtools, "view -T #{reference} -C -o #{cram_dir.tumor["#{TEXAS_SAMPLE}.cram"]} #{ bam_dir.tumor["#{TEXAS_SAMPLE}.bam"]} ")
    CMD.cmd(:samtools, "view -T #{reference} -C -o #{cram_dir.normal["#{TEXAS_SAMPLE}_normal.cram"]} #{ bam_dir.normal["#{TEXAS_SAMPLE}_normal.bam"]} ")

    fastq_dir = file('WGS').FASTQ
    Open.mkdir fastq_dir.tumor
    Open.mkdir fastq_dir.normal
    CMD.cmd(:samtools, "fastq -1 #{fastq_dir.tumor["tumor_read1.fq.gz"]} -2 #{fastq_dir.tumor["tumor_read2.fq.gz"]} -0 /dev/null -s /dev/null #{fastq.path}")
    CMD.cmd(:samtools, "fastq -1 #{fastq_dir.normal["normal_read1.fq.gz"]} -2 #{fastq_dir.normal["normal_read2.fq.gz"]} -0 /dev/null -s /dev/null #{fastq_normal.path}")

    Open.link ref.file('hg38'), File.join(files_dir, 'mini-reference')
    Dir.glob(files_dir + "/*") * "\n"
  end

  dep :texas_benchmark
  input :file_prefix, :string, "Prefix for files"
  task :bundle_texas => :text do |prefix|
    bench = step(:texas_benchmark)

    Open.write(file("README.md"), Rbbt.doc["README.real.md"].read)

    date = Time.now.strftime "%Y%m%d"
    %w(normal tumor).each do |type|
      sample_prefix = type == "normal" ? TEXAS_SAMPLE_PREFIX_NORMAL : TEXAS_SAMPLE_PREFIX_TUMOR
      bench.file("WGS").BAM[type].glob("*.bam").each do |file|
        new_file = [prefix || sample_prefix, clean_name, type, date] * "." + '.bam'
        Open.cp file, file('WGS').BAM[type][new_file]
      end
      bench.file("WGS").CRAM[type].glob("*.cram").each do |file|
        new_file = [prefix || sample_prefix, clean_name, type, date] * "." + '.cram'
        Open.cp file, file('WGS').CRAM[type][new_file]
      end
      bench.file("WGS").FASTQ[type].glob("*_read1.fq.gz").each do |file|
        new_file = [prefix || sample_prefix, clean_name, type, date] * "." + '_read1.fq.gz'
        Open.cp file, file('WGS').FASTQ[type][new_file]
      end
      bench.file("WGS").FASTQ[type].glob("*_read2.fq.gz").each do |file|
        new_file = [prefix || sample_prefix, clean_name, type, date] * "." + '_read2.fq.gz'
        Open.cp file, file('WGS').FASTQ[type][new_file]
      end
    end

    Open.cp bench.file('reference'), file('mini-reference')

    file('WGS').CRAM.glob('*/*.cram').each do |f|
      CMD.cmd(:samtools, "index #{f}")
    end

    file('WGS').BAM.glob('*/*.bam').each do |f|
      CMD.cmd(:samtools, "index #{f}")
    end

    reference = file('mini-reference')["hg38.fa.gz"]

    GATK.prepare_FASTA reference, file('mini-reference')
    BWA.prepare_FASTA reference, file('mini-reference')
    Samtools.prepare_FASTA reference, file('mini-reference')

    file('mini-reference').known_sites.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF vcf, file('mini-reference').known_sites
    end

    %w(tumor normal).each do |type|
      HTS.workdir = file('tmp')

      f1 = file('WGS').FASTQ[type].glob("*_read1.fq.gz").first
      f2 = file('WGS').FASTQ[type].glob("*_read2.fq.gz").first

      job = HTS.job(:BAM, nil,
              :skip_rescore => true,
              :fastq1 => f1,
              :fastq2 => f2,
              :reference => reference)

      job.produce

      sample_prefix = type == "normal" ? TEXAS_SAMPLE_PREFIX_NORMAL : TEXAS_SAMPLE_PREFIX_TUMOR
      new_file = [prefix || sample_prefix, clean_name, type, date, "mini-aln"] * "."

      if type == 'tumor'
        target_bam = file('WGS')["BAM-mini"][type][new_file] + '.bam'
        target = file('WGS')["CRAM-mini"][type][new_file] + '.cram'
      else
        target_bam = file('WGS')["BAM-mini"][type][new_file] + '.bam'
        target = file('WGS')["CRAM-mini"][type][new_file] + '.cram'
      end

      Open.mkdir File.dirname(target)

      Open.link job.path, target_bam
      CMD.cmd(:samtools, "index #{target_bam}")

      CMD.cmd(:samtools, "view -T #{reference} -C -o #{target} #{ job.path } ")
      CMD.cmd(:samtools, "index #{target}")
    end

    Open.cp recursive_inputs[:regions_to_slice], file('inputs')["regions.bed"]

    Open.rm_rf file('tmp')

    "DONE"
  end
end

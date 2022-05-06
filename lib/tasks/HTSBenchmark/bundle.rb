module HTSBenchmark

  helper :prepare_FASTA do |file, target|
    BWA.prepare_FASTA(file, target)
    Samtools.prepare_FASTA(file, target)
    file = GATK.prepare_FASTA(file, target)

    dir = File.dirname(file)

    alt = file + '.alt'
    if file.include?("hg38") && ! File.exists?(alt) 
      Open.write(alt, '')
    end

    Dir.glob(File.join(dir, '*.gz.*')).each do |source|
      target = source.sub('.gz','')
      next if File.exists?(target)
      Open.link source, target 
    end
    
    intervals = file + '.byNS.interval_list'
    args = {}
    args["REFERENCE"] = file
    args["OUTPUT"] = intervals
    args["MAX_TO_MERGE"] = GATKShard::GAP_SIZE
    GATK.run_log("ScatterIntervalsByNs", args)

    header = []
    chrs_intervals = {}
    Open.read(intervals).split("\n").each do |line|
      if line =~ /^@/
        header << line
      else
        next if line =~ /Nmer$/
        chr = line.split(/\s/).first
        chrs_intervals[chr] ||= []
        chrs_intervals[chr] << line
      end
    end

    dir = Path.setup(file + '.interval_lists')
    chrs_intervals.each do |chr,lines|
      ifile = dir[chr + '.byNS.interval_list']
      Open.write(ifile, (header + lines) * "\n" + "\n")
    end

    uncompressed = file.sub('.gz', '')

    if ! File.exists?(uncompressed)
      CMD.cmd("zcat '#{file}' > #{uncompressed}")
    end

    file
  end

  dep :simulate_tumor_normal_diploid, :do_vcf => true
  task :bundle_tumor_normal => :array do
    samples = step(:simulate_tumor_normal_diploid)

    neat_tumor, neat_normal = rec_dependencies.select{|d| d.task_name == :NEAT_simulate_DNA }

    log :truth, "Preparing truth set"

    truth_dir = file('truth')
    Open.mkdir truth_dir

    Open.link neat_tumor.file('output').glob("*.vcf.gz").first, truth_dir
    Open.link neat_normal.file('output').glob("*.vcf.gz").first, truth_dir

    truth_dir.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, truth_dir)
    end

    CMD.cmd(:bcftools, "isec -C '#{truth_dir.glob("*tumor.vcf.gz").first}' '#{truth_dir.glob("*normal.vcf.gz").first}' -p #{truth_dir}")

    truth_dir.glob("*.txt").each{|f| Open.rm f }

    truth_dir.glob("*.vcf").each do |f| 
      CMD.cmd("sed -i 's/WP=0$/WP=0\\/1/;s/WP=1$/WP=1\\/0/' '#{f}'")
      CMD.cmd("bgzip '#{f}'")
      Open.mv f + '.gz', truth_dir['somatic.vcf.gz']
      GATK.prepare_VCF(truth_dir["somatic.vcf.gz"], truth_dir) 
    end

    log :FASTQ, "Preparing FASTQ files"

    prefix_normal = nil
    prefix_tumor = nil
    date = Time.now.strftime "%Y%m%d"
    %w(normal tumor).each do |type|
      neat_job = type == 'tumor' ? neat_tumor : neat_normal

      neat_job.file('output').glob("*_read1.fq.gz").each do |file|
        new_file = [clean_name, type, date] * "." + '_read1.fq.gz'
        Open.cp file, file('WGS').FASTQ[type][new_file]
      end

      neat_job.file('output').glob("*_read2.fq.gz").each do |file|
        new_file = [clean_name, type, date] * "." + '_read2.fq.gz'
        Open.cp file, file('WGS').FASTQ[type][new_file]
      end

      neat_job.file('output').glob("*.bam").each do |file|
        new_file = [clean_name, type, date] * "." + '.bam'
        Open.cp file, file('WGS').BAM[type][new_file]
      end
    end

    log :reference, "Preparing reference"

    ref_dir = file('reference')
    miniref = step(:miniref)
    miniref.files_dir.glob("*/*.fa*").each do |file|
      Open.link file, ref_dir[File.basename(file)] 
    end

    miniref.files_dir.glob("*/known_sites/*").each do |file|
      Open.link file, ref_dir.known_sites[File.basename(file)] 
    end

    prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.known_sites.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir.known_sites)
    end

    Dir.glob(self.files_dir + "**/*")
  end

  dep :population, :bundle => true
  task :bundle_population => :array do
    samples = step(:population)

    log :FASTQ, "Preparing FASTQ files"

    fastq_dir = file('WGS/FASTQ')
    Dir.glob(samples.files_dir + "/tumor_*.gz").each do |file|
      Open.link file, fastq_dir.tumor[File.basename(file)]
    end

    Dir.glob(samples.files_dir + "/normal_*.gz").each do |file|
      Open.link file, fastq_dir.normal[File.basename(file)]
    end

    log :truth_set, "Preparing truth set"

    Open.write file('truth/somatic_mutations.list'), step(:population_genotypes).file('total_mutations').read.gsub(/copy-\d+_/,'')
    Open.link Sequence.job(:mutations_to_vcf, nil, :organism => "Hsa/may2017", :positions => file('truth/somatic_mutations.list')).recursive_clean.produce.path, file('truth/somatic_mutations_all.vcf')
    HTSBenchmark.minify_vcf file('truth/somatic_mutations_all.vcf'), file('truth/somatic.vcf'), step(:miniref_sizes).load
    CMD.cmd_log('gzip', file('truth/somatic.vcf'))

    Open.write file('truth/evolution.yaml'), TSV.get_stream(recursive_inputs[:evolution])

    log :reference, "Preparing reference"

    ref_dir = file('reference')
    miniref = step(:miniref)
    miniref.files_dir.glob("*/*.fa*").each do |file|
      Open.link file, ref_dir[File.basename(file)] 
    end

    miniref.files_dir.glob("*/known_sites/*").each do |file|
      Open.link file, ref_dir.known_sites[File.basename(file)] 
    end

    prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.known_sites.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir.known_sites)
    end

    Dir.glob(self.files_dir + "**/*")
  end

  dep :simulate_population, :bundle => true
  dep_task :bundle_simulate_population, HTSBenchmark, :bundle_population, :evolution => :placeholder do |jobname,options,dependencies|
    simulate_population = dependencies.flatten.first
    {:inputs => options.merge("HTSBenchmark#population" => simulate_population, :not_overriden => true), :jobname => jobname}
  end
end
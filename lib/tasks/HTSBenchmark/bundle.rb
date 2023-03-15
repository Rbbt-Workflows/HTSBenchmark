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

  dep :contaminated_population
  task :bundle_population => :array do
    samples = step(:contaminated_population)

    log :FASTQ, "Preparing FASTQ files"

    fastq_dir = file('WGS/FASTQ')
    Dir.glob(samples.files_dir + "/tumor_*.gz").each do |file|
      Open.link file, fastq_dir.tumor[File.basename(file)]
    end

    Dir.glob(samples.files_dir + "/normal_*.gz").each do |file|
      Open.link file, fastq_dir.normal[File.basename(file)]
    end

    log :truth_set, "Preparing truth set"

    chr_sizes = step(:miniref_sizes).load
    Open.open(file('inputs/regions.bed'), :mode => 'w') do |bed|
      chr_sizes.each do |chr,size|
        chr = "chr" + chr unless chr.include? 'chr'
        bed.puts [chr, 0, size ] * "\t" 
      end
    end

    inserted_mutations = []
    self.rec_dependencies.each do |dep|
      next unless dep.task_name.to_s == 'clone'
      vcf = dep.step("NEAT_simulate_DNA").file('output').glob("*.vcf.gz").first
      mutations = Sequence.job(:genomic_mutations, nil, :vcf_file => vcf).recursive_clean.run
      mutations = mutations.collect{|m| m.sub(/^copy-\d+_/,'') }
      inserted_mutations.concat mutations
    end


    somatic_mutations = step(:population_genotypes).file('total_mutations').read.gsub(/copy-\d+_/,'').split("\n")
    somatic_mutations = inserted_mutations & somatic_mutations

    HTSBenchmark.minify_mutations somatic_mutations, file('truth/somatic_mutations.list.all'), step(:miniref_sizes).load
    CMD.cmd("grep -v ':SV:' '#{file('truth/somatic_mutations.list.all')}' > '#{file('truth/somatic_mutations.list')}'", :nofail => true)
    vcf_job = Sequence.job(:mutations_to_vcf, nil, :organism => "Hsa/may2017", :positions => file('truth/somatic_mutations.list'))
    vcf_job.recursive_clean.produce
    CMD.cmd_log("gzip #{vcf_job.path } -c > #{file('truth/somatic.vcf.gz')}")

    HTSBenchmark.minify_mutations step(:genotype_germline_hg38).path, file('truth/germline_mutations.list'), step(:miniref_sizes).load
    vcf_job = Sequence.job(:mutations_to_vcf, nil, :organism => "Hsa/may2017", :positions => file('truth/germline_mutations.list'))
    vcf_job.recursive_clean.produce
    CMD.cmd_log("gzip #{vcf_job.path } -c > #{file('truth/germline.vcf.gz')}")

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

  #dep :simulate_population
  #dep_task :bundle_simulate_population_old, HTSBenchmark, :bundle_population, :evolution => :placeholder, :not_overriden => true do |jobname,options,dependencies|
  #  population = dependencies.flatten.first
  #  {:inputs => options.merge("HTSBenchmark#contaminated_population" => population), :jobname => jobname}
  #end

  dep :simulate_evolution
  dep_task :bundle_simulate_population, HTSBenchmark, :bundle_population, :evolution => :placeholder, :not_overriden => true do |jobname,options,dependencies|
    simevo = dependencies.flatten.first
    simevo.produce
    simevo.join unless simevo.done?
    {:inputs => options.merge(:evolution => Open.read(simevo.path)), :jobname => jobname}
  end


  dep :simulate_signature_population
  dep_task :bundle_simulate_signature_population, HTSBenchmark, :bundle_population, :evolution => :placeholder, :not_overriden => true do |jobname,options,dependencies|
    simevo = dependencies.flatten.first
    simevo.produce
    simevo.join unless simevo.done?
    {:inputs => options.merge(:evolution => Open.read(simevo.path)), :jobname => jobname}
  end
end

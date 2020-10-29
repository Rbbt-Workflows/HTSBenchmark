module HTSBenchmark
  dep :benchmark
  dep :somatic_variant_calling, :variant_caller => 'combined_caller_vcfs'
  task :bundle => :text do
    benchmark = step(:benchmark)

    germline, somatic = self.rec_dependencies.select{|dep| dep.task_name.to_s == 'minify_vcf'}.collect{|dep| dep.path}

    bundle = file('bundle')
    stage = benchmark.file('stage')

    bam = bam_normal = reference = tumor = normal = nil
    vcf_files = []
    Misc.in_dir stage do

      bam = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 BAM -ck HTS_high -pf", :log => true, :pipe => true).read.strip
      bam_normal = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 BAM_normal -ck HTS_high -pf", :log => true, :pipe => true).read.strip

      orig_reference =  stage.share.organisms.Hsa.hg38["hg38.fa.gz"]
      reference = BWA.prepare_FASTA orig_reference
      reference = GATK.prepare_FASTA orig_reference
      reference = Samtools.prepare_FASTA orig_reference

      stage.share.organisms.Hsa.hg38.known_sites.glob("*.vcf.gz").each do |file|
        vcf_files << GATK.prepare_VCF(file)
      end
      vcf_files
    end

    Open.ln_s File.dirname(reference), bundle.reference

    Open.mkdir bundle.known_sites
    vcf_files.each do |vcf|
      Dir.glob(File.dirname(vcf) + '/*').each do |file|
        Open.ln_s file, bundle.known_sites[File.basename(file)]
      end
    end

    stage.share.data.studies.Miniref.glob("*/Miniref/*.gz").each do |file|
      Open.ln_s file, bundle.FASTQ.tumor[File.basename(file)]
    end

    stage.share.data.studies.Miniref.glob("*/Miniref_normal/*.gz").each do |file|
      Open.ln_s file, bundle.FASTQ.normal[File.basename(file)]
    end

    Open.mkdir bundle.BAM
    Open.ln_s bam, bundle.BAM["tumor.bam"]
    Open.ln_s bam_normal, bundle.BAM["normal.bam"]

    Open.mkdir bundle.truth
    Open.ln_s germline, bundle.truth["germline.vcf"]
    Open.ln_s somatic, bundle.truth["somatic.vcf"]

    "DONE"
  end
end

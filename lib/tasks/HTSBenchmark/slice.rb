module Sample
  dep :BAM
  dep_task :slice_BAM, HTS, :extract_BAM_region_with_mates_samtools, :bam => :BAM

  dep :BAM
  dep_task :slice_BAM_normal, HTS, :extract_BAM_region_with_mates_samtools, :bam => :BAM_normal
end

module HTSBenchmark

  input :regions_to_slice, :file, "Regions to slice", nil, :required => true, :nofile => true
  dep Sample, :BAM
  dep Sample, :BAM_normal
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
  task :slice_bundle => :text do |regions_to_slice|

    fastq_tumor, fastq_normal = rec_dependencies.select{|d| d.task_name == :revert_BAM }

    log :FASTQ, "Preparing FASTQ files"

    date = Time.now.strftime "%Y%m%d"

    Open.cp fastq_tumor.path, file('WGS').BAM.tumor[[clean_name, 'tumor', date] * "." + '.ubam']
    Open.cp fastq_normal.path, file('WGS').BAM.normal[[clean_name, 'normal', date] * "." + '.ubam']

    log :reference, "Preparing reference"

    ref_dir = file('reference')
    slicedref = step(:sliceref)
    slicedref.files_dir.glob("*/*.fa*").each do |file|
      Open.link file, ref_dir[File.basename(file)] 
    end

    slicedref.files_dir.glob("*/known_sites/*").each do |file|
      Open.link file, ref_dir.known_sites[File.basename(file)] 
    end

    prepare_FASTA(ref_dir.glob("*.fa.gz").first, ref_dir)

    ref_dir.known_sites.glob("*.vcf.gz").each do |vcf|
      GATK.prepare_VCF(vcf, ref_dir.known_sites)
    end

    Open.cp regions_to_slice, file('inputs')["regions.bed"]

    Dir.glob(self.files_dir + "**/*")
  end
end

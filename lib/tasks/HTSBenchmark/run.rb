require 'rbbt/workflow/util/orchestrator'

module HTSBenchmark

  helper :find_bundle do |bundle|
    if Step === bundle
      bundle = bundle.files_dir 
    else
      bundle = bundle + '.files' if File.exists?(bundle + '.files')
    end

    if File.exists?(bundle) || bundle.include?("/")
      bundle =  Path.setup(bundle)
    else
      bundle = Rbbt.share.data.HTSBenchmark[bundle].find
    end
    
    bundle
  end

  helper :bundle_study_and_sample do |bundle|
    bundle = HTSBenchmark.helpers[:find_bundle].call bundle
    study_name = "HTSB-St-#{Misc.digest(bundle)[0..5]}"
    sample_name = "HTSB-Sa-#{Misc.digest(bundle)[0..5]}"

    [study_name, sample_name]
  end

  input :bundle, :file, "Bundle directory", nil, :nofile => true, :required => true
  input :input_type, :select, "Input type FASTQ, CRAM, BAM or realign", "FASTQ", :select_options => %w(FASTQ BAM CRAM BAM-re CRAM-re)
  task :stage_bundle => :string do |bundle,input_type|
    work = file('work')

    bundle = find_bundle bundle
    study_name, sample_name = bundle_study_and_sample bundle

    Open.link bundle.truth, work.truth
    Open.link bundle.inputs, work.inputs
    Open.link bundle.reference, work.reference

    case input_type.to_s
    when "FASTQ"
      Open.link bundle.WGS.FASTQ.tumor.glob("*_read1.fq.gz").first, work.share.data.studies[study_name].WGS[sample_name]["tumor_read1.fq.gz"]
      Open.link bundle.WGS.FASTQ.tumor.glob("*_read2.fq.gz").first, work.share.data.studies[study_name].WGS[sample_name]["tumor_read2.fq.gz"]
      Open.link bundle.WGS.FASTQ.normal.glob("*_read1.fq.gz").first, work.share.data.studies[study_name].WGS[sample_name + '_normal']["normal_read1.fq.gz"]
      Open.link bundle.WGS.FASTQ.normal.glob("*_read2.fq.gz").first, work.share.data.studies[study_name].WGS[sample_name + '_normal']["normal_read2.fq.gz"]
    when "BAM"
      Open.link bundle.WGS.BAM.tumor, work.share.data.studies[study_name].WGS[sample_name]
      Open.link bundle.WGS.BAM.normal, work.share.data.studies[study_name].WGS[sample_name + "_normal"]
    when "CRAM"
      Open.link bundle.WGS.CRAM.tumor, work.share.data.studies[study_name].WGS[sample_name]
      Open.link bundle.WGS.CRAM.normal, work.share.data.studies[study_name].WGS[sample_name + "_normal"]
    else
      Open.link bundle.WGS[type].tumor, work.share.data.studies[study_name].WGS[sample_name]
      Open.link bundle.WGS[type].normal, work.share.data.studies[study_name].WGS[sample_name + "_normal"]
    end

    [study_name, sample_name] * ":"
  end

  dep :stage_bundle, :compute => :produce
  input :add_workflows, :array, "Workflows to add", ["HTS", "Sequence"]
  input :variant_caller, :string, "Caller to run", "mutect2"
  input :reference_type, :select, "Use original reference or bundle", "bundle", :select_options => %w(bundle hg38 b37)
  dep HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder, :referece => :placeholder do |jobname,options,dependencies|
    stage = dependencies.flatten.first

    study_name, sample_name = HTSBenchmark.helpers[:bundle_study_and_sample].call options[:bundle]

    Workflow.require_workflow "Sample"
    workflows = options[:add_workflows]
    ref_type = options[:reference_type]
    variant_caller = options[:variant_caller]
    wfs = workflows.collect{|workflow| Workflow.require_workflow workflow }

    ref_type = options[:reference_type]
    if ref_type == 'bundle'
      options[:reference] = stage.file("work").reference["hg38.fa.gz"]
    else
      options[:reference] = ref_type
      regions = stage.file('work').inputs["regions.bed"]
      options[:interval_list] = regions if regions.exists?
    end

    %w(normal tumor).collect do |type|
      name = type == 'normal' ? sample_name + '_normal' : sample_name
      f1 = stage.file('work').share.data.studies[study_name].WGS[name]["#{type}_read1.fq.gz"]
      f2 = stage.file('work').share.data.studies[study_name].WGS[name]["#{type}_read2.fq.gz"]
      {:inputs => options.merge(:fastq1 => f1, :fastq2 => f2, :sample_name => name)}
    end
  end
  dep Sample, :mutect2 do |jobname,options,dependencies|
    stage = dependencies.flatten.first

    Workflow.require_workflow "Sample"
    workflows = options[:add_workflows]
    ref_type = options[:reference_type]
    variant_caller = options[:variant_caller]
    wfs = workflows.collect{|workflow| Workflow.require_workflow workflow }

    study_name, sample_name = HTSBenchmark.helpers[:bundle_study_and_sample].call options[:bundle]

    work = stage.file('work')

    normal, tumor = dependencies.flatten.select{|d| d.task_name.to_s == "BAM"}

    ref_type = options[:reference_type]
    if ref_type == 'bundle'
      options[:reference] = stage.file("work").reference["hg38.fa.gz"]
      options[:panel_of_normals] = 'none'
    else
      options[:reference] = ref_type
      regions = stage.file('work').inputs["regions.bed"]
      options[:interval_list] = regions if regions.exists?
    end

    options["Sample#BAM"] = tumor
    options["Sample#BAM_normal"] = normal
    options[:not_overriden] = true

    options[:type_of_sequencing] = 'WES'
    options[:organism] = "Hsa"


    {:task => variant_caller, :inputs => options, :jobname => sample_name }
  end
  task :run_bundle => :text do |workflow,caller,ref_type|
    bundle, normal, tumor, job = dependencies

    job.join
    if ref_type == 'bundle' && bundle.file("inputs")["regions.bed"].exists?
      Open.cp job.path, file('sliced_variants.vcf')
      bed = bundle.file("inputs")["regions.bed"]
      HTSBenchmark.restore_sliced_vcf job.path, self.tmp_path, bed
    else
      Open.link job.path, self.tmp_path
    end
    Open.link normal.path, file('normal.bam')
    Open.link tumor.path, file('tumor.bam')
    nil
  end

  dep :run_bundle
  dep_task :eval_bundle, HTS, :vcfeval, :input_vcf => :run_bundle, :truth_vcf => :placeholder do |jobname,options,dependencies|

    bundle = HTSBenchmark.helpers[:find_bundle].call options[:bundle]
    study_name, sample_name = HTSBenchmark.helpers[:bundle_study_and_sample].call bundle


    options[:truth_vcf] = File.join(bundle, 'truth/somatic.vcf.gz')

    options[:bed_regions] = File.join(bundle, 'inputs/regions.bed')

    {:inputs => options}
  end

  #{{{ DEVELOPMENT
  #
  dep :bundle_tumor_normal
  dep_task :run_tumor_normal, HTSBenchmark, :run_bundle, :bundle => :bundle_tumor_normal

  dep :bundle_tumor_normal
  dep_task :eval_tumor_normal, HTSBenchmark, :eval_bundle, :bundle => :bundle_tumor_normal
  
  dep :bundle_population
  dep_task :run_population, HTSBenchmark, :run_bundle, :bundle => :bundle_population

  dep :bundle_population
  dep_task :eval_population, HTSBenchmark, :eval_bundle, :bundle => :bundle_population

  dep :bundle_simulate_population
  dep_task :run_simulate_population, HTSBenchmark, :run_bundle, :bundle => :bundle_population

  dep :bundle_simulate_population
  dep_task :eval_simulate_population, HTSBenchmark, :eval_bundle, :bundle => :bundle_population

  #}}} DEVELOPMENT

  #dep :stage_bundle, :compute => :produce
  #input :workflows, :array, "Workflows to add", ["HTS", "Sequence"]
  #input :variant_caller, :string, "Caller to run", "mutect2"
  #input :reference_type, :select, "Use original reference or bundle", "bundle", :select_options => %w(bundle hg38 b37)
  #input :contain_stage, :boolean, "Contain jobs in staged bundle job", true
  #dep Sample, :mutect2 do |jobname,options,dependencies|
  #  stage = dependencies.flatten.first

  #  bundle = options[:bundle]

  #  if Step === bundle
  #    bundle = bundle.files_dir 
  #  else
  #    bundle = bundle + '.files' if File.exists?(bundle + '.files')
  #  end

  #  if File.exists?(bundle)
  #    bundle =  Path.setup(bundle)
  #  else
  #    bundle = Rbbt.share.data.HTSBenchmark[bundle].find
  #  end
  #  
  #  study_name = "HTSB-St-#{Misc.digest(bundle)[0..5]}"
  #  sample_name = "HTSB-Sa-#{Misc.digest(bundle)[0..5]}"

  #  Workflow.require_workflow "Sample"
  #  workflows = options[:workflows]
  #  ref_type = options[:reference_type]
  #  variant_caller = options[:variant_caller]
  #  wfs = workflows.collect{|workflow| Workflow.require_workflow workflow }
  #  
  #  work = stage.file('work')

  #  if options[:contain_stage]
  #    wfs.each{|wf| wf.workdir = work.var.jobs[wf.to_s] }
  #    Sample.workdir = work.var.jobs.Sample
  #  end

  #  Misc.in_dir work do
  #    Sample.study_repo = work.share.data.studies
  #    Sample.study_repo.resource = nil
  #    $dont_persist_studies = true
  #    options = {}
  #    if ref_type == 'bundle'
  #      options[:reference] = stage.file("work").reference["hg38.fa.gz"]
  #    else
  #      options[:reference] = ref_type
  #      regions = stage.file('work').inputs["regions.bed"]
  #      options[:interval_list] = regions if regions.exists?
  #    end
  #    options[:panel_of_normals] = 'none'
  #    options[:type_of_sequencing] = 'WES'
  #    options[:organism] = "Hsa"
  #    {:task => variant_caller, :jobname => sample_name, :inputs => options}
  #  end
  #end
  #task :run_bundle_old => :text do |workflow,caller,ref_type|
  #  bundle, job = dependencies

  #  job.join
  #  if ref_type == 'bundle' && bundle.file("inputs")["regions.bed"].exists?
  #    Open.cp job.path, file('sliced_variants.vcf')
  #    bed = bundle.file("inputs")["regions.bed"]
  #    HTSBenchmark.restore_sliced_vcf job.path, self.tmp_path, bed
  #  else
  #    Open.link job.path, self.tmp_path
  #  end
  #  nil
  #end

end

require 'rbbt/workflow/util/orchestrator'

module HTSBenchmark

  input :workflows, :array, "Workflows to add", ["HTS", "Sequence"]
  input :task, :string, "Task to run", "mutect2"
  input :bundle, :file, "Bundle directory", nil, :nofile => true
  input :input_type, :select, "Input type FASTQ, CRAM, BAM or realign", "FASTQ", :select_options => %w(FASTQ BAM CRAM BAM-re CRAM-re)
  input :reference_type, :select, "Use original reference or bundle", "bundle", :select_options => %w(bundle hg38 b37)
  task :run_bundle => :binary do |workflows,task,bundle,type,ref_type|
    Workflow.require_workflow "Sample"
    wfs = workflows.collect{|workflow| Workflow.require_workflow workflow }
    
    work = file('work')

    wfs.each{|wf| wf.workdir = work.var.jobs[wf.to_s] }
    Sample.workdir = work.var.jobs.Sample

    bundle = bundle + '.files' if File.exists?(bundle + '.files')

    bundle =  Path.setup(bundle)


    study_name = "HTSB-St-#{Misc.digest(bundle)[0..5]}"
    sample_name = "HTSB-Sa-#{Misc.digest(bundle)[0..5]}"

    case type.to_s
    when "FASTQ"
      Open.link bundle.WGS.FASTQ.tumor, work.share.data.studies[study_name].WGS[sample_name]
      Open.link bundle.WGS.FASTQ.normal, work.share.data.studies[study_name].WGS[sample_name + "_normal"]
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

    Misc.in_dir work do
      Sample.study_repo = work.share.data.studies
      Sample.study_repo.resource = nil
      $dont_persist_studies = true
      options = {}
      if ref_type == 'bundle'
        options[:reference] = bundle.reference.glob("*.fa.gz").first
        options[:skip_rescore] = true if File.size(options[:reference]) < 1_000_000
      else
        options[:reference] = ref_type
      end
      options[:panel_of_normals] = 'none'
      job = Sample.job(task, sample_name, options)
      #Workflow::Orchestrator.process job
      job.produce
      
      if ref_type == 'bundle' && bundle.inputs["regions.bed"].exists?
        Open.cp job.path, file('sliced_variants.vcf')
        bed = bundle.inputs["regions.bed"]
        HTSBenchmark.restore_sliced_vcf job.path, self.tmp_path, bed
      else
        Open.link job.path, self.tmp_path
      end
    end
    nil
  end

  dep :run_bundle
  dep_task :eval_bundle, HTS, :vcfeval, :input_vcf => :run_bundle, :truth_vcf => :placeholder do |jobname,options,dependencies|
    options[:truth_vcf] = File.join(options[:bundle], 'truth/somatic.vcf.gz')
    {:inputs => options}
  end
end

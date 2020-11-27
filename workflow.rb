require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

Workflow.require_workflow "Sample"
Workflow.require_workflow "HTS"

module HTSBenchmark
  extend Workflow

  #CALLERS = %w(mutect2 strelka varscan_somatic muse somatic_sniper)

  #dep HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder do |jobname,options|
  #  inputs = case jobname
  #           when 'tumor'
  #             {:fastq1 => Rbbt.data.gatk["tumor.1.fastq"].find, :fastq2 => Rbbt.data.gatk["tumor.2.fastq"].find}
  #           when 'normal'
  #             {:fastq1 => Rbbt.data.gatk["normal.1.fastq"].find, :fastq2 => Rbbt.data.gatk["normal.2.fastq"].find}
  #           end
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.GATK.' + jobname}
  #end
  #extension :bam
  #task :BAM => :binary do
  #  Open.rm self.path
  #  Open.ln_s step(:BAM_rescore).path, self.path
  #  nil
  #end

  #input :bamfile, :file, "BAM file", nil, :nofile => true
  #task :sortBAM => :binary do |bamfile|
  #  args= {}

  #  args["--INPUT"] = bamfile
  #  args["--OUTPUT"] = self.tmp_path
  #  args["--SORT_ORDER"] = 'coordinate'
  #  #args.merge(:pipe => true)
  #  CMD.cmd("gatk SortSam", args)
  #  Open.mv self.tmp_path, self.path
  #end


  #input :bamfile, :file, "BAM file", nil, :nofile => true
  #input :newSM , :string, "New sample name", nil
  #task :editBAMSM => :binary do |bamfile,newSM|
  #  args = {}
  #  line = CMD.cmd("samtools view -H #{bamfile} | grep '\@RG'", args).read.split
  #  values = Hash[line.map{|i| [i.split(':')[0],i.split(':')[1] ] }]
  #  args = {}
  #  args["-I"] = bamfile
  #  args["-O"] = self.tmp_path
  #  args["-ID"] = values["ID"]
  #  args["-LB"] = values["LB"]
  #  args["-PL"] = values["PL"]
  #  args["-PU"] = values["PU"]
  #  args["-SM"] = newSM
  #  CMD.cmd("gatk AddOrReplaceReadGroups",args)
  #  Open.mv self.tmp_path, self.path
  #end

  #dep :BAM, :jobname => 'normal'
  #extension :vcf
  #dep_task :haplotype, HTS, :haplotype, :BAM => :BAM

  ##dep HTS, :revert_BAM do |jobname,options|
  ##  inputs = {:bam_file => Dir.glob(Rbbt.share.data.studies.BAMSURGEON_TEST_MUTECT["WES.orig"].bamsurgeon_test_mutect_normal["*.bam"].find)}
  ##  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
  ##end
  #dep HTS, :BAM_rescore_realign, :compute => :produce do |jobname, options|
  #  inputs = {:bam_file => Dir.glob(Rbbt.share.data.studies.BAMSURGEON_TEST_MUTECT["WES.orig"].bamsurgeon_test_mutect_normal["*.bam"].find)}
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
  #end
  #extension :bam
  #dep_task :test_muts, HTS, :BAMSurgeon_add_snvs do |jobname, options, deps|
  #  #bamfile = Samtools.prepare_BAM(deps.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign' }.first.produce.path)
  #  bamfile = deps.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign' }.first.path
  #  inputs = {:varfile => Rbbt.data.bam_surgeon["test_snvs.txt"].find,
  #            :bamfile => bamfile,
  #            :picardjar => Rbbt.software.opt.PicardTools.produce.glob("PicardTools.jar").first}
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
  #end

  #dep :test_muts
  #extension :bam
  #dep_task :sorted_test_muts, HTS, :sort_BAM, :BAM => :test_muts

  #dep :sorted_test_muts
  #dep_task :synthetic_BAM, HTSBenchmark, :editBAMSM, :bamfile => :test_muts, :newSM => "synthetic_tumor"

  #dep :sorted_test_muts
  #dep_task :synthetic_mutect, HTS, :mutect2, :tumor => :sorted_test_muts do |jobname,options,deps|
  #  normal = deps.flatten.select{|d| d.task_name == :BAM_rescore_realign}.first

  #  {:inputs => options.merge(:normal => normal), :jobname => 'HTSBenchmark.' + jobname}
  #end


  #dep :test_muts
  #dep :sortBAM do |jobname,options,deps|
  #  bamfile = deps.flatten.select{|d| d.task_name.to_s == 'test_muts' }.first.path
  #  inputs = {:bamfile => bamfile}
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.sort_tumor.' + jobname}
  #end
  #dep :editBAMSM do |jobname,options,deps|
  #  bamfile = deps.flatten.select{|d| d.task_name.to_s == 'sortBAM' }.first.path 
  #  newSM = "synthetic_tumor"
  #  inputs = {:bamfile => bamfile, :newSM => newSM}
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark' + jobname}
  #end
  #dep :sortBAM do |jobname,options,deps|
  #  bamfile = deps.first.dependencies.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign'}.first.path
  #  inputs = {:bamfile => bamfile}
  #  {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.sort_normal.' + jobname}
  #end
  #task :prepared_bams => :binary do |jobname,options,deps|
  #  Open.rm self.path 
  #  Open.ln_s step(:editBAMSM).path, self.path
  #  nil
  #end


  #dep :prepared_bams
  #dep HTS, :mutect2, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
  #  tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
  #  normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
  #  {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
  #end
  #extension :vcf
  #task :mutect_benchmark => :binary do |jobname,options,deps|
  #  Open.rm self.path
  #  Open.ln_s step(:mutect2).path, self.path
  #  nil
  #end


  #dep :prepared_bams
  #dep HTS, :varscan_somatic, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
  #  tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
  #  normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
  #  {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
  #end
  #extension :vcf
  #task :varscan_somatic_benchmark => :binary do |jobname,options,deps|
  #  Open.rm self.path
  #  Open.ln_s step(:varscan_somatic).path, self.path
  #  nil
  #end

  #dep :prepared_bams
  #dep HTS, :somatic_sniper, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
  #  tumor = deps.first.info[:dependencies].select {|a| a.include? :editBAMSM}.first[1]
  #  normal = deps.first.info[:dependencies].last[1]
  #  {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
  #end
  #extension :vcf
  #task :somatic_sniper_benchmark => :binary do |jobname,options,deps|
  #  Open.rm self.path
  #  Open.ln_s step(:somatic_sniper).path, self.path
  #  nil
  #end

  #dep :prepared_bams
  #dep HTS, :muse, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
  #  tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
  #  normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
  #  {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
  #end
  #extension :vcf
  #task :muse_benchmark => :binary do |jobname,options,deps|
  #  Open.rm self.path
  #  Open.ln_s step(:muse).path, self.path
  #  nil
  #end

  #CALLERS.each do |snv_caller|
  #  bed_task = (snv_caller + '_bed').to_sym
  #  synthetic_bam_task = ("synthetic_" + snv_caller).to_sym
  #  synthetic_bam_realign_task = ("synthetic_" + snv_caller + "_realign").to_sym

  #  dep Sample, snv_caller
  #  dep Sequence, :expanded_vcf, :vcf_file => snv_caller.to_sym
  #  task bed_task => :text do
  #    TSV.traverse step(:expanded_vcf), :into => :stream do |mutation,values,fields|
  #      values = NamedArray.setup(values, fields)
  #      mutation = values["Original"].first
  #      af = values[clean_name + ':AF']
  #      chr, pos, ref, mut_allele = mutation.split(":")
  #      next unless ref.length == 1
  #      next unless %w(C T A G).include? mut_allele
  #      [chr, pos, pos, af, mut_allele] * "\t"
  #    end
  #  end

  #  dep Sample, :BAM_normal
  #  dep bed_task
  #  extension :bam
  #  dep_task synthetic_bam_task, HTS, :BAMSurgeon_add_snvs, :bamfile => :BAM_normal, :varfile => bed_task, :picardjar => Rbbt.software.opt.PicardTools.produce.glob("PicardTools.jar").first do |jobname,options,deps|
  #    options = Sample.add_sample_options jobname, options
  #    {:inputs => options, :jobname => jobname}
  #  end

  #  dep synthetic_bam_task
  #  dep_task synthetic_bam_realign_task, HTS, :BAM_rescore_realign, :bam_file => synthetic_bam_task do |jobname,options|
  #    options = Sample.add_sample_options jobname, options
  #    {:inputs => options, :jobname => jobname}
  #  end

  #  CALLERS.each do |snv_caller2|
  #    synthetic_bam_caller_task = ("synthetic_" + snv_caller + "_" + snv_caller2).to_sym

  #    dep synthetic_bam_realign_task
  #    dep_task synthetic_bam_caller_task, HTS, snv_caller2.to_sym, :tumor => synthetic_bam_realign_task do |jobname, options, deps|
  #      normal = deps.flatten.select{|d| d.task_name == :BAM_normal }.first
  #      options = Sample.add_sample_options jobname, options
  #      {:inputs => options.merge(:normal => normal), :jobname => jobname}
  #    end
  #  end
  #end



end
require 'tasks/simulate/genotypes'
require 'tasks/simulate/NEAT'
require 'tasks/miniref.rb'
require 'tasks/benchmark.rb'
require 'tasks/bundle.rb'
require 'tasks/evaluate/vcfeval'
require 'tasks/evaluate/IGV'

#require 'rbbt/knowledge_base/HTSBenchmark'
#require 'rbbt/entity/HTSBenchmark'
#require 'germline'

require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/HTSBenchmark'

Workflow.require_workflow "HTS"
module HTSBenchmark
  extend Workflow

   CALLERS = %w(mutect2 strelka varscan_somatic muse somatic_sniper)

   dep HTS, :BAM, :fastq1 => :placeholder, :fastq2 => :placeholder do |jobname,options|
     inputs = case jobname
              when 'tumor'
                {:fastq1 => Rbbt.data.gatk["tumor.1.fastq"].find, :fastq2 => Rbbt.data.gatk["tumor.2.fastq"].find}
              when 'normal'
                {:fastq1 => Rbbt.data.gatk["normal.1.fastq"].find, :fastq2 => Rbbt.data.gatk["normal.2.fastq"].find}
              end
     {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.GATK.' + jobname}
   end
   extension :bam
   task :BAM => :binary do
     Open.rm self.path
     Open.ln_s step(:BAM_rescore).path, self.path
     nil
   end

   input :bamfile, :file, "BAM file", nil, :nofile => true
   task :sortBAM => :binary do |bamfile|
      args= {}
	
      args["--INPUT"] = bamfile
      args["--OUTPUT"] = self.tmp_path
      args["--SORT_ORDER"] = 'coordinate'
      #args.merge(:pipe => true)
      CMD.cmd("gatk SortSam", args)
      Open.mv self.tmp_path, self.path
   end


   input :bamfile, :file, "BAM file", nil, :nofile => true
   input :newSM , :string, "New sample name", nil
   task :editBAMSM => :binary do |bamfile,newSM|
      args = {}
      line = CMD.cmd("samtools view -H #{bamfile} | grep '\@RG'", args).read.split
      values = Hash[line.map{|i| [i.split(':')[0],i.split(':')[1] ] }]
      args = {}
	  args["-I"] = bamfile
      args["-O"] = self.tmp_path
	  args["-ID"] = values["ID"]
	  args["-LB"] = values["LB"]
	  args["-PL"] = values["PL"]
	  args["-PU"] = values["PU"]
	  args["-SM"] = newSM
	  CMD.cmd("gatk AddOrReplaceReadGroups",args)
 	  Open.mv self.tmp_path, self.path
   end

   dep :BAM, :jobname => 'normal'
   extension :vcf
   dep_task :haplotype, HTS, :haplotype, :BAM => :BAM

   dep HTS, :revert_BAM do |jobname,options|
      inputs = {:bam_file => Dir.glob(Rbbt.share.data.studies.BAMSURGEON_TEST_MUTECT["WES.orig"].bamsurgeon_test_mutect_normal["*.bam"].find)}
      {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
   end
   dep HTS,:BAM_rescore_realign
   dep HTS, :BAMSurgeon_add_snvs do |jobname, options, deps|
      bamfile = Samtools.prepare_BAM(deps.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign' }.first.path)
      inputs = {:varfile => Rbbt.data.bam_surgeon["test_snvs.txt"].find,
                :bamfile => bamfile,
                :picardjar => Rbbt.software.opt.PicardTools["picard-2.20.5.jar"].find}

      {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
   end
   extension :bam
   task :test_muts => :binary do |jobname,options,deps|
      Open.rm self.path
      Open.ln_s step(:BAMSurgeon_add_snvs).path, self.path
      nil
   end

   dep :test_muts
   dep :sortBAM do |jobname,options,deps|
      bamfile = deps.flatten.select{|d| d.task_name.to_s == 'test_muts' }.first.path
      inputs = {:bamfile => bamfile}
      {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.sort_tumor.' + jobname}
   end
   dep :editBAMSM do |jobname,options,deps|
      bamfile = deps.flatten.select{|d| d.task_name.to_s == 'sortBAM' }.first.path 
      newSM = "synthetic_tumor"
      inputs = {:bamfile => bamfile, :newSM => newSM}
      {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark' + jobname}
   end
   dep :sortBAM do |jobname,options,deps|
      bamfile = deps.first.dependencies.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign'}.first.path
      iii bamfile
      inputs = {:bamfile => bamfile}
      {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.sort_normal.' + jobname}
   end
   task :prepared_bams => :binary do |jobname,options,deps|
      Open.rm self.path 
      Open.ln_s step(:editBAMSM).path, self.path
	  nil
   end


   dep :prepared_bams
   dep HTS, :mutect2, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
      tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
      normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
      {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
   end
   extension :vcf
   task :mutect_benchmark => :binary do |jobname,options,deps|
      Open.rm self.path
      Open.ln_s step(:mutect2).path, self.path
      nil
   end

   
   dep :prepared_bams
   dep HTS, :varscan_somatic, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
      tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
      normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
      {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
   end
   extension :vcf
   task :varscan_somatic_benchmark => :binary do |jobname,options,deps|
      Open.rm self.path
      Open.ln_s step(:varscan_somatic).path, self.path
      nil
   end

   dep :prepared_bams
   dep HTS, :somatic_sniper, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
      tumor = deps.first.info[:dependencies].select {|a| a.include? :editBAMSM}.first[1]
      normal = deps.first.info[:dependencies].last[1]
      {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
   end
   extension :vcf
   task :somatic_sniper_benchmark => :binary do |jobname,options,deps|
      Open.rm self.path
      Open.ln_s step(:somatic_sniper).path, self.path
      nil
   end

   dep :prepared_bams
   dep HTS, :muse, :tumor => :placeholder, :normal => :placeholder do |jobname,options,deps|
      tumor = deps.flatten.last.info[:dependencies].select{|d| d.include?:editBAMSM}.first.last
      normal =deps.flatten.last.info[:dependencies].select{|d| d.include?:sortBAM}.select{|d| d.last.include?("sort_normal")}.first.last
      iii normal
      {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
   end
   extension :vcf
   task :muse_benchmark => :binary do |jobname,options,deps|
      Open.rm self.path
      Open.ln_s step(:muse).path, self.path
      nil
   end
end
#require 'HTSBenchmark/tasks/basic.rb'

#require 'rbbt/knowledge_base/HTSBenchmark'
#require 'rbbt/entity/HTSBenchmark'


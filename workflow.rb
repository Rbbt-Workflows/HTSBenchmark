require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

Workflow.require_workflow "Sample"
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

  dep HTS, :BAM_rescore_realign, :compute => :produce do |jobname, options|
    inputs = {:bam_file => Dir.glob(Rbbt.share.data.studies.BAMSURGEON_TEST_MUTECT["WES.orig"].bamsurgeon_test_mutect_normal["*.bam"].find)}
    {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
  end
  extension :bam
  dep_task :test_muts, HTS, :BAMSurgeon_add_snvs do |jobname, options, deps|
    bamfile = deps.flatten.select{|d| d.task_name.to_s == 'BAM_rescore_realign' }.first.path
    inputs = {:varfile => Rbbt.data.bam_surgeon["test_snvs.txt"].find,
              :bamfile => bamfile,
              :picardjar => Rbbt.software.opt.PicardTools.produce.glob("PicardTools.jar").first}

    {:inputs => options.merge(inputs), :jobname => 'HTSBenchmark.' + jobname}
  end

  dep :test_muts
  extension :bam
  dep_task :sorted_test_muts, HTS, :sort_BAM, :BAM => :test_muts

  dep :sorted_test_muts
  dep_task :synthetic_BAM, HTSBenchmark, :editBAMSM, :bamfile => :test_muts, :newSM => "synthetic_tumor"

  dep :sorted_test_muts
  dep_task :synthetic_mutect, HTS, :mutect, :tumor => :sorted_test_muts do |jobname,options,deps|
    normal = deps.flatten.select{|d| d.task_name == :BAM_rescore_realign}.first

    {:inputs => options.merge(:normal => normal), :jobname => 'HTSBenchmark.' + jobname}
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

  A = 0; C = 1; G = 2; T = 3
  CALLERS2 = %w(mutect2 strelka varscan_somatic muse somatic_sniper)
  CALLERS2.each do |snv_caller|
    bed_task = (snv_caller + '_bed').to_sym
    synthetic_bam_task = ("synthetic_" + snv_caller).to_sym
    synthetic_bam_realign_task = ("synthetic_" + snv_caller + "_realign").to_sym
    synthetic_bam_editBAMSM_task = ("synthetic_" + snv_caller + "_editBAMSM").to_sym
    synthetic_bam_sortBAM_task = ("synthetic_" + snv_caller + "_sortBAM").to_sym

    dep Sample, snv_caller, :somatic_score => 40
    dep Sequence, :expanded_vcf, :vcf_file => snv_caller.to_sym
    task bed_task => :text do
      TSV.traverse step(:expanded_vcf), :into => :stream do |mutation,values,fields|
        values = NamedArray.setup(values, fields)
        mutation = values["Original"].first
        af = nil
        case snv_caller
        when "muse"
          alt = values["TUMOR:AD"].first.split(",").last.to_f
          dp = values["TUMOR:DP"].first.to_f
          af = (alt/dp).to_s
        when "somatic_sniper"
          alt = values["Original"].first[-1]
          bcount_field = fields.select{|f| f.include? 'BCOUNT'}.last 
          bcount = values[bcount_field].first
          total_bcount = bcount.split(',').map(&:to_i).inject(0,:+)
          alt_bcount = bcount.split(',')[Object.module_eval(alt)]
          af = (alt_bcount.to_f/total_bcount.to_f).to_s
        when "mutect2"
          af_field = fields.select{|f| f.include? 'AF'}.first
          af = values[af_field]
        when "varscan_somatic"
          af = values["TUMOR:FREQ"].first[0..-2].to_f / 100.0
        when "strelka"
          next unless values["Filter"].first == "PASS" 
          alt = values["Original"].first[-1]
		  total_dp = values["TUMOR:DP"].first.to_f - values["TUMOR:FDP"].first.to_f
		  alt_count = values["TUMOR:" + alt + "U"].first.split(',').first.to_f                 
		  af = (alt_count/total_dp).to_s 
        end
        chr, pos, ref, mut_allele = mutation.split(":")
        next unless ref.length == 1
        next unless %w(C T A G).include? mut_allele
        [chr, pos, pos, af, mut_allele] * "\t"
      end
    end

    dep Sample, :BAM_normal
    dep bed_task
    extension :bam
    dep_task synthetic_bam_task, HTS, :BAMSurgeon_add_snvs, :bamfile => :BAM_normal, :varfile => bed_task, :picardjar => Rbbt.software.opt.PicardTools.produce.glob("PicardTools.jar").first do |jobname,options,deps|
      options = Sample.add_sample_options jobname, options
      {:inputs => options, :jobname => jobname}
    end
    
    dep synthetic_bam_task 
    dep_task synthetic_bam_sortBAM_task, HTS, :sort_BAM, :BAM => synthetic_bam_task
    
    dep synthetic_bam_sortBAM_task
    dep_task synthetic_bam_editBAMSM_task, HTSBenchmark, :editBAMSM, :bamfile => synthetic_bam_sortBAM_task, :newSM => "synthetic_tumor_" + snv_caller.to_s
  

    dep synthetic_bam_editBAMSM_task
    dep_task synthetic_bam_realign_task, HTS, :BAM_rescore_realign, :bam_file => synthetic_bam_editBAMSM_task
     
    CALLERS.each do |snv_caller2|
      synthetic_bam_caller_task = ("synthetic_" + snv_caller + "_" + snv_caller2).to_sym
      
      dep Sample, :BAM_normal 
      dep synthetic_bam_realign_task
      dep_task synthetic_bam_caller_task, HTS, snv_caller2.to_sym, :tumor => synthetic_bam_realign_task do |jobname, options, deps|
        normal = deps.flatten.select{|d| d.task_name == :BAM_normal }.first
        options = Sample.add_sample_options jobname, options
        {:inputs => options.merge(:normal => normal), :jobname => jobname}
      end
    end
  end




end

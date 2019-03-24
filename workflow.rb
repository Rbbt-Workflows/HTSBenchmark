require 'rbbt-util'
require 'rbbt/workflow'

Misc.add_libdir if __FILE__ == $0

#require 'rbbt/sources/HTSBenchmark'

Workflow.require_workflow "HTS"
module HTSBenchmark
  extend Workflow

  def self.af_not_in_resource
    0.0000025
  end

  def self.germline_resource
    GATK.known_sites.b37["af-only-gnomad.vcf.gz"].produce.find
	end

  CALLERS = %w(mutect2_clean strelka varscan_somatic varscan_somatic_alt)

  dep HTS, :BAM_rescore, :fastq1 => :placeholder, :fastq2 => :placeholder do |jobname,options|
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

  CALLERS.each do |var_caller|
    dep :BAM, :compute => :bootstrap do |jobname, options|
      %w(tumor normal).collect do |type|
        {:inputs => options, :jobname => type}
      end
    end
    extension :vcf
    dep_task var_caller, HTS, var_caller, :tumor => :BAM, :normal => :BAM,
      :germline_resource => HTSBenchmark.germline_resource,
      :af_not_in_resource => HTSBenchmark.af_not_in_resource do |jobname,options,deps|
        normal = deps.flatten.select{|d| d.task_name.to_s == 'BAM' && d.clean_name.include?('normal')}.first
        tumor = deps.flatten.select{|d| d.task_name.to_s == 'BAM' && d.clean_name.include?('tumor')}.first

        {:inputs => options.merge({:normal => normal, :tumor => tumor}), :jobname => jobname}
    end
  end

  dep :strelka, :compute => :bootstrap, :tumor => :placeholder, :normal => :placeholder do |jobname, options|
    CALLERS.collect do |var_caller|
      {:task => var_caller, :inputs => options, :jobname => jobname}
    end
  end
  task :caller_comparison => :tsv do
    caller_variants = {}
    dependencies.each do |dep|
      var_caller = dep.task_name

      variants = TSV.traverse dep, :into => [] do |chr, values|
        pos, id, ref, alt, *rest = values

        [chr, pos, alt] * ":"
      end
      caller_variants[var_caller] = variants
    end

    var_callers = caller_variants.keys.sort
    tsv = TSV.setup({}, :key_field => "Caller", :fields => ["Unique"] + var_callers, :type => :list, :cast => :to_i)
    var_callers.each do |c1|
      m1 = caller_variants[c1]
      matches = var_callers.collect do |c2|
        m2 = caller_variants[c2]
        m1 & m2
      end
      unique = m1 - matches.flatten
      tsv[c1] = [unique.length] + matches.collect{|m| m.length}
    end

    tsv
  end
end

#require 'HTSBenchmark/tasks/basic.rb'

#require 'rbbt/knowledge_base/HTSBenchmark'
#require 'rbbt/entity/HTSBenchmark'


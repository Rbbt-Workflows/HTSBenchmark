require 'rbbt-util'
require 'rbbt/workflow'
require 'rbbt/workflow/remote_workflow'

module HTSBenchmark
  GERMLINE_RESULT = 0.992642
  Rbbt.claim Rbbt.share.data.studies.GIAB, :proc do |directory|
	  url_mother    = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam"
	  url_father    = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam"
    url_son       = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam" 
	  url_truthset  = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015/HG002-HG003-HG004.jointVC.filter.vcf"
	  url_intervals = "ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/OsloUniversityHospital_Exome_GATK_jointVC_11242015/wex_Agilent_SureSelect_v05_b37.baits.slop50.merged.list"
    Open.mkdir directory.data
    Open.mkdir directory["WES.orig"]

	  `wget #{url_mother} -O #{directory}/data/mother.bam`
    `wget #{url_father} -O #{directory}/data/father.bam`
    `wget #{url_son} -O #{directory}/data/son.bam`
    `wget #{url_truthset} -O #{directory}/data/truthset.vcf`
    `wget #{url_intervals} -O #{directory}/data/intervals.list`
    
    Path.setup(directory)
    Open.mkdir directory["WES.orig"].Mother
    Open.mkdir directory["WES.orig"].Father
    Open.mkdir directory["WES.orig"].Son

    Open.ln_s directory.data["mother.bam"], directory["WES.orig"].Mother["mother.bam"] unless File.exists?(directory["WES.orig"].Mother["mother.bam"].find)
    Open.ln_s directory.data["father.bam"], directory["WES.orig"].Father["father.bam"] unless File.exists?(directory["WES.orig"].Father["father.bam"].find)
    Open.ln_s directory.data["son.bam"], directory["WES.orig"].Son["son.bam"] unless File.exists?(directory["WES.orig"].Father["father.bam"].find)
    nil
  end

  dep Sample, :BAM do |jobname,options|
    Rbbt.share.data.studies.GIAB.produce.find
    ["Mother", "Father", "Son"].collect do |s|
      {:inputs => options, :jobname => s}
    end
  end
  dep HTS, :haplotype do |jobname,options,dependencies|
    options[:BAM_list] = dependencies.flatten
    options[:interval_list] = Rbbt.share.data.studies.GIAB.data["intervals.list"].find
    {:inputs => options}
  end
  task :hap_py do
    CMD.cmd_log("export HGREF=#{Organism["Hsa"].b37["b37.fa"].find}; hap.py #{step(:haplotype).path} #{Rbbt.share.data.studies.GIAB.data["truthset.vcf"].find} -o test_GIAB 2>/dev/null > #{self.tmp_path}")
    Open.ln_s self.tmp_path, self.path
  end

  dep :hap_py
  task :test_GIAB do
    tsv = TSV.open(step(:hap_py), :sep => '\t', :header_hash => '')
    result = tsv.flatten[6].split(' ')[11].to_f
    raise unless result == GERMLINE_RESULT
  end
end

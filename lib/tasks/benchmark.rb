module HTSBenchmark
  
  dep :bundle
  input :type, :select, "Type of benchmark", 'miniref', :select_options => %w(miniref fullref golden)
  task :benchmark => :text do  |type|
    bundle = step('bundle').file('bundle')
    stage = file('stage')

    type = type.to_s
    sample = case type
             when "miniref"
               "ARGO-Miniref"
             when "fullref"
               "ARGO-fullref"
             when "golden"
               "ARGO-golden"
             end

    set_info :sample, sample

    if type == "miniref" || type == 'golden'
      Misc.in_dir stage.share.organisms.Hsa.hg38 do
        CMD.cmd("ln -s '#{bundle.reference}'/* . ")
        CMD.cmd("ln -s '#{bundle.known_sites}' . ")
      end
    end

    Open.write(stage.share.data.studies["ARGO-Benchmark"].options.reference, "hg38")

    if type == "golden"
      Misc.in_dir stage.share.data.studies["ARGO-Benchmark"].WES[sample] do
        CMD.cmd("ln -s '#{bundle.golden_BAM["tumor_golden.bam"]}' .")
      end
      Misc.in_dir stage.share.data.studies["ARGO-Benchmark"].WES["#{sample}_normal"] do
        CMD.cmd("ln -s '#{bundle.golden_BAM["normal_golden.bam"]}' .")
      end
    else
      Misc.in_dir stage.share.data.studies["ARGO-Benchmark"].WES do
        CMD.cmd("ln -s '#{bundle.FASTQ.tumor}' #{sample}")
        CMD.cmd("ln -s '#{bundle.FASTQ.normal}' #{sample}_normal")
      end
    end

    "DONE"
  end

  dep :benchmark
  input :variant_caller, :select, "Variant caller to use", "mutect2", :select_options => %w(mutect2 varscan somatic_sniper muse pindel delly svaba strelka combined_caller_vcfs consensus_somatic_variants)
  input :min_callers, :integer, "Mininum number of callers for consensus_somatic_variants", 2
  extension :vcf
  task :somatic_variant_calling => :text do |variant_caller,min_callers|
    benchmark = step(:benchmark)
    stage = benchmark.file('stage')

    sample = benchmark.info[:sample]

    Misc.in_dir benchmark.file('stage') do
      dir = benchmark.file('stage').share.data.studies["ARGO-Benchmark"].WES

      bam = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 BAM -pf -prov --workdir_all #{stage.workdir}").read.strip
      bam_normal = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 BAM_normal -pf -prov --workdir_all #{stage.workdir}").read.strip

      if Persist.newer?(bam, dir, true)
        Log.warn "Recursive clean of variant calling"
        clean_str = "-rcl"
      else
        clean_str = "-cl"
      end

      mnp = 1

      CMD.cmd_log("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 #{variant_caller} --max_mnp_distance=#{mnp} --min_callers #{min_callers} -ck HTS_high #{clean_str} --update -pf --workdir_all #{stage.workdir}")

      vcf = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 #{variant_caller} --max_mnp_distance=#{mnp} --min_callers #{min_callers} -pf --workdir_all #{stage.workdir}", :log => true).read.strip

      Open.rm file('BAM.bam')
      Open.rm file('BAM_normal.bam')

      Open.link bam, file('BAM.bam')
      Open.link bam_normal, file('BAM_normal.bam')

      Open.cp vcf, self.tmp_path
    end
    nil
  end

  dep :benchmark
  extension :vcf
  task :germline_variant_calling => :text do
    benchmark = step(:benchmark)
    stage = benchmark.file('stage')

    sample = benchmark.info["sample"]

    Misc.in_dir benchmark.file('stage') do
      dir = benchmark.file('stage').share.data.studies["ARGO-Benchmark"].WES

      bam_normal = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 BAM_normal -pf -prov --workdir_all #{stage.workdir}").read.strip

      if Persist.newer?(bam, dir, true)
        Log.warn "Recursive clean of variant calling"
        clean_str = "-rcl"
      else
        clean_str = "-cl"
      end

      vcf = CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn #{sample} --log 0 haplotype -ck HTS_high #{clean_str} --update -pf --workdir_all #{stage.workdir}", :log => true).read.strip

      Open.rm file('BAM_normal.bam')

      Open.link bam_normal, file('BAM_normal.bam')

      Open.cp vcf, self.tmp_path
    end
  end

end

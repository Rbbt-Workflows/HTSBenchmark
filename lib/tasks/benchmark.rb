module HTSBenchmark
  dep :simulate_germline_hg38_vcf
  dep :simulate_somatic_hg38_vcf
  dep :miniref, :vcf => :simulate_somatic_hg38_vcf, :jobname => 'hg38'
  dep :minify_vcf, :vcf_file => :simulate_germline_hg38_vcf, :sizes => :miniref_sizes, :jobname => "germline"
  dep :minify_vcf, :vcf_file => :simulate_somatic_hg38_vcf, :sizes => :miniref_sizes, :jobname => "somatic"
  dep :NEAT_genreads, :somatic => :placeholder, :germline => :placeholder, :reference => :miniref do |jobname,options,dependencies|
    germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}
    options[:germline] = germline
    options[:somatic] = somatic
    {:inputs => options, :jobname => jobname}
  end
  task :benchmark => :text do
    stage = file('stage')

    Misc.in_dir stage.share.organisms.Hsa do
      CMD.cmd("ln -s '#{step(:miniref).file('hg38')}' . ")
    end

    Misc.in_dir stage.share.data.studies.Miniref do
      CMD.cmd("ln -s '#{step(:NEAT_genreads).file('output')}' data ")
    end

    Open.write(stage.share.data.studies.Miniref.options.reference, "hg38")

    Misc.in_dir stage.share.data.studies.Miniref.WES.Miniref do
      CMD.cmd("ln -s ../../data/tumo* .")
    end

    Misc.in_dir stage.share.data.studies.Miniref.WES.Miniref_normal do
      CMD.cmd("ln -s ../../data/normal* .")
    end

    "DONE"
  end

end

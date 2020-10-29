module HTSBenchmark

  dep :benchmark
  input :variant_caller, :select, "Variant caller to use", "mutect2", :select_options => %w(mutect2 varscan somatic_sniper muse strelka)
  input :min_callers, :integer, "Mininum number of callers", 2
  extension :vcf
  task :somatic_variant_calling => :text do |variant_caller,min_callers|
    benchmark = step(:benchmark)

    Misc.in_dir benchmark.file('stage') do
      CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 #{variant_caller} --min_callers #{min_callers} -ck HTS_high -cl", :log => true, :pipe => true)
    end
  end

  dep :benchmark
  extension :vcf
  task :germline_variant_calling => :text do
    benchmark = step(:benchmark)

    Misc.in_dir benchmark.file('stage') do
      CMD.cmd("env PWD=#{benchmark.file('stage')} rbbt task Sample -W HTS -jn Miniref --log 0 haplotype -ck HTS_high", :log => true, :pipe => true)
    end
  end

  input :input_vcf, :file
  input :truth_vcf, :file
  dep Sequence, :genomic_mutations, :vcf_file => :input_vcf
  dep Sequence, :genomic_mutations, :vcf_file => :truth_vcf
  task :compare_vcf => :yaml do |input,truth|
    inputm, truthm = dependencies.collect{|dep| dep.load.collect{|m| m.sub('chr', '')} }  
    common = inputm & truthm
    input = inputm - truthm
    truth = truthm - inputm
    Open.write(file('common'), common * "\n")
    Open.write(file('input'), input * "\n")
    Open.write(file('truth'), truth * "\n")
    {
      :common => common.length,
      :input => input.length,
      :truth => truth.length,
    }

  end


  dep :germline_variant_calling, :compute => :produce
  dep_task :compare_germline_vcf, HTSBenchmark, :compare_vcf, :input_vcf => :germline_variant_calling, :truth_vcf => :placeholder do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    benchmark = variants.step(:benchmark)
    neat = variants.step(:NEAT_genreads)
    germline, somatic = benchmark.dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

    options[:truth_vcf] = germline
    options[:input_vcf] = variants
    {:inputs => options, :jobname => "Miniref"}
  end

  dep :somatic_variant_calling, :compute => :produce
  dep_task :compare_somatic_vcf, HTSBenchmark, :compare_vcf, :input_vcf => :somatic_variant_calling, :truth_vcf => :placeholder do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    benchmark = variants.step(:benchmark)
    neat = benchmark.step(:NEAT_genreads)
    germline, somatic = benchmark.dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

    options[:truth_vcf] = somatic
    options[:input_vcf] = variants
    {:inputs => options, :jobname => "Miniref"}
  end


  dep :germline_variant_calling, :compute => :produce
  dep_task :hap_py, HTS, :hap_py, :input_vcf => :germline_variant_calling, :truth_vcf => :placeholder, :reference => :placeholder do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    neat = variants.step(:NEAT_genreads)
    benchmark = variants.step(:benchmark)
    germline, somatic = benchmark.dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

    options[:truth_vcf] = germline
    options[:input_vcf] = variants

    options[:reference] = benchmark.file('stage/share/organisms/Hsa/hg38/hg38.fa.gz')
    {:inputs => options, :jobname => "Miniref"}
  end

  dep :somatic_variant_calling, :compute => :produce
  dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :somatic_variant_calling, :truth_vcf => :placeholder, :reference => :placeholder do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    neat = variants.step(:NEAT_genreads)
    benchmark = variants.step(:benchmark)
    germline, somatic = benchmark.dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

    options[:truth_vcf] = somatic
    options[:input_vcf] = variants

    options[:reference] = benchmark.file('stage/share/organisms/Hsa/hg38/hg38.fa.gz')
    {:inputs => options, :jobname => "Miniref"}
  end
end

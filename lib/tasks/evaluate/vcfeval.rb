module HTSBenchmark
  dep :somatic_variant_calling
  dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :somatic_variant_calling, :truth_vcf => :placeholder, :reference => 'hg38' do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    bundle = variants.step(:bundle)
    benchmark = variants.step(:benchmark)

    germline = bundle.file('bundle/truth/germline.vcf')
    somatic = bundle.file('bundle/truth/somatic.vcf')
    #germline, somatic = bundle.rec_dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

    options[:truth_vcf] = somatic
    options[:input_vcf] = variants

    miniref = benchmark.file('stage/share/organisms/Hsa/hg38/hg38.fa.gz')
    options[:reference] = miniref if miniref.exists?
    {:inputs => options, :jobname => "Miniref"}
  end

  dep :vcfeval
  task :fp_positions => :array do
    fp = step(:vcfeval).file('output/fp.vcf.gz')
    Sequence.job(:genomic_mutations, nil, :vcf_file => fp).run
  end

  dep :vcfeval
  task :fn_positions => :array do
    fn = step(:vcfeval).file('output/fn.vcf.gz')
    Sequence.job(:genomic_mutations, nil, :vcf_file => fn).run
  end
end

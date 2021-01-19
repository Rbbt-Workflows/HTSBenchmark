module HTSBenchmark

  dep :bundle
  dep :somatic_variant_calling
  task :fixed_truth_vcf => :text do
    bundle = step(:bundle)
    input = step(:somatic_variant_calling).path
    truth = bundle.file('bundle/truth/somatic.vcf')

    Misc.in_dir files_dir do
      truth_sorted = File.join('.', "truth.vcf")

      truth_sorted_orig = File.join('.', "truth.orig.vcf")
      truth_sorted_tmp = File.join('.', "truth.tmp.vcf")
      truth_sorted_tmp2 = File.join('.', "truth.tmp2.vcf")

      truth_io = TSV.get_stream truth
      Open.write(truth_sorted_orig, truth_io)

      CMD.cmd("echo '##fileformat=VCFv4.2' > #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##' #{truth_sorted_orig} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##contig' #{input} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##contig' #{input} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep '##FORMAT' #{input} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("echo '##FORMAT=<ID=,Number=R,Type=Integer,Description=>' >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("echo '##FORMAT=<ID=GT,Number=R,Type=String,Description=>' >> #{truth_sorted_tmp}", :nofail => true)

      ##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the
      CMD.cmd("grep '#CHR' #{truth_sorted_orig} >> #{truth_sorted_tmp}", :nofail => true)
      CMD.cmd("grep -v '#' #{truth_sorted_orig} |grep -v _alt| grep -v _random >> #{truth_sorted_tmp}", :nofail => true)

      Open.open(truth_sorted_tmp) do |io|
        Open.open(truth_sorted_tmp2, :mode => 'w') do |file|
          TSV.traverse io, :type => :array do |line|
            value = case line
                    when /^##/
                      line
                    when /^#/
                      line + "\t" + "FORMAT" + "\t"  + clean_name
                    else
                      next if line.split("\t")[4].split(",").include? line.split("\t")[3]
                      line = (line.split("\t")[0..4] + ["", ".",""]) * "\t"
                      if line =~ /^chr/ || ! reference.include?("hg38")
                        line + "\t" + "GT" + "\t"  + "0/1"
                      else
                        "chr" + line + "\t" + "GT" + "\t"  + "0/1"
                      end
                    end
            file.puts value
          end
        end
      end
    end
    Open.mv file('truth.tmp2.vcf'), self.tmp_path
    nil
  end

  #dep :somatic_variant_calling
  #dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :somatic_variant_calling, :truth_vcf => :placeholder, :reference => 'hg38' do |jobname, options,dependencies|
  #  variants = dependencies.flatten.first
  #  bundle = variants.step(:bundle)
  #  benchmark = variants.step(:benchmark)

  #  germline = bundle.file('bundle/truth/germline.vcf')
  #  somatic = bundle.file('bundle/truth/somatic.vcf')
  #  #germline, somatic = bundle.rec_dependencies.select{|dep| dep.task_name.to_s == "minify_vcf"}

  #  options[:truth_vcf] = somatic
  #  options[:input_vcf] = variants

  #  miniref = benchmark.file('stage/share/organisms/Hsa/hg38/hg38.fa.gz')
  #  options[:reference] = miniref if miniref.exists?
  #  {:inputs => options, :jobname => "Miniref"}
  #end

  dep :somatic_variant_calling
  dep :fixed_truth_vcf
  dep_task :vcfeval, HTS, :vcfeval, :input_vcf => :somatic_variant_calling, :truth_vcf => :fixed_truth_vcf, :reference => 'hg38' do |jobname, options,dependencies|
    variants = dependencies.flatten.first
    benchmark = variants.step(:benchmark)

    miniref = benchmark.file('stage/share/organisms/Hsa/hg38/hg38.fa.gz')
    options[:reference] = miniref if miniref.exists?
    {:inputs => options, :jobname => "Miniref"}
  end

  dep :vcfeval
  dep_task :fp_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fp.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end

  dep :vcfeval
  dep_task :fn_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fn.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end
end

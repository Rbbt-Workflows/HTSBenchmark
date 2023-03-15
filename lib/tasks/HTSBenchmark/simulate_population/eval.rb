module HTSBenchmark
  dep :choose_signature_mutations, :jobname => "Default"
  input :clone, :string, "Clone name"
  task :chosen_mutations => :array do |clone|
    tsv = dependencies.first.load
    if clone.nil?
      tsv.values.flatten
    else
      raise RbbtException, "Clone not found #{clone}" unless tsv.include? clone
      tsv[clone]
    end
  end

  dep :genotype_somatic_hg38, :jobname => 'Default', :mutations_per_MB => 0
  dep :merge_signatures
  dep :chosen_mutations, :clone_signatures => :merge_signatures, :mutations => :genotype_somatic_hg38
  dep_task :predict_clonal_signatures, MutationSignatures, :assign_signatures,
    :mutations => :chosen_mutations, :organism => HTSBenchmark.organism, :snvs_only => true, :watson => true


  dep :predict_clonal_signatures, :clone => :placeholder do |jobname,options|
    clones = options[:signature_proportions].fields
    clones.collect do |clone|
      {:inputs => options.merge(:clone => clone), :jobname => clone}
    end
  end
  task :eval_subclonal_signatures => :tsv do
    Step.wait_for_jobs dependencies
    predicted_signatures = dependencies.
      select{|d| d.task_name == :predict_clonal_signatures }.
      inject(nil) do |acc,d|
        tsv = d.load
        tsv[d.clean_name] = tsv.delete("Sample")
        acc = acc.nil? ? tsv : acc.merge!(tsv)
      end

    signature_proportions = recursive_inputs[:signature_proportions].transpose.to_list

    signatures = signature_proportions.fields
    predicted_signatures = predicted_signatures.slice signatures

    res = TSV.setup({}, :key_field => "Clone", :fields => %w(cosine similarity), :type => :single, :cast => :to_f)

    clones = predicted_signatures.keys
    clones.each do |clone|
      cosine = signature_proportions[clone].zip(predicted_signatures[clone]).collect do |prop,pred|
        prop.to_f * pred.to_f
      end
      res[clone] = Misc.sum(cosine)
    end

    res
  end

  dep :genotype_somatic_hg38, :jobname => 'Default', :mutations_per_MB => 0
  dep :chosen_mutations, :mutations => :genotype_somatic_hg38
  task :signature_chr_sizes => :tsv do
    chr_sizes = TSV.setup({}, "Chromosome~Size#:type=:single#:cast=:to_f")
    mutations = dependencies.last.load
    mutations.each do |mutation|
      chr, pos, alt = mutation.split(":")
      chr_sizes[chr] = [chr_sizes[chr], pos.to_i].compact.max
    end
    total = Misc.sum(chr_sizes.values)
    chr_sizes["Total MB"] = (mbs = total.to_f / 1_000_000).round(2)
    chr_sizes["Total mutations"] = mutations.length
    chr_sizes["mutations per MB"] = (mutations.length.to_f / mbs).round(2)
    chr_sizes
  end

  dep :signature_chr_sizes
  dep :eval_subclonal_signatures
  task :eval_signatures => :tsv do
    Step.wait_for_jobs dependencies
    tsv = TSV.setup({}, :key_field => "Statistic", :fields => ["Value"], :type => :single)
    chr_sizes = step(:signature_chr_sizes).load
    subclonal = step(:eval_subclonal_signatures).load

    tsv["Total MB"] = chr_sizes["Total MB"]
    tsv["Total mutations"] = chr_sizes["Total mutations"]
    tsv["mutations per MB"] = chr_sizes["mutations per MB"]
    tsv.merge!(subclonal)
  end

end

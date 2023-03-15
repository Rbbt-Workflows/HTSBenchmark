module HTSBenchmark
  input :signature_proportions, :tsv, "Signature activity proportion", nil, :required => true
  task :merge_signatures => :tsv do |signature_proportions|

    signatures = Rbbt.share.databases.MutationSignatures.COSMIC_SBS.tsv
    signature_codes = signature_proportions.keys & signatures.fields

    clones = signature_proportions.fields

    acc = nil
    clones.each do |clone|
      log :clone, "Processing #{clone}"
      clone_proportions = signature_proportions.column(clone).to_single
      clone_proportions.unnamed = true
      proportions = clone_proportions.values_at(*signature_codes)

      res = signatures.slice(signature_codes).to_single{|vs| Misc.sum(vs.zip(proportions).collect{|p| p[0] * p[1] }) }
      res.fields = [clone]
      res.unnamed = true

      clone_proportions.each do |triplet,value|
        next if signature_codes.include? triplet
        res[triplet] ||= 0
        res[triplet] += value.to_f
      end

      acc = acc ? acc.attach(res) : res
    end

    acc
  end

  dep :merge_signatures
  dep MutationSignatures, :context, :organism => HTSBenchmark.organism
  input :place_all, :boolean, "Make sure that all mutations are placed to some cluster", false
  task :match_mutation_signatures => :tsv do |place_all|
    signatures = step(:merge_signatures).load

    set_info :clones, signatures.fields

    dumper = TSV::Dumper.new :key_field => "Genomic Mutation", :fields => %w(Clones), :type => :flat
    dumper.init

    context_io = TSV.get_stream step(:context)
    io = Misc.sort_genomic_locations_strict(context_io)
    TSV.traverse io, :into => dumper, :bar => self.progress_bar("Matching mutation signatures") do |mutation,context|
      mutation = mutation.first if Array === mutation
      change = mutation.split(":")[2]

      if context[1] == "G" or context[1] == "A"
        #context = Bio::Sequence::NA.new(context).complement.upcase.reverse
        #change = Bio::Sequence::NA.new(change).complement.upcase
        context = context.chars.values_at(2,1,0).collect{|b| Misc::BASE2COMPLEMENT[b] } 
        change = Misc::BASE2COMPLEMENT[change]
      end

      context = "#{context[0]}[#{context[1]}>#{change}]#{context[2]}"

      context_probs = signatures[context]

      next if context_probs.nil?

      context_probs = [context_probs] unless Array === context_probs

      if place_all
        matches = []
        max = context_probs.max
        context_probs.each_with_index{|p,i| matches << i if p == max }
      else
        matches = []
        context_probs.each_with_index{|p,i| matches << i if p >= rand }

        next if matches.empty?
      end

      [mutation, matches]
    end
    dumper
  end

  dep :match_mutation_signatures, :compute => :produce
  input :clone_chr_mutation_counts, :array, "Number of mutations per chromosome for each clone", nil
  input :mutations_per_MB, :float, "Density of mutations in mutations per MB", 5
  task :choose_signature_mutations => :tsv do |clone_chr_mutations,mutations_per_MB|

    clone_names = step(:match_mutation_signatures).info[:clones]
    clone_chr_mutations = [100] * clone_names.length if clone_chr_mutations.nil? || clone_chr_mutations.empty?
    clone_chr_mutations = clone_chr_mutations.collect{|n| n.to_i }

    min_pos = clone_chr_mutations.collect{|n| (1_000_000 * n) / (mutations_per_MB / clone_names.length)  }

    clone_muts = TSV.setup({}, :key_field => "Clone", :fields => ["Genomic Mutation"], :type => :flat)

    last = nil
    chr_clone_muts = nil
    TSV.traverse step(:match_mutation_signatures), :bar => self.progress_bar("Choosing clone mutations") do |mutation,clones|
      mutation = mutation.first if Array === mutation
      chr, pos, change = mutation.split(":")

      last ||= chr
      chr_clone_muts ||= {}

      if last != chr
        clone_names.zip(clone_chr_mutations) do |clone_name,lim|
          chr_clone_muts[clone_name] = chr_clone_muts[clone_name].shuffle[0..lim-1] if chr_clone_muts[clone_name].length > lim
        end

        chr_clone_muts.each do |clone_name, muts|
          clone_muts[clone_name] ||= []
          clone_muts[clone_name].concat muts
        end

        last, chr_clone_muts = chr, {}, {}
      end

      clones_missing_mutations = clone_names.zip(clone_chr_mutations, min_pos).
        select{|name,lim,m_pos| 
          muts = chr_clone_muts[name] || []
          last_pos = muts.last.split(":")[1].to_i if muts.any?
          muts.nil? || muts.length < lim || (last_pos || 0) < m_pos
        }

      next unless clones_missing_mutations.any?

      clone = clones.shuffle.first.to_i
      clone_name = clone_names[clone]

      chr_clone_muts[clone_name] ||= []
      chr_clone_muts[clone_name] << mutation

    end

    clone_names.zip(clone_chr_mutations) do |clone_name,lim|
      chr_clone_muts[clone_name] = chr_clone_muts[clone_name].shuffle[0..lim-1] if chr_clone_muts[clone_name].length > lim
    end

    chr_clone_muts.each do |clone_name, muts|
      clone_muts[clone_name] ||= []
      clone_muts[clone_name].concat muts
    end

    clone_muts
  end
end

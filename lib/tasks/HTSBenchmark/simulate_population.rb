module HTSBenchmark

  dep :genotype_somatic_hg38
  dep :SV_somatic_hg38_PCAWG
  dep :miniref_sizes, :mutations => :genotype_somatic_hg38
  input :parents, :array, "Array of clone parents of all clones except founder", [0,0,1,2,3]
  input :fractions, :array, "Array of clone cell fractions in population", [0.1,0.2,0.3,0.4,0.5,0.6]
  input :mutation_fractions, :array, "Array of fractions of mutations assigned to each clone", [0.1,0.2,0.3,0.4,0.5,0.6]
  input :SV_fractions, :array, "Array of fractions of SVs assigned to each clone", [0.1,0.2,0.3,0.4,0.5,0.6]
  task :simulate_evolution => :yaml do |parents, fractions, mut_fractions, sv_fractions|
    num_clones = fractions.length

    parents = [nil] + parents if parents.length < fractions.length

    sizes = step(:miniref_sizes).load

    mutations = step(:genotype_somatic_hg38).load.shuffle

    svs = step(:SV_somatic_hg38_PCAWG).load

    clean_svs = []
    svs.values.shuffle.each do |values|
      type, chr, start, eend, target_chr, target_start, target_end = values
      target_chr = nil if target_chr.empty?
      target_end = nil if target_end.empty?
      next if sizes[chr].nil?
      next if sizes[chr] < eend.to_i
      next if target_chr && sizes[target_chr] < target_start.to_i
      next if target_chr && target_end && sizes[target_chr] < target_end.to_i

      clean_svs << values
    end

    mutations = mutations.select do |mut|
      chr, pos, alt = mut.split(":")
      pos.to_i <= sizes[chr]
    end

    parents = parents.collect{|p| String === p && p =~ /^\d+$/ ?  p.to_i : p }

    fractions = fractions.collect{|f| f.to_f }
    mut_fractions = mut_fractions.collect{|f| f.to_f }
    sv_fractions = sv_fractions.collect{|f| f.to_f }
    
    total_fractions = Misc.sum(fractions)
    fractions = fractions.collect{|f| f / total_fractions }

    total_mut_fractions = Misc.sum(mut_fractions)
    mut_fractions = mut_fractions.collect{|f| f / total_mut_fractions }

    sv_fractions = sv_fractions.collect{|v| v.to_f}
    total_sv_fractions = Misc.sum(sv_fractions)
    sv_fractions = sv_fractions.collect{|f| f / total_sv_fractions } unless total_sv_fractions == 0.0

    total_mutations = mutations.length
    total_svs = clean_svs.length
    current_mut = current_sv = 0
    evolution = []

    parents.zip(fractions, mut_fractions, sv_fractions).each do |parent,fraction,mut_fraction,sv_fraction|
      if sv_fraction > 0
        final_sv = current_sv + (total_svs * sv_fraction).ceil
        csvs = clean_svs[current_sv..final_sv-1]
        current_sv = final_sv
      else
        csvs = []
      end

      final_mut = current_mut + (total_mutations * mut_fraction).ceil
      cmuts = mutations[current_mut..final_mut-1]
      current_mut = final_mut

      evolution << {"parent": parent, "fraction": fraction, "mutations": cmuts, "SVs": csvs}
    end

    evolution
  end

  dep :simulate_evolution, :compute => :produce
  dep_task :simulate_population, HTSBenchmark, :contaminated_population, :evolution => :placeholder do |jobname,options,dependencies|
    simevo = dependencies.flatten.first
    simevo.produce
    evo = simevo.load
    {:inputs => options.merge(:evolution => evo.to_yaml)}
  end
end

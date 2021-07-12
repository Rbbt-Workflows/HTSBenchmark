module HTSBenchmark

  input :evolution, :text, "Evolution"
  input :germline, :array, "Germline mutations"
  task :population_genotypes => :array do |evolution,germline|
    evolution = YAML.load(evolution)

    ancestry = {}
    evolution.each_with_index do |info,i|
      id = info["id"] || "clone-#{i}"
      clone_dir = file("clone_#{i}")

      if parent = info["parent"]
        if String === parent
          parent_info = evolution.select{|info| info["id"] == parent }.first
        else
          parent_info = evolution[parent]
        end

        parent_number = evolution.index parent_info

        ancestry[i] ||= ancestry[parent_number] + [parent_number]

        Open.write(clone_dir.parent, parent_number.to_s)
      else
        ancestry[i] = []
        Open.write(clone_dir.parent, "none")
      end
    end

    germline_mutations = germline.collect{|mutation| HTSBenchmark.haploid_mutation(mutation) }

    all_mutations = []
    evolution.each_with_index do |info,i|
      clone_dir = file("clone_#{i}")
      next unless info["fraction"].to_f > 0

      mutations = info["mutations"]

      mutations = mutations.collect{|mutation| HTSBenchmark.haploid_mutation(mutation) }

      all_mutations += mutations
      Open.write(clone_dir.mutations, mutations * "\n")

      svs = TSV.setup([], :key_field => "SV ID", :fields => ["Type", "Chromosome", "Start", "End", "Target chromosome", "Target start", "Target end"], :type => :list)

      info["SVs"].each_with_index do |values,i|
        values = HTSBenchmark.haploid_SV(values)
        id = Misc.digest(values.compact * ":" + ".#{i}")
        svs[id] = values
      end if info["SVs"]

      svs.each do |id,values|
        type, chr, start, eend, target_chr, target_start, target_end = values
        all_mutations << [chr, start.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-source"] * ":"
        all_mutations << [chr, eend.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-source-end"] * ":"
        target_chr = chr if target_chr == 'same' || target_chr == 'cis'
        all_mutations << [target_chr, target_start.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-target"] * ":" if target_start
      end

      Open.write(clone_dir.SVs, svs.to_s)
    end

    Open.write(file("total_mutations"), all_mutations.compact.uniq * "\n")
    Open.write(file("total_mutations_clean"), all_mutations.uniq.compact.collect{|m| m.sub(/^copy-\d+_chr/, '') } * "\n")

    evolution.length.times do |clone_number|
      clone_ancestry = ancestry[clone_number]
      clone_svs = TSV.open(file("clone_#{clone_number}").SVs)
      clone_mutations = Open.read(file("clone_#{clone_number}").mutations).split("\n")

      all_mutations = clone_mutations
      clone_ancestry.each do |ancestor_number|
        ancestor_svs = TSV.open(file("clone_#{ancestor_number}").SVs)
        ancestor_mutations = TSV.open(file("clone_#{ancestor_number}").mutations)
        mutation_translations = HTSBenchmark.transpose_mutations(ancestor_svs, all_mutations)
        all_mutations = mutation_translations.values.flatten.uniq
        mutation_translations = HTSBenchmark.transpose_mutations(ancestor_svs, ancestor_mutations)
        all_mutations += mutation_translations.values.collect{|v| v.shuffle.first }
        all_mutations = all_mutations.uniq
        clone_svs = HTSBenchmark.transpose_SVs(ancestor_svs, clone_svs, false)

        germline_mutation_translations = HTSBenchmark.transpose_mutations(ancestor_svs, germline_mutations)
        germline_mutations += germline_mutation_translations.values.collect{|v| v.shuffle.first }
      end

      all_mutations.uniq!

      Open.write(file("clone_#{clone_number}").all_mutations, all_mutations * "\n")
      Open.write(file("clone_#{clone_number}").germline_mutations, germline_mutations * "\n")
      Open.write(file("clone_#{clone_number}").transposed_SVs, clone_svs.to_s)
    end

    Dir.glob(files_dir + "/clone_*")
  end

  dep :genotype_germline_hg38, :jobname => "Default", :compute => :produce
  dep :population_genotypes, :compute => :produce, :germline => :genotype_germline_hg38
  dep :miniref_sizes, :mutations => :placeholder do |jobname,options,dependencies|
    population_genotypes = dependencies.select{|d| d.task_name == :population_genotypes }.first
    total = population_genotypes.file('total_mutations_clean')
    {:inputs => options.merge(:mutations => total)}
  end
  dep :simulate_normal, :mutations => :genotype_germline_hg38, "HTSBenchmark#miniref_sizes" => :miniref_sizes, :jobname => "Default"
  dep :clone, :germline => :placeholder, :somatic => :placeholder, :reference => :placeholder, :compute => :produce, :depth => 90 do |jobname,options,dependencies|
    normal = dependencies.flatten.select{|dep| dep.task_name === :simulate_normal }.first
    population_genotypes = dependencies.flatten.select{|dep| dep.task_name === :population_genotypes }.first
    germline = dependencies.flatten.select{|dep| dep.task_name === :genotype_germline_hg38 }.first
    evolution = YAML.load(options[:evolution])

    clone_jobs = []
    evolution.each_with_index do |info,i|
      id = info["id"] || "clone-#{i}"

      if parent = info["parent"]
        if String === parent
          parent_info = evolution.select{|info| info["id"] == parent }.first
        else
          parent_info = evolution[parent]
        end

        parent_number = evolution.index parent_info

        parent_job = clone_jobs[parent_number]
        clone_jobs 

        parent_reference = parent_job.step(:SV_reference)
        parent_somatic, parent_germline = parent_job.rec_dependencies.select{|d| d.task_name == :SV_mutations }

        clone_options = options.merge(
          :reference => parent_reference,
          :somatic => population_genotypes.file("clone_#{i}").all_mutations, 
          :germline => population_genotypes.file("clone_#{i}").germline_mutations, 
          :svs => population_genotypes.file("clone_#{i}").transposed_SVs, 
        )
        job = HTSBenchmark.job(:clone, jobname + "-clone_#{i}", clone_options)
      else
        clone_options = options.merge(
          :reference => normal.step(:miniref_ploidy), 
          :somatic => population_genotypes.file("clone_#{i}").all_mutations, 
          :germline => population_genotypes.file("clone_#{i}").germline_mutations, 
          :svs => population_genotypes.file("clone_#{i}").transposed_SVs,
        )
        job = HTSBenchmark.job(:clone, jobname + "-clone_#{i}", clone_options)
      end
      clone_jobs << job
    end
    clone_jobs
  end
  dep :merge_clones do |jobname,options,dependencies|
    evolution = YAML.load(options[:evolution])
    samples = dependencies.flatten.compact.select{|d| d.task_name.to_s == 'clone' }
    fractions = evolution.collect{|info| info["fraction"].to_f }
    samples.each do |s| s.init_info end
    {:inputs => {:clones => samples.collect{|s| s.path }, :fractions => fractions } }
  end
  task :population => :array do
    Step.wait_for_jobs dependencies
    normal = step(:simulate_normal)
    tumor = step(:merge_clones)

    nr1, nr2 = normal.load
    tr1, tr2 = tumor.load

    Open.link(nr1, file('normal_read1.fq.gz'))
    Open.link(nr2, file('normal_read2.fq.gz'))
    Open.link(tr1, file('tumor_read1.fq.gz'))
    Open.link(tr1, file('tumor_read2.fq.gz'))

    Dir.glob(files_dir + "/*.fq.gz")
  end

end

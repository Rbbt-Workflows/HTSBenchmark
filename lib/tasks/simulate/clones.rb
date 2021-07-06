
module HTSBenchmark

  input :clones, :array, "Array of NEAT job paths", nil, :nofile => true
  input :fractions, :array, "Array of clone cellular fraction", nil, :nofile => true
  task :NEAT_merge_clones => :text do |clones, fractions|

    Open.mkdir files_dir
    ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each do |filename|
      output = file(filename)

      sout = Misc.open_pipe false, false do |sin|

        clones.zip(fractions).each_with_index do |v,i|
          clone, fraction = v
          fraction = fraction.to_f
          clone_step = Step.new clone 

          input = clone_step.file('output')[filename]


          skip = nil
          rnd = Random.new 1234
          TSV.traverse input, :type => :array, :bar => "Processing #{ Misc.fingerprint [clone_step, filename] }" do |line|
            if line =~ /^@(?:tumor|normal)/
              skip = rnd.rand > fraction 
              next if skip
              line = "@clone_#{i}-" + line[1..-1]
            else
              next if skip
            end

            sin.puts line
          end
        end

        sin.close

      end

      CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
    end
    Dir.glob(files_dir + "/*")
  end

  helper :subclones do |mutations,driver,info|
    genotype = [driver]
    info["num"].to_i.times do 
      genotype << mutations.gets.strip
    end

    clones = [genotype]
    info.each do |driver,subclones|
      next if driver == "num"
      clones += subclones(mutations, driver, subclones).collect{|list| genotype + list }
    end

    clones
  end


  input :evolution, :text, "Evolution"
  task :evolution => :text do |evolution|
    evolution = YAML.load(evolution)

    evolution.each_with_index do |info,i|
      id = info["id"] || "clone-#{i}"
      real_mutations = info["mutations"] || []

      if parent = info["parent"]
        if String === parent
          parent_info = evolution.select{|info| info["id"] == parent }.first
        else
          parent_info = evolution[parent]
        end

        real_mutations += parent_info["real_mutations"]
      end

      info["real_mutations"] = real_mutations
    end

    all_mutations = []
    evolution.each_with_index do |info,i|
      next unless info["fraction"].to_f > 0

      mutations = info["mutations"]
      all_mutations += mutations
      Open.write(file("mutations/clone_#{i}_mutations"), mutations * "\n")

      svs = TSV.setup(info["SVs"] || [], :type => :list)

      svs.each do |id,values|
        type, chr, start, eend, target_chr, target_start, target_end = values
        all_mutations << [chr, start.to_i, "SV", Misc.digest(values) + '-source'] * ":"
        all_mutations << [target_chr, target_start.to_i, "SV", Misc.digest(values) + '-target'] * ":" if target_start
      end
      Open.write(file("mutations/clone_#{i}_SVs"), svs.to_s)
    end

    Open.write(file("mutations/total"), all_mutations.uniq * "\n")

    file('mutations').glob("*")
  end

  dep :evolution, :compute => :produce
  dep_task :total_vcf, Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname, options, dependencies|
    dep = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution'}.first
    {:inputs => {:positions => dep.file('mutations/total') } }
  end

  dep :parent_clone do
    []
  end
  input :mutations, :array
  task :transpose_mutations => :array do |mutations|
    if dependencies.empty?
      mutations
    else
      mutation_translations = step(:parent_clone).file('mutation_translations.tsv').tsv
      TSV.traverse mutations, :type => :array do |mutation|
        mutation_translation[mutation]
      end
    end
  end

  dep :parent_clone do [] end
  input :SVs, :tsv
  task :transpose_SVs => :array do |svs|
    if dependencies.empty?
      mutations
    else
      mutation_translations = step(:parent_clone).file('mutation_translations.tsv').tsv
      new = svs.annotate({})
      TSV.traverse svs do |id,values|
        type, chr, start, eend, target_chr, target_start, target_end = values

        if mutation_translations.include?(id + '-source')
          new_source_chr, new_source_start = mutation_translations[id + "-source"]
          new_source_start = new_source_start.to_i
          new_source_end = eend.to_i + (new_source_start - start.to_i)
          chr = new_source_chr
          start = new_source_start
          eend = new_source_eend
        end

        if mutation_translations.include?(id + '-target')
          new_target_chr, new_target_start = mutation_translations[id + "-target"]
          new_target_start = new_target_start.to_i
          new_target_end = eend.to_i + (new_target_start - start.to_i)
          chr = new_target_chr
          start = new_target_start
          eend = new_target_eend
        end

        new[id] = [type, chr, start, eend, target_chr, target_start, target_end]
      end
      new
    end
  end


  dep :evolution
  dep :parent_clone do [] end
  input :clone_number, :integer
  dep :transpose_mutations, :mutations => :placeholder do |jobname,options,dependencies|
    evolution = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    clone_number = options[:clone_number]
    ['somatic', 'germline'].collect do |type|
      file = evolution.file("clone_#{clone_number}_mutations")
      {:inputs => options.merge(:mutations => file), :jobname => jobname + "-#{type}"}
    end
  end
  dep :transpose_SVs, :SVs => :placeholder do |jobname,options,dependencies|
    evolution = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    clone_number = options[:clone_number]
    file = evolution.file("clone_#{clone_number}_SVs")
    {:inputs => options.merge(:mutations => file), :jobname => jobname}
  end
  dep :total_vcf
  dep :miniref_sizes, :vcf => :total_vcf
  dep_task :clone, HTSBenchmark, :simulated_SV_sample, "HTSBenchmark#miniref_sizes" => :miniref_sizes do |jobname,options,dependencies|
    clone_somatic_mutations, clone_germline_mutations = dependencies.flatten.select{|d| d.task_name.to_s == 'transpose_mutations' }
    clone_SVs = dependencies.flatten.select{|d| d.task_name.to_s == 'transpose_SVs' }.first

    parent = dependencies.flatten.select{|d| d.task_name.to_s == 'parent_clone' }.first
    parent_miniref = parent.step(:SV_miniref) if parent

    inputs = options.dup
    inputs["SVs"] = clone_SVs
    inputs["HTSBenchmark#simulate_somatic_hg38"] = clone_somatic_mutations
    inputs["HTSBenchmark#simulate_germline_hg38"] = clone_germline_mutations
    inputs["HTSBenchmark#miniref_ploidy"] = parent_miniref if parent_miniref
    {:inputs => inputs} 
  end
  dep_task :parent_clone, HTSBenchmark, :clone

  dep :evolution, :compute => :produce
  dep :total_vcf
  dep :miniref_sizes, :vcf => :total_vcf
  dep :clone, :compute => [:bootstrap, 10], :clone_number => :placeholder do |jobname,options,dependencies|
    evolution_dep = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    miniref_sizes = dependencies.flatten.select{|d| d.task_name.to_s == 'miniref_sizes' }.first
    evolution = options[:evolution]
    evolution = YAML.load(evolution)
    num_clones = evolution.length
    options[:all_mutations] = evolution_dep.file('total_mutations')

    clone_deps = []
    evolution.each_with_index do |clone_info,clone_number| 
      clone_inputs = options.dup
      clone_inputs[:clone_number] = clone_number
      parent = clone_info[:parent]
      parent = parent.to_i - 1 unless parent.nil?
      mutation_file = evolution_dep.file("mutations/clone_#{clone_number}_mutations")
      clone_inputs["HTSBenchmark#simulate_somatic_hg38"] = mutation_file
      clone_inputs["HTSBenchmark#miniref_sizes"] = miniref_sizes
      clone_inputs["HTSBenchmark#parent_clone"] = clone_deps[parent] if parent && clone_deps[parent]
      clone_name = jobname + "-clone_#{clone_number}"
      clone_job = HTSBenchmark.job(:clone, clone_name, clone_inputs)
      clone_deps << clone_job
    end

    clone_deps
  end
  dep_task :population, HTSBenchmark, :NEAT_merge_clones, :clones => :placeholder, :fractions => :placeholder do |jobname,options,dependencies|
    evolution = YAML.load(options[:evolution])
    samples = dependencies.flatten.compact.select{|d| d.task_name.to_s == 'clone' }
    indices = samples.collect{|d| d.step(:simulate_somatic_hg38).clean_name.split("-").last.scan(/\d+/).first.to_i }
    fractions = evolution.values_at(*indices).collect{|info| info["fraction"].to_f }
    {:inputs => {:clones => samples.collect{|s| s.path }, :fractions => fractions } }
  end

  #######


  dep :evolution, :compute => :produce
  dep :total_vcf
  dep :miniref_sizes, :vcf => :total_vcf
  dep :simulated_sample, :compute => [:bootstrap, 10] do |jobname,options,dependencies|
    dep = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    miniref_sizes = dependencies.flatten.select{|d| d.task_name.to_s == 'miniref_sizes' }.first
    dep.file('mutations').glob("clone_*").collect do |mutation_file| 
      inputs = options.dup
      inputs["HTSBenchmark#simulate_somatic_hg38"] = mutation_file
      inputs["HTSBenchmark#miniref_sizes"] = miniref_sizes
      {:inputs => inputs} 
    end
  end
  dep_task :clonal, HTSBenchmark, :NEAT_merge_clones, :clones => :placeholder, :fractions => :placeholder do |jobname,options,dependencies|
    evolution = YAML.load(options[:evolution])
    samples = dependencies.flatten.compact.select{|d| d.task_name.to_s == 'simulated_sample' }
    indices = samples.collect{|d| d.step(:simulate_somatic_hg38).clean_name.split("-").last.scan(/\d+/).first.to_i }
    fractions = evolution.values_at(*indices).collect{|info| info["fraction"].to_f }
    {:inputs => {:clones => samples.collect{|s| s.path }, :fractions => fractions } }
  end

  dep :simulate_somatic_hg38
  input :cell_fractions, :array, "Clonal fractions", [0.6, 0.2, 0.2]
  input :mutation_fractions, :array, "Clonal fraction of mutations", [0.2, 0.4, 0.4]
  input :clone_parents, :array, "Clone parents", [1, 2]
  task :simulate_evolution => :yaml do |cf, mf, cp|
    evolution = []
    cf.zip([nil] + cp).each do |_cf,_cp|
      evolution << {"fraction" => _cf, "parent" => _cp, "mutations" => []}
    end

    acc = [0]
    current = 0
    mf.each do |f|
      current += f.to_f
      acc << current
    end

    TSV.traverse step(:simulate_somatic_hg38), :type => :array do |mutation|
      r = rand
      acc.each_with_index do |v,i|
        if v > r
          evolution[i]["mutations"] << mutation
          break
        end
      end
    end

    evolution
  end

  dep :simulate_evolution
  dep_task :simulate_clonal, HTSBenchmark, :clonal, :evolution => :placeholder do |jobname,inputs,dependencies|
    inputs["evolution"] = dependencies.flatten.first.path.read
    {:inputs => inputs}
  end
end

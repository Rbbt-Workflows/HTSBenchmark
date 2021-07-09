
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
          clone_step = Step === clone ? clone : Step.new(clone)

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

      svs = TSV.setup([], :key_field => "SV ID", :fields => ["Type", "Chromosome", "Start", "End", "Target chromosome", "Target start", "Target end"], :type => :list)

      info["SVs"].each_with_index do |values,i|
        id = Misc.digest(values.compact * ":" + ".#{i}")
        svs[id] = values
      end if info["SVs"]

      svs.each do |id,values|
        type, chr, start, eend, target_chr, target_start, target_end = values
        all_mutations << [chr, start.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-source"] * ":"
        target_chr = chr if target_chr == 'same' || target_chr == 'cis'
        all_mutations << [target_chr, target_start.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-target"] * ":" if target_start
      end

      Open.write(file("mutations/clone_#{i}_SVs"), svs.to_s)
    end

    Open.write(file("mutations/total"), all_mutations.uniq * "\n")
    Open.write(file("mutations/total_clean"), all_mutations.uniq.collect{|m| m.sub(/^copy-\d+_chr/, '') } * "\n")

    file('mutations').glob("*")
  end

  dep :evolution, :compute => :produce, :jobname => "Default"
  dep_task :total_vcf, Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname, options, dependencies|
    dep = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution'}.first
    {:inputs => {:positions => dep.file('mutations/total_clean') } }
  end

  dep :parent_clone do
    []
  end
  input :mutations, :array
  input :add_somatic, :boolean, "Add somatic mutations from parent clone", false
  input :just_mutations, :boolean, "Only return mutations without reference", false
  task :transpose_mutations => :tsv do |mutations,add_somatic,just_mutations|
    
    mutations = Open.read(mutations).split("\n") unless Array === mutations

    ancestry = []
    ancestor = dependencies.first
    while ancestor
      ancestry << ancestor
      begin
        ancestor = ancestor.step(:parent_clone)
      rescue
        ancestor = nil
      end
    end

    while ancestor = ancestry.pop
      mutation_translations = ancestor.step(:add_SV_to_reference).file('mutation_translations.tsv').tsv

      mutations = TSV.traverse mutations, :type => :array, :into => [] do |mutation|
        res = mutation_translations[mutation]
        res = res.shuffle.first
        Log.warn mutation if res.nil?
        next if res.nil?
        res
      end

    end

    if dependencies.any?
      if add_somatic
        mutations += dependencies.first.step(:add_SV_to_reference).file('somatic').list 
      end

      if just_mutations
        Open.write(self.tmp_path, mutations.uniq * "\n")
        nil
      else
        mut_ref =  dependencies.first.step(:add_SV_to_reference).file('all_mutations.reference').tsv 
        mut_ref.select(mutations.uniq)
      end
    else
      if just_mutations
        Open.write(self.tmp_path, mutations.uniq * "\n")
        nil
      else
        translations = {}
        mutations.uniq.each do |m|
          mf = m.sub(/^copy-\d+_chr/,'')
          translations[mf] = m
        end
        tsv = Sequence.job(:reference, nil, :positions => translations.keys, :organism => HG38_ORGANISM).exec
        new = tsv.annotate({})
        tsv.each do |k,v|
          k = translations[k]
          new[k] = v
        end
        new
      end
    end
  end

  dep :transpose_mutations
  dep_task :transpose_mutations_vcf, Sequence, :mutations_to_vcf, "Sequence#reference" => :transpose_mutations

  dep :parent_clone, :compute => :produce do [] end
  input :SVs, :tsv
  task :transpose_SVs => :tsv do |svs|
    svs = TSV.open(svs) unless TSV === svs
    ancestry = []
    ancestor = dependencies.first
    while ancestor
      ancestry << ancestor
      begin
        ancestor = ancestor.step(:parent_clone)
      rescue
        ancestor = nil
      end
    end

    while ancestor = ancestry.pop
      mutation_translations = ancestor.step(:add_SV_to_reference).file('mutation_translations.tsv').tsv
      new_svs = svs.annotate({})
      TSV.traverse svs do |id,values|
        type, chr, start, eend, target_chr, target_start, target_end = values

        target_chr = chr if target_chr == 'cis' || target_chr == 'same'

        orig_chr = chr
        orig_target_chr = target_chr

        #chr = chr.sub(/^copy-\d+_chr/, '')
        #target_chr = target_chr.sub(/^copy-\d+_chr/, '')

        source_id = [chr, start,"SV",id] * ":" + '-source'
        if mutation_translations.include?(source_id)
          new_list = mutation_translations[source_id]
          new_list_match = new_list.select{|chr,source| chr == orig_chr }
          new_list = new_list_match if new_list_match.any?

          new = new_list.shuffle.first
          new_source_chr, new_source_start = new
          new_source_start = new_source_start.to_i
          new_source_end = eend.to_i + (new_source_start - start.to_i)
          chr = new_source_chr
          start = new_source_start
          eend = new_source_end
        end

        target_id = [target_chr, target_start,"SV",id] * ":" + '-target'
        if mutation_translations.include?(target_id)
          new_list = mutation_translations[target_id]
          new_list_match = new_list.select{|chr,source| chr == orig_target_chr }
          new_list = new_list_match if new_list_match.any?

          new = new_list.shuffle.first
          new_target_chr, new_target_start = new
          new_target_start = new_target_start.to_i
          new_target_end = eend.to_i + (new_target_start - start.to_i)
          chr = new_target_chr
          start = new_target_start
          eend = new_target_end
        end

        new_svs[id] = [type, chr, start, eend, target_chr, target_start, target_end]
      end
      svs = new_svs
    end

    svs
  end

  # CLONE

  dep :evolution, :jobname => "Default", :compute => :produce
  dep :parent_clone do [] end
  input :clone_number, :integer
  dep :transpose_mutations_vcf, :mutations => :placeholder do |jobname,options,dependencies|
    evolution = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    clone_number = options[:clone_number]
    file = evolution.file("mutations/clone_#{clone_number}_mutations")
    {:inputs => options.merge(:mutations => file), :jobname => jobname + "-somatic" }
  end
  dep :simulate_germline_hg38, :jobname => "Default"
  dep :transpose_mutations_vcf, :mutations => :simulate_germline_hg38, :add_somatic => false
  dep :transpose_mutations, :mutations => :placeholder, :add_somatic => true, :just_mutations => true do |jobname,options,dependencies|
    evolution = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    clone_number = options[:clone_number]
    file = evolution.file("mutations/total")
    {:inputs => options.merge(:mutations => file), :jobname => jobname + "-total"}
  end
  dep :transpose_SVs, :SVs => :placeholder do |jobname,options,dependencies|
    evolution = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    clone_number = options[:clone_number]
    file = evolution.file("mutations/clone_#{clone_number}_SVs")
    {:inputs => options.merge(:SVs => file), :jobname => jobname}
  end
  dep :total_vcf
  dep :miniref_sizes, :vcf => :total_vcf
  dep_task :clone, HTSBenchmark, :simulated_SV_sample, "HTSBenchmark#miniref_sizes" => :miniref_sizes, :haploid => true, :only_tumor => true do |jobname,options,dependencies|
    clone_somatic_mutations, clone_germline_mutations = dependencies.flatten.select{|d| d.task_name.to_s == 'transpose_mutations_vcf' }
    all_mutations = dependencies.flatten.select{|d| d.task_name.to_s == 'transpose_mutations' }.first
    clone_SVs = dependencies.flatten.select{|d| d.task_name.to_s == 'transpose_SVs' }.last

    parent = dependencies.flatten.select{|d| d.task_name.to_s == 'parent_clone' }.first

    parent_miniref = parent.step(:SV_miniref) if parent

    inputs = options.dup
    inputs["SVs"] = clone_SVs
    inputs["all_mutations"] = all_mutations
    inputs["HTSBenchmark#simulate_somatic_hg38_vcf"] = clone_somatic_mutations
    inputs["HTSBenchmark#simulate_germline_hg38_vcf"] = clone_germline_mutations
    inputs["HTSBenchmark#miniref_ploidy"] = parent_miniref if parent_miniref
    {:inputs => inputs} 
  end
  dep_task :parent_clone, HTSBenchmark, :clone


  # POPULATION

  dep :evolution, :compute => :produce, :jobname => "Default"
  dep :total_vcf
  dep :miniref_sizes, :vcf => :total_vcf
  dep :simulate_germline_hg38_vcf, :jobname => "Default"
  dep :simulate_somatic_hg38_vcf
  dep :miniref_ploidy, :vcf => :total_vcf, :compute => :produce, :jobname => 'hg38'
  dep :minify_vcf, :vcf_file => :simulate_germline_hg38_vcf, :sizes => :placeholder, :jobname => "germline" do |jobname,options,dependencies|
    miniref_sizes = dependencies.flatten.select{|d| d.task_name.to_s == 'miniref_sizes' }.first
    {:inputs => options.merge(:sizes => miniref_sizes), :jobname => jobname }
  end
  dep :vcf_ploidy, :haploid_reference => :miniref_ploidy do |jobname,options,dependencies|
    if ! options[:only_tumor]
      germline = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}.first
      {:inputs => options.merge(:vcf => germline.path), :jobname => germline.clean_name}
    else
      []
    end
  end
  dep HTSBenchmark, :NEAT_genreads, :somatic => :placeholder, :germline => :placeholder, :reference => :miniref_ploidy, :produce_samples => :normal, :haploid => true do |jobname,options,dependencies|
    if ! options[:only_tumor]
      germline = dependencies.flatten.select{|dep| dep.task_name.to_s == 'vcf_ploidy'}.first
      options[:germline] = germline
      options[:somatic] = nil
      {:inputs => options, :jobname => jobname + '-normal'}
    else
      []
    end
  end
  dep :clone, :compute => [:bootstrap, 10], :clone_number => :placeholder do |jobname,options,dependencies|
    evolution_dep = dependencies.flatten.select{|d| d.task_name.to_s == 'evolution' }.first
    miniref_sizes = dependencies.flatten.select{|d| d.task_name.to_s == 'miniref_sizes' }.first
    evolution = options[:evolution]
    evolution = YAML.load(evolution)
    num_clones = evolution.length

    clone_deps = []
    clone_names = {}
    evolution.each_with_index do |clone_info,clone_number| 
      clone_inputs = IndiferentHash.setup(options.dup)
      clone_inputs[:all_mutations] = evolution_dep.file('mutations/total')
      clone_inputs[:clone_number] = clone_number
      name = clone_info["id"]
      clone_names[name] = clone_number if name
      parent = clone_info["parent"]
      parent = clone_names[parent] if clone_names.include?(parent)
      parent = parent.to_i unless parent.nil?
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
  dep HTSBenchmark, :NEAT_merge_clones, :clones => :placeholder, :fractions => :placeholder do |jobname,options,dependencies|
    evolution = YAML.load(options[:evolution])
    samples = dependencies.flatten.compact.select{|d| d.task_name.to_s == 'clone' || d.task_name.to_s == 'parent_clone' }
    samples = samples.uniq
    indices = samples.collect{|d| d.step(:simulate_somatic_hg38_vcf).clean_name.split("-").last.scan(/\d+/).first.to_i }
    fractions = evolution.values_at(*indices).collect{|info| info["fraction"].to_f }
    samples.each do |s| s.init_info end
    {:inputs => {:clones => samples.collect{|s| s.path }, :fractions => fractions } }
  end
  task :population => :text do
    normal = step(:NEAT_genreads)
    tumor = step(:NEAT_merge_clones)

    Open.link(normal.file('output/normal_read1.fq.gz'), file('normal_read1.fq.gz'))
    Open.link(normal.file('output/normal_read2.fq.gz'), file('normal_read2.fq.gz'))
    Open.link(tumor.file('tumor_read1.fq.gz'), file('tumor_read1.fq.gz'))
    Open.link(tumor.file('tumor_read2.fq.gz'), file('tumor_read2.fq.gz'))

    Dir.glob(files_dir + "/*.fq.gz") * "\n"
  end

end

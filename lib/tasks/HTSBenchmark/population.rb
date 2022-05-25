module HTSBenchmark

  input :evolution, :text, "Evolution"
  input :germline, :array, "Germline mutations"
  task :population_genotypes => :array do |evolution,germline|
    if Step === evolution
      evolution = evolution.load
    elsif String === evolution
      evolution = YAML.load(evolution)
    else
      evolution = evolution
    end

    log :ancestry, "Preparing ancestry"
    ancestry = {}
    evolution.each_with_index do |info,i|
      IndiferentHash.setup(info)
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

    log :mutations, "Preparing germline mutations"
    germline_mutations = germline.collect{|mutation| HTSBenchmark.haploid_mutation(mutation) }

    log :mutations, "Preparing clone mutations"
    all_mutations = []
    evolution.each_with_index do |info,i|
      IndiferentHash.setup(info)
      clone_dir = file("clone_#{i}")

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
        target_start = nil if target_start == ""
        all_mutations << [target_chr, target_start.to_i, "SV", Misc.digest(values.compact * ":") + ".#{i}-target"] * ":" if target_start 
      end

      Open.write(clone_dir.SVs, svs.to_s)
    end

    Open.write(file("total_mutations"), all_mutations.compact.uniq * "\n")
    Open.write(file("total_mutations_clean"), all_mutations.uniq.compact.collect{|m| m.sub(/^copy-\d+_chr/, '') } * "\n")

    bar = self.progress_bar("Processing clones", :max => evolution.length)
    bar.init
    evolution.length.times do |clone_number|

      clone_svs = TSV.open(file("clone_#{clone_number}").SVs)
      clone_mutations = Open.read(file("clone_#{clone_number}").mutations).split("\n")
      clone_ancestry = ancestry[clone_number]

      ancestry_mutations = []
      clone_transposed_mutations = clone_mutations
      clone_transposed_svs = clone_svs
      clone_ancestry.each do |ancestor_number|

        ancestor_directory = file("clone_#{ancestor_number}")
        ancestor_transposed_svs = ancestor_directory.transposed_SVs.tsv
        ancestor_transposed_mutations = ancestor_directory.transposed_mutations.list

        ancestry_mutations = 
          HTSBenchmark.transpose_mutations(ancestor_transposed_svs, ancestry_mutations).
          values.flatten.compact

        ancestry_mutations += 
          HTSBenchmark.transpose_mutations(ancestor_transposed_svs, ancestor_transposed_mutations).
          values.collect{|v| v.shuffle.first}.compact

        clone_transposed_mutations = 
          HTSBenchmark.transpose_mutations(ancestor_transposed_svs, clone_transposed_mutations).
          values.collect{|v| v.shuffle.first}.compact

        clone_transposed_svs = HTSBenchmark.transpose_SVs(ancestor_transposed_svs, clone_transposed_svs, false)
      end

      if clone_ancestry.any?
        ancestor_directory = file("clone_#{clone_ancestry.last}")
        ancestor_transposed_svs = ancestor_directory.transposed_SVs.tsv
        ancestor_transposed_germline = ancestor_directory.transposed_germline.list
        clone_transposed_germline = HTSBenchmark.transpose_mutations(ancestor_transposed_svs, ancestor_transposed_germline).values.flatten.compact
      else
        clone_transposed_germline = germline_mutations
      end

      Open.write(file("clone_#{clone_number}").transposed_mutations, (ancestry_mutations + clone_transposed_mutations) * "\n")
      Open.write(file("clone_#{clone_number}").transposed_SVs, clone_transposed_svs.to_s)
      Open.write(file("clone_#{clone_number}").transposed_germline, clone_transposed_germline * "\n")

      bar.tick
    end
    bar.remove

    Dir.glob(files_dir + "/clone_*")
  end

  input :bundle, :boolean, "Production run, do VCFs for miniref bundle", false
  input :tumor_depth, :integer, "Depth to sequence tumor", 90
  input :normal_depth, :integer, "Depth to sequence normal", 30
  dep :genotype_germline_hg38, :jobname => "Default", :compute => :produce
  dep :population_genotypes, :compute => :produce, :germline => :genotype_germline_hg38
  dep :miniref_sizes, :mutations => :placeholder do |jobname,options,dependencies|
    population_genotypes = dependencies.select{|d| d.task_name == :population_genotypes }.first
    total = population_genotypes.file('total_mutations_clean')
    {:inputs => options.merge(:mutations => total)}
  end
  dep :simulate_normal, 
    :mutations => :genotype_germline_hg38, 
    "HTSBenchmark#miniref_sizes" => :miniref_sizes, 
    "HTSBenchmark#genotype_somatic_hg38" => :genotype_germline_hg38, 
    :not_overriden => true,
    :depth => :normal_depth,
    :jobname => "Default"
  dep :clone, :germline => :placeholder, :somatic => :placeholder, :svs => :placeholder, :reference => :placeholder, :compute => :produce, :depth => :placeholder do |jobname,options,dependencies|
    normal = dependencies.flatten.select{|dep| dep.task_name === :simulate_normal }.first
    population_genotypes = dependencies.flatten.select{|dep| dep.task_name === :population_genotypes }.first
    germline = dependencies.flatten.select{|dep| dep.task_name === :genotype_germline_hg38 }.first
    evolution = Step === options[:evolution] ? options[:evolution].load : YAML.load(options[:evolution])

    options[:depth] = options[:tumor_depth]

    clone_jobs = []
    evolution.each_with_index do |info,i|
      IndiferentHash.setup(info)
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
          :somatic => population_genotypes.file("clone_#{i}").transposed_mutations, 
          :germline => population_genotypes.file("clone_#{i}").transposed_germline, 
          :svs => population_genotypes.file("clone_#{i}").transposed_SVs, 
        )
        job = HTSBenchmark.job(:clone, jobname + "-clone_#{i}", clone_options)
      else
        clone_options = options.merge(
          :reference => normal.step(:miniref_ploidy), 
          :somatic => population_genotypes.file("clone_#{i}").transposed_mutations, 
          :germline => population_genotypes.file("clone_#{i}").transposed_germline, 
          :svs => population_genotypes.file("clone_#{i}").transposed_SVs,
        )
        job = HTSBenchmark.job(:clone, jobname + "-clone_#{i}", clone_options)
      end
      clone_jobs << job
    end
    clone_jobs
  end
  dep :merge_clones, :clones => :placeholder, :fractions => :placeholder do |jobname,options,dependencies|
    evolution = Step === options[:evolution] ? options[:evolution].load : YAML.load(options[:evolution])
    samples = dependencies.flatten.compact.select{|d| d.task_name.to_s == 'clone' }
    fractions = evolution.collect{|info| (info[:fraction] || info["fraction"]).to_f }
    samples.each do |s| s.init_info end
    {:inputs => {:clones => samples, :fractions => fractions } }
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
    Open.link(tr2, file('tumor_read2.fq.gz'))

    Dir.glob(files_dir + "/*.fq.gz")
  end

  input :tumor_depth, :integer, "Depth to sequence tumor", 90
  input :normal_depth, :integer, "Depth to sequence normal", 30
  dep :population
  dep :simulate_normal, 
    :mutations => :genotype_germline_hg38, 
    "HTSBenchmark#miniref_sizes" => :miniref_sizes, 
    "HTSBenchmark#genotype_somatic_hg38" => :genotype_germline_hg38, 
    :not_overriden => true,
    :depth => :tumor_depth,
    :jobname => "Default" 
  input :normal_in_tumor_contamination, :float, "Proportion of normal contamination in tumor", 0.1
  input :tumor_in_normal_contamination, :float, "Proportion of tumor contamination in normal", 0
  task :contaminated_population => :array do |tumor_depth,normal_depth,normal_in_tumor_contamination, tumor_in_normal_contamination|
    Open.mkdir files_dir

    normal = dependencies.select{|d| d.task_name == :simulate_normal }.first
    tumor = step(:merge_clones)

    if normal_in_tumor_contamination > 0
      ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)

        clones = [normal, tumor]
        fractions = [normal_in_tumor_contamination, 1 - normal_in_tumor_contamination]

        sout = Misc.open_pipe false, false do |sin|

          clones.zip(fractions).each_with_index do |v,ci|
            clone, fraction = v
            fraction = fraction.to_f
            clone_step = Step === clone ? clone : Step.new(clone)

            input = clone_step.load[i]

            skip = nil
            rnd = Random.new 1234
            TSV.traverse input, :type => :array, :bar => "Processing #{ Misc.fingerprint [input, filename] }" do |line|
              if line =~ /^@.*(clone|normal)_/
                rand = rnd.rand
                skip = rand > fraction 
                next if skip
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
    else
      ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)
        Open.link tumor.load[i], output
      end
    end

    if tumor_in_normal_contamination > 0
      ['normal_read1.fq.gz', 'normal_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)

        clones = [tumor, normal]
        fractions = [tumor_in_normal_contamination, (1 - tumor_in_normal_contamination) * (normal_depth.to_f / tumor_depth.to_f)]

        sout = Misc.open_pipe false, false do |sin|

          clones.zip(fractions).each_with_index do |v,ci|
            clone, fraction = v
            fraction = fraction.to_f
            clone_step = Step === clone ? clone : Step.new(clone)

            input = clone_step.load[i]

            skip = nil
            rnd = Random.new 1234
            TSV.traverse input, :type => :array, :bar => "Processing #{ Misc.fingerprint [input, filename] }" do |line|
              if line =~ /^@.*(clone|normal)_/
                rand = rnd.rand
                skip = rand > fraction 
                next if skip
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
    else
      ['normal_read1.fq.gz', 'normal_read2.fq.gz'].each_with_index do |filename,i|
        output = file(filename)
        Open.link normal.load[i], output
      end
    end
    Dir.glob(files_dir + "/*")
  end

end

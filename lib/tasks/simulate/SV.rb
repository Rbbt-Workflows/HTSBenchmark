require 'SV'

module HTSBenchmark

  input :reference, :file, "Reference file", nil, :nofile => true
  input :germline_mutations, :array, "Germline mutations"
  input :somatic_mutations, :array, "Somatic mutations"
  input :SVs, :tsv, "List of SVs"
  input :all_mutations, :array, "All mutations"
  task :add_SV_to_reference => :text do |reference,germline,somatic,svs,all_mutations|

    somatic = somatic.collect{|m| m.start_with?('chr') ? m : 'chr' + m }
    germline = germline.collect{|m| m.start_with?('chr') ? m : 'chr' + m }

    all_mutations = somatic + germline if all_mutations.nil? || all_mutations.empty?

    reference_io, mutation_translations = HTSBenchmark.apply_SVs(reference, all_mutations, svs)

    Open.write(file('mutation_translations.tsv'), mutation_translations.to_s)

    Open.mkdir files_dir 
    target_reference = file(File.basename(reference).split("_").first) + '.fa' 
    Open.write(target_reference, reference_io)
    #CMD.cmd('bgzip', " -c > #{target_reference}", :in => reference_io)

    Open.write(file('somatic'), somatic.collect{|m| mutation_translations[m].collect{|m| m.sub(/^chr/, '') } }.flatten.compact * "\n")
    Open.write(file('germline'), germline.collect{|m| mutation_translations[m].collect{|m| m.sub(/^chr/, '') } }.flatten.compact * "\n")

    CMD.cmd(:samtools, "faidx #{target_reference}")

    ['somatic', 'germline'].each do |type|
      TmpFile.with_file do |ranges|
        mutations = Open.read(file(type)).split("\n")

        mutation_ranges = mutations.collect do |m| 
          chr, pos, alt = m.split(":")
          if alt["-"]
            pos = pos.to_i - 1
            eend = pos + alt.length
          else
            eend = pos
          end
          chr + ":" + pos.to_s + "-" + eend.to_s
        end
        
        Open.write(ranges, mutation_ranges * "\n")

        pos_reference = {}
        TSV.traverse CMD.cmd(:samtools, "faidx #{target_reference} -r #{ranges} 2> /dev/null | tr '\\n' '\\t' | tr '>' '\\n'", :pipe => true), :type => :array do |line|
          pos_info, ref = line.split("\t")
          next if ref.nil?
          chr, range = pos_info.split(":")
          pos = range.split("-").first
          pos_reference[[chr, pos]] = ref
        end

        reference = TSV.setup({}, :key_field => "Genomic Mutation", :fields => ["Reference"], :type => :single)
        mutations.each do |mutation|
          chr, pos, alt = mutation.split(":")
          if alt["-"]
            pos = (pos.to_i + 1).to_s
          end
          ref = pos_reference[[chr, pos]]
          next if ref.nil?
          reference[mutation] = ref
        end

        target = file(type + '.reference')
        Open.write(target, reference.to_s)
      end
    end

    CMD.cmd('bgzip', "#{target_reference}")

    Dir.glob(files_dir + "/*.fa.gz").first
  end

  dep :simulate_germline_hg38
  dep :simulate_somatic_hg38
  dep :simulate_somatic_hg38_vcf
  dep :miniref_ploidy, :vcf => :simulate_somatic_hg38_vcf, :jobname => 'hg38'
  dep_task :SV_miniref, HTSBenchmark, :add_SV_to_reference, :reference => :miniref_ploidy, :germline_mutations => :simulate_germline_hg38, :somatic_mutations => :simulate_somatic_hg38

  dep :SV_miniref
  task :SV_miniref_sizes => :json do
    sizes = {}
    chr = nil
    TSV.traverse Dir.glob(step(:SV_miniref).files_dir + "/*.fa.gz").first, :type => :array do |line|
      if m = line.match(/^>([^\s]*)/)
        chr = m.captures[0]
      else
        sizes[chr] = line.length
      end
    end
    sizes
  end

  dep :SV_miniref, :compute => :produce
  dep :simulate_somatic_hg38_vcf
  dep Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname,options,dependencies|
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.first
    file = sv_miniref.file('germline.reference')
    {:inputs => options.merge("Sequence#reference" => file), :jobname => jobname + '-germline' }
  end
  dep Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname,options,dependencies|
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.first
    file = sv_miniref.file('somatic.reference')
    {:inputs => options.merge("Sequence#reference" => file), :jobname => jobname + '-somatic' }
  end
  dep :SV_miniref_sizes, :SV_miniref => :SV_miniref
  dep :minify_vcf, :vcf_file => :placeholder, :sizes => :SV_miniref_sizes, :jobname => "germline" do |jobname,options,dependencies|
    vcf_file = dependencies.flatten.select{|dep| dep.task_name.to_s == 'mutations_to_vcf' }.first
    {:inputs => options.merge(:vcf_file => vcf_file) }
  end
  dep :minify_vcf, :vcf_file => :placeholder, :sizes => :SV_miniref_sizes, :jobname => "somatic" do |jobname,options,dependencies|
    vcf_file = dependencies.flatten.select{|dep| dep.task_name.to_s == 'mutations_to_vcf' }.last
    {:inputs => options.merge(:vcf_file => vcf_file) }
  end
  dep_task :simulated_SV_sample, HTSBenchmark, :NEAT_genreads, :somatic => :placeholder, :germline => :placeholder, :reference => :placeholder do |jobname,options,dependencies|
    germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.first
    reference = sv_miniref.files_dir + "/hg38.fa.gz"
    options[:reference] = reference
    options[:germline] = germline
    options[:somatic] = somatic
    {:inputs => options, :jobname => jobname}
  end
end

require 'SV'

module HTSBenchmark

  input :reference, :file, "Reference file", nil, :nofile => true
  input :germline_mutations, :array, "Germline mutations"
  input :somatic_mutations, :array, "Somatic mutations"
  input :SVs, :tsv, "List of SVs"
  input :all_mutations, :array, "All mutations"
  extension 'fa.gz'
  task :add_SV_to_reference => :binary do |reference,germline,somatic,svs,all_mutations|

    somatic = somatic.collect{|m| m.start_with?('chr') ? m : 'chr' + m }
    germline = germline.collect{|m| m.start_with?('chr') ? m : 'chr' + m }

    if all_mutations.nil? || all_mutations.empty?
      all_mutations = somatic + germline
    else
      all_mutations = all_mutations.keys
      all_mutations += germline + somatic
    end

    all_mutations = all_mutations.collect{|m| m.sub(/^chr/, '') }.uniq

    if false && svs.empty?
      mutation_translations = TSV.setup({}, "Genomic Mutation~Translation", :type => :flat)
      all_mutations.each do |mutation|
        mutation_translations[mutation] = [mutation]
      end
      reference_io = Open.open(reference)
    else
      reference_io, mutation_translations = HTSBenchmark.apply_SVs(reference, all_mutations, svs)
    end

    Open.write(file('mutation_translations.tsv'), mutation_translations.to_s)

    Open.mkdir files_dir 
    target_reference = file('hg38.fa')
    Open.write(target_reference, reference_io)
    #CMD.cmd('bgzip', " -c > #{target_reference}", :in => reference_io)

    Open.write(file('somatic'), somatic.collect{|m| m = m.sub(/^chr/,''); mutation_translations[m].collect{|m| m.sub(/^chr/, '') } }.flatten.compact * "\n")
    Open.write(file('germline'), germline.collect{|m| m = m.sub(/^chr/,''); mutation_translations[m].collect{|m| m.sub(/^chr/, '') } }.flatten.compact * "\n")
    Open.write(file('all_mutations'), all_mutations.collect{|m| m = m.sub(/^chr/,''); mutation_translations[m].collect{|m| m.sub(/^chr/, '') } }.flatten.compact * "\n")

    CMD.cmd(:samtools, "faidx #{target_reference}")

    ['somatic', 'germline', 'all_mutations'].each do |type|
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

    Open.link(Dir.glob(files_dir + "/*.fa.gz").first, self.tmp_path)
    nil
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
    TSV.traverse Dir.glob(dependencies.first.files_dir + "/*.fa.gz").first, :type => :array do |line|
      if m = line.match(/^>([^\s]*)/)
        chr = m.captures[0]
      else
        sizes[chr] = line.length
      end
    end
    sizes
  end

  input :only_tumor, :boolean, "Only simulate the tumor", false
  dep :simulate_germline_hg38_vcf, :jobname => "Default"
  dep :simulate_somatic_hg38_vcf, :jobname => "Default"
  dep :miniref, :vcf => :simulate_somatic_hg38_vcf, :jobname => "hg38"
  dep :miniref_ploidy, :vcf => :simulate_somatic_hg38_vcf, :jobname => "hg38"
  dep :minify_vcf, :vcf_file => :simulate_germline_hg38_vcf, :sizes => :miniref_sizes, :jobname => "germline"
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
  dep :SV_miniref, :compute => :produce
  dep Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname,options,dependencies|
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.last
    file = sv_miniref.file('germline.reference')
    {:inputs => options.merge("Sequence#reference" => file), :jobname => jobname + '-germline' }
  end
  dep Sequence, :mutations_to_vcf, :positions => :placeholder, :organism => HG38_ORGANISM do |jobname,options,dependencies|
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.last
    file = sv_miniref.file('somatic.reference')
    {:inputs => options.merge("Sequence#reference" => file), :jobname => jobname + '-somatic' }
  end
  dep :SV_miniref_sizes
  dep :minify_vcf, :vcf_file => :placeholder, :sizes => :SV_miniref_sizes, :jobname => "germline" do |jobname,options,dependencies|
    vcf_file = dependencies.flatten.select{|dep| dep.task_name.to_s == 'mutations_to_vcf' }[-2]
    {:inputs => options.merge(:vcf_file => vcf_file) }
  end
  dep :minify_vcf, :vcf_file => :placeholder, :sizes => :SV_miniref_sizes, :jobname => "somatic" do |jobname,options,dependencies|
    vcf_file = dependencies.flatten.select{|dep| dep.task_name.to_s == 'mutations_to_vcf' }[1]
    {:inputs => options.merge(:vcf_file => vcf_file) }
  end
  dep HTSBenchmark, :NEAT_genreads, :somatic => :placeholder, :germline => :placeholder, :reference => :placeholder, :produce_samples => :tumor do |jobname,options,dependencies|
    orig_germline, germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}
    sv_miniref = dependencies.flatten.select{|dep| dep.task_name.to_s == 'SV_miniref' }.first
    reference = Path.setup(sv_miniref.files_dir + "/hg38.fa.gz")
    options[:reference] = reference
    options[:germline] = germline
    options[:somatic] = somatic
    {:inputs => options, :jobname => jobname + '-tumor'}
  end
  task :simulated_SV_sample => :text do |only_tumor|

    if only_tumor
      tumor = dependencies.select{|d| d.task_name.to_s == 'NEAT_genreads' }.first
      Open.link(tumor.file('output/tumor_read1.fq.gz'), file('output/tumor_read1.fq.gz'))
      Open.link(tumor.file('output/tumor_read2.fq.gz'), file('output/tumor_read2.fq.gz'))
    else
      normal, tumor = dependencies.select{|d| d.task_name.to_s == 'NEAT_genreads' }
      Open.link(normal.file('output/normal_read1.fq.gz'), file('output/normal_read1.fq.gz'))
      Open.link(normal.file('output/normal_read2.fq.gz'), file('output/normal_read2.fq.gz'))
      Open.link(tumor.file('output/tumor_read1.fq.gz'), file('output/tumor_read1.fq.gz'))
      Open.link(tumor.file('output/tumor_read2.fq.gz'), file('output/tumor_read2.fq.gz'))
    end

    Dir.glob(files_dir + "/output/*.fq.gz") * "\n"
  end
end

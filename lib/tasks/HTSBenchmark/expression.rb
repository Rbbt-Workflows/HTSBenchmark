module HTSBenchmark

  dep :genotype_germline_hg38
  dep_task :germline_transcript_offsets, Sequence, :transcript_offsets, :compute => :produce, :positions => :genotype_germline_hg38, :organism => HTSBenchmark.organism, :vcf => false

  dep :genotype_germline_hg38
  dep :population_genotypes, :germline => :genotype_germline_hg38
  dep_task :somatic_transcript_offsets, Sequence, :transcript_offsets, :compute => :produce, :positions => :placeholder, :organism => HTSBenchmark.organism, :vcf => false do |jobname,options,dependencies|
    mutations = dependencies.flatten.last.file('total_mutations_clean')
    options[:positions] = mutations
    {:inputs => options}
  end

  dep :genotype_germline_hg38
  dep :population_genotypes, :germline => :genotype_germline_hg38
  dep_task :SV_partial_genes, Sequence, :genes_at_ranges, :compute => :produce, :ranges => :placeholder, :organism => HTSBenchmark.organism, :full_overlap => false do |jobname,options,dependencies|
    ranges = dependencies.flatten.last.file('SV_ranges')
    options[:ranges] = ranges
    {:inputs => options}
  end

  dep :genotype_germline_hg38
  dep :population_genotypes, :germline => :genotype_germline_hg38
  dep_task :SV_full_genes, Sequence, :genes_at_ranges, :compute => :produce, :ranges => :placeholder, :organism => HTSBenchmark.organism, :full_overlap => true do |jobname,options,dependencies|
    ranges = dependencies.flatten.last.file('SV_ranges')
    options[:ranges] = ranges
    {:inputs => options}
  end

  dep :germline_transcript_offsets
  dep :somatic_transcript_offsets
  dep :SV_partial_genes
  dep :SV_full_genes
  task :expression_consequence => :tsv do
    Step.wait_for_jobs dependencies
    transcript_offsets_job_germline, transcript_offsets_job, genes_partial_job, genes_full_job = dependencies
    pop_job = transcript_offsets_job.step(:population_genotypes)

    transcript_offsets_germline = transcript_offsets_job_germline.load
    transcript_offsets = transcript_offsets_job.load

    somatic_mutations = transcript_offsets.keys
    germline_mutations = transcript_offsets_germline.keys

    genes_partial = genes_partial_job.load
    genes_full = genes_full_job.load

    clone_dirs = Dir.glob(File.join(pop_job.files_dir, 'clone_*'))

    transcript_mutations = []
    gene_svs = []
    clone_dirs.each do |clone_dir|
      clone = File.basename(clone_dir)
      Path.setup(clone_dir)
      svs = clone_dir.all_SVs.tsv
      mutations = clone_dir.mutations.list # ToDo change to all mutations

      mutations.each do |hap_mut|
        copy, _sep, mut = hap_mut.partition("_")
        chr, pos, change = mut.split(":")
        mut.sub!(/^chr/,'')
        (transcript_offsets[mut] || []).each do |toff|
          transcript, offset, strand = toff.split(":")
          transcript_mutations << [clone, transcript, copy, offset, strand, change, mut]
        end
      end

      svs.each do |id, values|
        type, chr_str, start, eend = values

        copy, chr = chr_str.split("_")

        range = [chr, start, eend] * ":"

        full_genes = genes_full[range]
        partial_genes = genes_partial[range] - full_genes

        case type
        when "DEL"
          partial_genes.each do |gene|
            gene_svs << [clone, gene, copy, "DEL"]
          end
        when "INS", "INV"
          full_genes.each do |gene|
            gene_svs << [clone, gene, copy, "DUP"]
          end
        end
      end
    end

    germline_mutations.each do |mut|
      hap_mut = HTSBenchmark.haploid_mutation(mut)
      copy, _sep, mut = hap_mut.partition("_")
      chr, pos, change = mut.split(":")
      mut.sub!(/^chr/,'')
      (transcript_offsets_germline[mut] || []).each do |toff|
        transcript, offset, strand = toff.split(":")
        transcript_mutations << ['germline', transcript, copy, offset, strand, change, mut]
      end
    end


    dumper = TSV::Dumper.new :key_field => "Entity", :fields => ["Copy", "Clone", "Change", "Location"], :type => :list
    dumper.init
    t = Thread.new do 
      transcript_mutations.each do |clone, transcript, copy, offset, strand, change, mut|
        dumper.add transcript, [copy, clone, mut, [offset, strand] * ":"]
      end
      gene_svs.each do |clone, gene, copy, type|
        dumper.add gene, [copy, clone, type, nil]
      end
      dumper.close
    end

    stream = dumper.stream
    ConcurrentStream.setup stream, :threads => t
  end

end

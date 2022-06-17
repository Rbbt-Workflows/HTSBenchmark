module HTSBenchmark

  dep :population_genotypes
  dep Sequence, :transcript_offsets, :compute => :produce, :positions => :placeholder, :organism => HTSBenchmark.organism do |jobname,options,dependencies|
    mutations = dependencies.flatten.first.file('total_mutations_clean')
    options[:positions] = mutations
    {:inputs => options}
  end
  dep Sequence, :genes_at_ranges, :compute => :produce, :ranges => :placeholder, :organism => HTSBenchmark.organism do |jobname,options,dependencies|
    ranges = dependencies.flatten.first.file('SV_ranges')
    options[:ranges] = ranges
    {:inputs => options}
  end
  dep Sequence, :genes_at_ranges, :compute => :produce, :ranges => :placeholder, :full_overlap => true, :organism => HTSBenchmark.organism do |jobname,options,dependencies|
    ranges = dependencies.flatten.first.file('SV_ranges')
    options[:ranges] = ranges
    {:inputs => options}
  end
  task :expression_consequence => :text do
    Step.wait_for_jobs dependencies
    pop_job, transcript_offsets_job, genes_partial_job, genes_full_job = dependencies

    transcript_offsets = transcript_offsets_job.load
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

    sout = Misc.open_pipe do |sin|
      transcript_mutations.each do |clone, transcript, copy, offset, strand, change, mut|
        sin.puts [transcript, copy, clone, mut, [offset, strand] * ":"] * "\t"
      end
      gene_svs.each do |clone, gene, copy, type|
        sin.puts [gene, copy, clone, type] * "\t"
      end
    end
    sout
  end

end

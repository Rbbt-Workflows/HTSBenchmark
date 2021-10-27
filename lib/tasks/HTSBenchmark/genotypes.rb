module HTSBenchmark
  
  HG19_ORGANISM="Hsa/feb2014"
  HG38_ORGANISM="Hsa/may2017"

  #{{{ GERMLINE

  task :genotype_germline_hg19_all_chr => :array do 
    Workflow.require_workflow "Genomes1000"
    TSV.traverse Genomes1000.rsids, :fields => ["Genomic Mutation", "EUR_AF"], :type => :list, :into => :stream, :bar => self.progress_bar("Selecting germline mutations") do |rsid,values|
      mutation, af = values
      next unless rand < af.to_f
      mutation
    end
  end

  input :chromosome, :string, "Chromosome to choose from", nil
  dep :genotype_germline_hg19_all_chr
  task :genotype_germline_hg19 => :array do |chr|
    chr = chr.sub('chr', '') if chr
    TSV.traverse step(:genotype_germline_hg19_all_chr), :type => :array, :into => :stream do |mutation|
      next if chr && mutation.split(":").first.sub('chr','') != chr
      mutation
    end
  end

  dep :genotype_germline_hg19
  dep_task :genotype_germline_hg38_lf, Sequence, :lift_over, :positions => :genotype_germline_hg19, :source => HG19_ORGANISM, :target => HG38_ORGANISM

  dep :genotype_germline_hg38_lf
  task :genotype_germline_hg38_lf_chr => :array do 
    chr = recursive_inputs[:chromosome]
    TSV.traverse step(:genotype_germline_hg38_lf), :type => :array, :into => :stream do |line|
      next if ! line.include?(":") || line.split(":").first.include?("_")
      next if chr && line.split(":").first != chr
      line
    end
  end

  dep :genotype_germline_hg38_lf_chr
  dep Sequence, :reference, :positions => :genotype_germline_hg38_lf_chr, :organism => HG38_ORGANISM, :vcf => false, :full_reference_sequence => false
  task :genotype_germline_hg38 => :array do 
    TSV.traverse step(:reference), :into => :stream do |mutation, reference|
      next if mutation.split(":")[2].split(",").include? reference
      mutation
    end
  end

  #{{{ SOMATIC
   
  input :study, :string, "PCAWG Study code", "Bladder-TCC"
  input :chromosome, :string, "Chromosome to choose from", nil
  input :max, :integer, "Maximum number of mutations", 20_000_000
  task :genotype_somatic_hg19_PCAWG => :array do |study,chr,max|
    chr = chr.sub('chr', '') if chr
    Workflow.require_workflow "PCAWG"
    mutations = Study.setup(study).genomic_mutations
    mutations = mutations.select{|mutation| mutation.split(":").first.sub('chr','') == chr}  if chr
    mutations = mutations.reject{|mutation| mutation.split(":").first.include? "_"}
    mutations = mutations.shuffle[0..max-1] if mutations.length > max
    mutations
  end

  dep :genotype_somatic_hg19_PCAWG
  dep_task :genotype_somatic_hg38_lf, Sequence, :lift_over, :positions => :genotype_somatic_hg19_PCAWG, :source => HG19_ORGANISM, :target => HG38_ORGANISM

  dep :genotype_somatic_hg38_lf
  task :genotype_somatic_hg38_lf_chr => :array do 
    chr = recursive_inputs[:chromosome]
    TSV.traverse step(:genotype_somatic_hg38_lf), :type => :array, :into => :stream do |line|
      next if ! line.include?(":") || line.split(":").first.include?("_")
      next if chr && line.split(":").first != chr
      line
    end
  end

  dep :genotype_somatic_hg38_lf_chr
  dep Sequence, :reference, :positions => :genotype_somatic_hg38_lf_chr, :organism => 'Hsa/may2017', :vcf => false, :full_reference_sequence => false
  task :genotype_somatic_hg38 => :array do 
    TSV.traverse step(:reference), :into => :stream do |mutation, reference|
      next if mutation.split(":")[2].split(",").include? reference
      mutation
    end
  end

  helper :target_loc do |schr,bps,types,total|
    trans_freq = types["TRA"].to_f / total
    trans = rand < trans_freq

    tchr = nil
    if trans
      pos = rand(total)
      current = 0
      bps.each do |chr,l|
        current += l.length
        tchr = chr
        break if current > pos
      end
    else
      tchr = schr
    end

    tpos = bps[tchr][rand(bps[tchr].length)]

    [tchr, tpos]
  end

  input :study, :string, "PCAWG Study code", "Bladder-TCC"
  input :chromosome, :string, "Chromosome to choose from", nil
  input :max, :integer, "Maximum number of mutations", 20_000_000
  task :SV_somatic_hg19_PCAWG => :tsv do |study,chr,max|
    chr = chr.sub('chr', '') if chr
    Workflow.require_workflow "PCAWG"

    bps = {}
    types = {}
    Study.setup(study).donors.each do |donor|
      donor.SV.each do |sv|
        chr, pos, type, chr2, pos2 = sv.split(":")
        bps[chr] ||= []
        bps[chr2] ||= []
        bps[chr] << pos.to_i
        bps[chr2] << pos2.to_i
        bps[chr2] << pos2.to_i + 1
        types[type] ||= 0
        types[type] += 1
      end
      donor.cnvs.each do |cnv|
        chr, pos, pos2 = cnv.split(":")
        bps[chr] ||= []
        bps[chr] << pos.to_i
        bps[chr] << pos2.to_i
      end
    end

    total = Misc.sum(types.values)

    svs = TSV.setup([], :key_field => "SV ID", :fields => ["Type", "Chromosome", "Start", "End", "Target chromosome", "Target start", "Target end"], :type => :list)
    Study.setup(study).donors.each do |donor|
      donor.SV.each do |sv|
        chr, pos, type, chr2, pos2 = sv.split(":")
        pos = pos.to_i

        tchr, tpos = target_loc chr, bps, types, total
        bps[tchr] ||= []
        bps[tchr] << tpos.to_i

        bps[chr] ||= []
        bps[chr] << pos.to_i
        bps[chr] << pos2.to_i


        id = Misc.digest([sv, tchr, tpos] * ":")
        tpos = tpos.to_i

        prob_conseq_dup = 0.7
        case type
        when "DUP"
          if rand < prob_conseq_dup
            svs[id] = ["INS", chr, pos, pos2, tchr, tpos]
          else
            svs[id] = ["INS", chr, pos, pos2, chr, pos]
          end
        when "DEL"
          svs[id] = ["DEL", chr, pos, pos2]
        when "h2hINV"
          svs[id] = ["INV", chr, pos, pos2, chr, pos]
        when "t2tINV"
          svs[id] = ["INV", chr, pos, pos2, chr, pos2.to_i + 1]
        when "TRA"
          closest = bps[chr].sort.select{|p| p > pos }.first
          svs[id] = ["INS", chr, pos, closest, chr2, pos2]
          svs[id] = ["DEL", chr, pos, closest]
        end
      end
    end
    
    positions = []
    bps.each do |chr,list|
      positions += list.sort.collect{|p| [chr, p] * ":" } 
    end

    Open.write(file('positions'), positions.uniq * "\n")

    svs
  end


  dep :SV_somatic_hg19_PCAWG
  dep Sequence, :lift_over, :positions => :placeholder, :source => HG19_ORGANISM, :target => HG38_ORGANISM do |jobname, options,dependencies|
    positions = dependencies.flatten.first.file('positions')
    {:inputs => options.merge(:positions => positions)}
  end
  task :SV_somatic_hg38_PCAWG => :tsv do 
    positions = step(:SV_somatic_hg19_PCAWG).file('positions').list
    lifted = step(:lift_over).load
    translations = Misc.zip2hash(positions, lifted)

    parser = TSV::Parser.new step(:SV_somatic_hg19_PCAWG)
    dumper = TSV::Dumper.new parser.options
    dumper.init
    TSV.traverse parser, :into => dumper do  |id, values|
      type, chr, start, eend, target_chr, target_start, target_end = values

      ps = [chr, start] * ":"
      nps = translations[ps]
      next if nps.nil? || nps.empty?
      chr, start = nps.split(":")

      pe = [chr, eend] * ":"
      npe = translations[pe]
      next if npe.nil? || npe.empty?
      chr, eend = npe.split(":")

      if target_chr
        tps = [target_chr, target_start] * ":"
        tnps = translations[tps]
        next if tnps.nil? || tnps.empty?
        target_chr, target_start = tnps.split(":")

        if target_end
          tpe = [target_chr, target_end] * ":"
          tnpe = translations[tpe]
          next if tnpe.nil? || tnpe.empty?
          target_chr, target_eend = tnpe.split(":")
        end
      end
      
      values = type, chr, start, eend, target_chr, target_start, target_end 
      [id, values]
    end
  end

end

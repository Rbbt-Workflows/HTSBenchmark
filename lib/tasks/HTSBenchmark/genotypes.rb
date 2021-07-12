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

end

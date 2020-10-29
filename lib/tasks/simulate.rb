module HTSBenchmark

  #{{{ GERMLINE

  input :chromosome, :string, "Chromosome to choose from", nil
  task :simulate_germline_hg19 => :array do |chr|
    chr = chr.sub('chr', '') if chr
    Workflow.require_workflow "Genomes1000"
    TSV.traverse Genomes1000.rsids, :fields => ["Genomic Mutation", "EUR_AF"], :type => :list, :into => :stream, :bar => self.progress_bar("Selecting germline mutations") do |rsid,values|
      mutation, af = values
      next if chr && mutation.split(":").first.sub('chr','') != chr
      next unless rand < af.to_f
      mutation
    end
  end

  dep :simulate_germline_hg19, :compute => :produce
  dep_task :simulate_germline_hg38, Sequence, :lift_over, :positions => :simulate_germline_hg19, :source => "Hsa/feb2014", :target => "Hsa/dec2016"

  dep :simulate_germline_hg38
  task :simulate_germline_hg38_clean => :array do 
    TSV.traverse step(:simulate_germline_hg38), :type => :array, :into => :stream do |line|
      next if ! line.include?(":") || line.split(":").first.include?("_")
      line
    end
  end

  dep :simulate_germline_hg19
  extension :vcf
  dep_task :simulate_germline_hg19_vcf, Sequence, :mutations_to_vcf, :positions => :simulate_germline_hg19, :organism => "Hsa/may2017"

  dep :simulate_germline_hg38_clean
  extension :vcf
  dep_task :simulate_germline_hg38_vcf, Sequence, :mutations_to_vcf, :positions => :simulate_germline_hg38_clean, :organism => "Hsa/may2017"

  #{{{ SOMATIC
   
  input :study, :string, "PCAWG Study code", "Bladder-TCC"
  input :chromosome, :string, "Chromosome to choose from", nil
  input :max, :integer, "Maximum number of mutations", 20_000_000
  task :simulate_somatic_hg19_PCAWG => :array do |study,chr,max|
    chr = chr.sub('chr', '') if chr
    Workflow.require_workflow "PCAWG"
    mutations = Study.setup(study).genomic_mutations
    mutations = mutations.select{|mutation| mutation.split(":").first.sub('chr','') == chr}  if chr
    mutations = mutations.reject{|mutation| mutation.split(":").first.include? "_"}
    mutations = mutations.shuffle[0..max-1] if mutations.length > max
    mutations
  end

  dep :simulate_somatic_hg19_PCAWG, :compute => :produce
  dep_task :simulate_somatic_hg38, Sequence, :lift_over, :positions => :simulate_somatic_hg19_PCAWG, :source => "Hsa/feb2014", :target => "Hsa/dec2016"

  dep :simulate_somatic_hg38
  task :simulate_somatic_hg38_clean => :array do 
    TSV.traverse step(:simulate_somatic_hg38), :type => :array, :into => :stream do |line|
      next if !line.include?(":") || line.split(":").first.include?("_")
      line
    end
  end


  dep :simulate_somatic_hg19_PCAWG
  dep_task :simulate_somatic_hg19_vcf, Sequence, :mutations_to_vcf, :positions => :simulate_somatic_hg19_PCAWG, :organism => "Hsa/may2017"

  dep :simulate_somatic_hg38_clean
  extension :vcf
  dep_task :simulate_somatic_hg38_vcf, Sequence, :mutations_to_vcf, :positions => :simulate_somatic_hg38_clean, :organism => "Hsa/may2017"


  input :somatic, :file, "Somatic VCF", nil, :nofile => true
  input :germline, :file, "Germline VCF", nil, :nofile => true
  input :reference, :select, "Reference code", "b37", :select_options => %w(b37 hg38 mm10), :nofile => true
  input :normal_depth, :integer, "Depth on normal sample", 30
  input :tumor_depth, :integer, "Depth on normal sample", 80
  task :NEAT_genreads => :text do |somatic,germline,reference,normal_depth,tumor_depth|
    input = file('input')
    Open.mkdir input

    Open.cp(somatic, input.somatic)
    Open.cp(germline, input.germline)

    CMD.cmd("cat #{input.somatic} > #{input.joined}")
    CMD.cmd("cat #{input.germline} |grep -v '#' >> #{input.joined}")

    reference_gunzip = input[File.basename(reference).sub('.gz','')]
    CMD.cmd(:zcat, "'#{reference}' > #{reference_gunzip}")

    output = file('output')
    Open.mkdir output

    if reference.include?("hg38")
      [input.somatic, input.germline, input.joined].each do |file|
        TmpFile.with_file do |tmpfile|
          Open.open(tmpfile, :mode => "w") do |sin|
            TSV.traverse file, :type => :array do |line|

              l = if line =~ /^(?:##)/ 
                    line
                  elsif line =~ /^#CHR/
                    line + "\t" + "FORMAT" + "\t" + "Sample"
                  else
                    line = "chr" + line unless line =~ /^chr/
                    parts = line.split("\t")[0..4]
                    parts[4] = parts[4].split(",").first if parts[4]
                    (parts + [".", "PASS", ".", "GT", "0|1"]) * "\t"
                  end

              sin.puts l
            end
          end
          Open.mv tmpfile, file
        end
      end
    end

    CMD.cmd_log("genReads.py", "-c #{normal_depth} -r '#{reference_gunzip}' -M 0 -R 100 --pe 100 10 -o '#{output.normal}' -v '#{input.germline}' --bam --vcf")
    CMD.cmd_log("genReads.py", "-c #{tumor_depth} -r '#{reference_gunzip}' -M 0 -R 100 --pe 100 10 -o '#{output.tumor}' -v '#{input.joined}' --bam --vcf")

    output.glob("*.fq").each do |file|
      CMD.cmd("bgzip #{file}")
    end
    
    "DONE"
  end
end

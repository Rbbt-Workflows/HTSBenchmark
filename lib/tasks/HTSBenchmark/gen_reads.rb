require 'tools/samtools'

module HTSBenchmark

  dep :mutations_to_reference
  input :reference, :binary, "Reference file", nil, :nofile => true
  input :depth, :integer, "Sequencing depth to simulate", 60
  input :haploid_reference, :boolean, "Reference is haploid (each chromosome copy separate)"
  input :sample_name, :string, "Sample name", nil, :jobname => true
  dep Sequence, :mutations_to_vcf, "Sequence#reference" => :mutations_to_reference, :mutations => :skip
  task :NEAT_simulate_DNA_single => :array do |reference,depth,haploid,sample_name|

    if haploid
      depth = (depth.to_f / 2).ceil
    end

    mutations_vcf = file('mutations.vcf')
    Open.write(mutations_vcf) do |sin|
      vcf = step(:mutations_to_vcf)
      #vcf = Sequence.job(:mutations_to_vcf, nil, "Sequence#reference" => step(:mutations_to_reference))
      #Misc.with_env "RBBT_UPDATE", "TRUE" do
      #  vcf.clean unless vcf.updated?
      #end
      TSV.traverse vcf, :type => :array do |line|
        l = if line =~ /^(?:##)/ 
            line
            elsif line =~ /^#CHR/
              line + "\t" + "FORMAT" + "\t" + "Sample"
            else
              line = "chr" + line unless line =~ /^chr/ || line =~ /^copy/
              parts = line.split("\t")[0..4]

              parts[4] = parts[4].split(",").first if parts[4]

              (parts + [".", "PASS", ".", "GT", (haploid ? "1|1" : (rand < 0.5 ? "0|1" : "1|0"))]) * "\t"
            end
        sin.puts l
      end
    end

    output = file('output')
    reference_gunzip = file('hg38.fa')
    CMD.cmd(:zcat, "'#{reference}' > #{reference_gunzip}")

    Open.mkdir output
    CMD.cmd_log("gen_reads.py", "-c #{depth} -r '#{reference_gunzip}' -p 2 -M 0 -R 101 --pe 100 10 -o '#{output[sample_name]}' -v '#{mutations_vcf}' --vcf --bam")

    output.glob("*.fq").each do |file|
      CMD.cmd("bgzip #{file}")
    end

    output.glob("*.fq.gz")
  end

  dep :mutations_to_reference
  input :reference, :binary, "Reference file", nil, :nofile => true
  input :depth, :integer, "Sequencing depth to simulate", 60
  input :haploid_reference, :boolean, "Reference is haploid (each chromosome copy separate)"
  input :sample_name, :string, "Sample name", nil, :jobname => true
  dep Sequence, :mutations_to_vcf, "Sequence#reference" => :mutations_to_reference, :mutations => :skip
  task :NEAT_simulate_DNA => :array do |reference,depth,haploid,sample_name|

    if haploid
      depth = (depth.to_f / 2).ceil
    end

    mutations_vcf = file('mutations.vcf')
    Open.write(mutations_vcf) do |sin|
      vcf = step(:mutations_to_vcf)
      #vcf = Sequence.job(:mutations_to_vcf, nil, "Sequence#reference" => step(:mutations_to_reference))
      #Misc.with_env "RBBT_UPDATE", "TRUE" do
      #  vcf.clean unless vcf.updated?
      #end
      TSV.traverse vcf, :type => :array do |line|
        l = if line =~ /^(?:##)/ 
            line
            elsif line =~ /^#CHR/
              line + "\t" + "FORMAT" + "\t" + "Sample"
            else
              line = "chr" + line unless line =~ /^chr/ || line =~ /^copy/
              parts = line.split("\t")[0..4]

              parts[4] = parts[4].split(",").first if parts[4]

              (parts + [".", "PASS", ".", "GT", (haploid ? "1|1" : (rand < 0.5 ? "0|1" : "1|0"))]) * "\t"
            end
        sin.puts l
      end
    end

    chr_output = file('chr_output')
    reference_gunzip = file('hg38.fa')
    io = CMD.cmd(:zcat, "'#{reference}'", :pipe => true)

    Open.mkdir chr_output

    chrs = []
    file = nil
    TSV.traverse io, :type => :array do |line|

      if m = line.match(/>(\w+)/)
        chr = m.captures[0]
        chrs << chr
        file.close if file
        file = Open.open(chr_output[chr].reference, :mode => 'w')
      end

      file.puts line
    end
    file.close


    cpus = config(:cpus, :genReads, :NEAT, :gen_reads)
    TSV.traverse chrs, :type => :array, :cpus => cpus, :bar => "Generating reads by chromosome" do |chr|
      Open.mkdir chr_output[chr]
      reference = chr_output[chr].reference
      CMD.cmd_log("gen_reads.py", "-c #{depth} -r '#{reference}' -p 2 -M 0 -R 101 --pe 100 10 -o '#{chr_output[chr][sample_name]}' -v '#{mutations_vcf}' --vcf --bam")
    end

    
    output = file('output')

    # Merge VCF
    vcf = Open.open(output[sample_name] + ".vcf", :mode => "w")

    header = true
    chr_output.glob("*/*.vcf.gz").each do |file|
      TSV.traverse file, :type => :array do |line|
        next if not header and line =~ /^#/
        next if line =~ /^##reference/
        vcf.puts line
      end
      header = false
    end
    vcf.close

    # Merge BAM
    bam = output[sample_name] + ".bam"
    bam_parts = chr_output.glob("*/*.bam")

    CMD.cmd(:samtools, "merge '#{bam}' #{bam_parts * " "}")

    # Merge FASTQ
    fq1 = output[sample_name] + "_read1.fq"
    chr_output.glob("*/*_read1.fq.gz").each do |file|
      CMD.cmd("zcat '#{file}' >> '#{fq1}'")
    end
    
    fq2 = output[sample_name] + "_read2.fq"
    chr_output.glob("*/*_read2.fq.gz").each do |file|
      CMD.cmd("zcat '#{file}' >> '#{fq2}'")
    end

    # Compress
    output.glob("*.fq").each do |file|
      CMD.cmd("bgzip '#{file}'")
    end

    CMD.cmd("bgzip #{output[sample_name] + ".vcf"}")

    # Cleanup parts
    FileUtils.rm_rf chr_output

    output.glob("*.fq.gz")
  end

end

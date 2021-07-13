require 'tools/samtools'

module HTSBenchmark


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

    output = file('output')
    reference_gunzip = file('hg38.fa')
    CMD.cmd(:zcat, "'#{reference}' > #{reference_gunzip}")

    Open.mkdir output
    CMD.cmd_log("gen_reads.py", "-c #{depth} -r '#{reference_gunzip}' -p 2 -M 0 -R 101 --pe 100 10 -o '#{output[sample_name]}' -v '#{mutations_vcf}' --vcf")

    output.glob("*.fq").each do |file|
      CMD.cmd("bgzip #{file}")
    end

    output.glob("*.fq.gz")
  end

end

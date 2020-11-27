module HTSBenchmark

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

    bam = output['tumor_golden.bam']
    bam_normal = output['normal_golden.bam']
    #CMD.cmd("samtools view -H #{bam} | sed \"s/SM:[^\\t]*/SM:NEAT_tumor/g\" | samtools reheader - #{bam} > #{bam}.tmp; mv #{bam}.tmp #{bam}")
    #CMD.cmd("samtools view -H #{bam_normal} | sed \"s|@SQ.*|@SQ\tSN:chr3\tLN:467844\tM5:6a81612c18fa6e87418e6229be2293ac\tUR:file:/home/mvazque2/.rbbt/var/fasta_indices/90c12e17bb9b921fd43605d53a28397a/hg38.fa.gz|\" | samtools reheader - #{bam_normal} > #{bam_normal}.tmp; mv #{bam_normal}.tmp #{bam_normal}")
    #
    CMD.cmd("samtools view -H #{bam_normal} | sed \"s/:NEAT/:NEAT_normal/g\" | samtools reheader - #{bam_normal} > #{bam_normal}.tmp; mv #{bam_normal}.tmp #{bam_normal}")

    CMD.cmd("samtools addreplacerg  -o #{bam_normal}.tmp -R 'NEAT_normal' #{bam_normal}  ; samtools view -b -S #{bam_normal}.tmp > #{bam_normal}; rm #{bam_normal}.tmp")
    CMD.cmd("samtools addreplacerg  -o #{bam}.tmp -R 'NEAT' #{bam}  ; samtools view -b -S #{bam}.tmp > #{bam}; rm #{bam}.tmp")
    
    "DONE"
  end
end

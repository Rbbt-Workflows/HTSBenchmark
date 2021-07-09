module HTSBenchmark


  input :somatic, :file, "Somatic VCF", nil, :nofile => true
  input :germline, :file, "Germline VCF", nil, :nofile => true
  input :reference, :file, "Reference file", nil, :nofile => true
  input :normal_depth, :integer, "Depth on normal sample", 30
  input :tumor_depth, :integer, "Depth on normal sample", 80
  input :haploid, :boolean, "Simulate with ploidy 1", false
  input :produce_samples, :select, "What samples to produce", :both, :select_options => %w(both tumor normal)
  task :NEAT_genreads => :text do |somatic,germline,reference,normal_depth,tumor_depth,haploid,produce_samples|
    input = file('input')
    Open.mkdir input

    Open.cp(somatic, input.somatic) if somatic
    Open.cp(germline, input.germline)

    CMD.cmd("cat #{input.germline} > #{input.joined}", :nofail => true)
    CMD.cmd("cat #{input.somatic} |grep -v '#' >> #{input.joined}", :nofail => true) if somatic

    #reference_gunzip = input[File.basename(reference).sub('.gz','')]
    reference_gunzip = input['hg38.fa']
    CMD.cmd(:zcat, "'#{reference}' > #{reference_gunzip}")

    output = file('output')
    Open.mkdir output

    if reference.include?("hg38")
      [input.somatic, input.germline, input.joined].each do |file|
        next unless File.exists?(file)
        TmpFile.with_file do |tmpfile|
          Open.open(tmpfile, :mode => "w") do |sin|
            TSV.traverse file, :type => :array do |line|

              l = if line =~ /^(?:##)/ 
                    line
                  elsif line =~ /^#CHR/
                    line + "\t" + "FORMAT" + "\t" + "Sample"
                  else
                    line = "chr" + line unless line =~ /^chr/ || line =~ /^copy/
                    parts = line.split("\t")[0..4]
                    parts[4] = parts[4].split(",").first if parts[4]
                    (parts + [".", "PASS", ".", "GT", haploid ? "1|1" : "0|1"]) * "\t"
                  end

              sin.puts l
            end
          end
          Open.mv tmpfile, file
        end
      end
    end

    if haploid
      normal_depth = (normal_depth.to_f / 2).ceil
      tumor_depth = (tumor_depth.to_f / 2).ceil
    end

    case produce_samples.to_s
    when 'both'
      produce_tumor = true
      produce_normal = true
    when 'normal'
      produce_normal = true
    when 'tumor'
      produce_tumor = true
    end

    cmd_tool = "gen_reads.py"
    if produce_tumor
      CMD.cmd_log(cmd_tool, "-c #{tumor_depth} -r '#{reference_gunzip}' -p 2 -M 0 -R 101 --pe 100 10 -o '#{output.tumor}' -v '#{input.joined}' --bam --vcf")
      bam = output['tumor_golden.bam']
      CMD.cmd("samtools addreplacerg  -o #{bam}.tmp -R 'NEAT' #{bam}  ; samtools view -b -S #{bam}.tmp > #{bam}; rm #{bam}.tmp")
    end

    if produce_normal
      CMD.cmd_log(cmd_tool, "-c #{normal_depth} -r '#{reference_gunzip}' -p 2 -M 0 -R 101 --pe 100 10 -o '#{output.normal}' -v '#{input.germline}' --bam --vcf")
      bam_normal = output['normal_golden.bam']
      CMD.cmd("samtools view -H #{bam_normal} | sed \"s/:NEAT/:NEAT_normal/g\" | samtools reheader - #{bam_normal} > #{bam_normal}.tmp; mv #{bam_normal}.tmp #{bam_normal}")
      CMD.cmd("samtools addreplacerg  -o #{bam_normal}.tmp -R 'NEAT_normal' #{bam_normal}  ; samtools view -b -S #{bam_normal}.tmp > #{bam_normal}; rm #{bam_normal}.tmp")
    end
    
    output.glob("*.fq").each do |file|
      CMD.cmd("bgzip #{file}")
    end

    "DONE"
  end

  dep :simulate_germline_hg38_vcf
  dep :simulate_somatic_hg38_vcf
  dep :miniref, :vcf => :simulate_somatic_hg38_vcf, :jobname => 'hg38'
  dep :minify_vcf, :vcf_file => :simulate_germline_hg38_vcf, :sizes => :miniref_sizes, :jobname => "germline"
  dep :minify_vcf, :vcf_file => :simulate_somatic_hg38_vcf, :sizes => :miniref_sizes, :jobname => "somatic"
  dep_task :simulated_sample, HTSBenchmark, :NEAT_genreads, :somatic => :placeholder, :germline => :placeholder, :reference => :miniref, :produce_samples => :both do |jobname,options,dependencies|
    germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}
    options[:germline] = germline
    options[:somatic] = somatic
    {:inputs => options, :jobname => jobname}
  end

  dep :simulate_germline_hg38_vcf
  dep :simulate_somatic_hg38_vcf
  dep :miniref_ploidy, :vcf => :simulate_somatic_hg38_vcf, :jobname => 'hg38'
  dep :minify_vcf, :vcf_file => :simulate_germline_hg38_vcf, :sizes => :miniref_sizes, :jobname => "germline"
  dep :minify_vcf, :vcf_file => :simulate_somatic_hg38_vcf, :sizes => :miniref_sizes, :jobname => "somatic"
  dep :vcf_ploidy, :haploid_reference => :miniref_ploidy do |jobname,options,dependencies|
    germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'minify_vcf'}
    [germline, somatic].collect do |step|
      {:inputs => options.merge(:vcf => step.path), :jobname => step.clean_name}
    end
  end
  dep_task :simulated_sample_ploidy, HTSBenchmark, :NEAT_genreads, :haploid => true, :somatic => :placeholder, :germline => :placeholder, :reference => :miniref_ploidy do |jobname,options,dependencies|
    germline, somatic = dependencies.flatten.select{|dep| dep.task_name.to_s == 'vcf_ploidy'}
    options[:germline] = germline
    options[:somatic] = somatic
    {:inputs => options, :jobname => jobname}
  end
end

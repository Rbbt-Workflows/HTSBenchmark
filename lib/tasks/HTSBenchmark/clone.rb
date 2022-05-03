module HTSBenchmark

  input :somatic, :array
  input :germline, :array
  dep :SV_mutations, :mutations => :somatic, :duplicate => true
  dep :SV_mutations, :mutations => :germline, :duplicate => true
  dep :merge_somatic_germline, :somatic => :placeholder, :germline => :placeholder do |jobname,options,dependencies|
    somatic, germline = dependencies.flatten.select{|d| d.task_name == :SV_mutations }
    {:inputs => options.merge(:somatic => somatic, :germline => germline), :jobname => jobname}
  end
  dep :SV_reference
  dep_task :clone, HTSBenchmark, :NEAT_simulate_DNA, :mutations => :merge_somatic_germline, :reference => :SV_reference, :haploid_reference => true


  input :clones, :array, "Array of NEAT job paths", nil, :nofile => true
  input :fractions, :array, "Array of clone cellular fraction", nil, :nofile => true
  task :merge_clones => :array do |clones, fractions|

    Open.mkdir files_dir
    ['tumor_read1.fq.gz', 'tumor_read2.fq.gz'].each_with_index do |filename,i|
      output = file(filename)

      sout = Misc.open_pipe false, false do |sin|

        clones.zip(fractions).each_with_index do |v,ci|
          clone, fraction = v
          fraction = fraction.to_f
          clone_step = Step === clone ? clone : Step.new(clone)

          input = clone_step.load[i]

          skip = nil
          rnd = Random.new 1234
          TSV.traverse input, :type => :array, :bar => "Processing #{ Misc.fingerprint [input, filename] }" do |line|
            if line =~ /^@.*clone_/
              rand = rnd.rand
              skip = rand > fraction 
              next if skip
            else
              next if skip
            end

            sin.puts line
          end
        end

        sin.close

      end

      CMD.cmd(:bgzip, "-c > #{output}", :in => sout)
    end
    Dir.glob(files_dir + "/*")
  end

end

module HTSBenchmark

  def self.annotate_BAM_with_readnames(bam, new, bar = nil, cpus = 3)
    sam = CMD.cmd(:samtools, "view -h #{bam}", :pipe => true, "--threads" => cpus)

    fixed = TSV.traverse sam, :type => :array, :into => :stream, :bar => bar  do |line|
      if line[0] == "@"
        line
      else
        name, *rest = line.split("\t")
        if m = name.match(/clone_(\d+)/)
          clone = m.captures.first
        end
        if m = name.match(/(copy-\d+_chr[0-9A-B])/i)
          chr_copy = m.captures.first
        end
        
        clone_copy = [clone, chr_copy] * "_" if clone && chr_copy
        rest << "XA:i:#{clone}" if clone
        rest << "XC:Z:#{chr_copy}" if chr_copy
        rest << "XP:Z:#{clone_copy}" if clone_copy
        name + "\t" + rest * "\t"
      end
    end

    new = File.expand_path(new)
    Open.mkdir File.dirname(new)
    CMD.cmd(:samtools, "view -h -b > '#{new}'", :in => fixed, "--threads" => cpus)
  end
end

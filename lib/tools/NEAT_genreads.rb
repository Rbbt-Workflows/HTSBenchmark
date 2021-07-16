require 'rbbt-util'
require 'rbbt/util/cmd'

module NEATGenreads

  Rbbt.claim Rbbt.software.opt.NEATGenreads, :install do
    {:git => "https://github.com/zstephens/neat-genreads.git",
     :commands => "chmod +x ${OPT_DIR}/NEATGenreads/gen_reads.py; [[ -f ${OPT_DIR}/bin/gen_reads.py ]] || ln -s ../NEATGenreads/genReads.py ${OPT_DIR}/bin/gen_reads.py"}
  end

  CMD.tool "gen_reads.py", Rbbt.software.opt.NEATGenreads, "gen_reads.py --help"

end

if __FILE__ == $0
  Log.with_severity 0 do
    ppp CMD.cmd("genReads.py", "-r #{Rbbt.share.organisms.Hsa.hg38["hg38.fa.gz"].find} ").read
  end
end

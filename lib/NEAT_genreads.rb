require 'rbbt-util'
require 'rbbt/util/cmd'

module NEATGenreads

  Rbbt.claim Rbbt.software.opt.NEATGenreads, :install do
    {:git => "https://github.com/zstephens/neat-genreads.git",
     :commands => "chmod +x ${OPT_DIR}/NEATGenreads/genReads.py; [[ -f ${OPT_DIR}/bin/genReads.py ]] || ln -s ../NEATGenreads/genReads.py ${OPT_DIR}/bin/genReads.py"}
  end

  CMD.tool "genReads.py", Rbbt.software.opt.NEATGenreads, "genReads.py --help"
end

if __FILE__ == $0
  Log.with_severity 0 do
    ppp CMD.cmd("genReads.py", "-r #{Rbbt.share.organisms.Hsa.hg38["hg38.fa.gz"].find} ").read
  end
end

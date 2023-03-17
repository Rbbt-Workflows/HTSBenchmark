require 'rbbt-util'
require 'rbbt/util/cmd'

module NEATGenreads

  Rbbt.claim Rbbt.software.opt.NEATGenreads, :install do
    commands=<<~EOF
chmod +x ${OPT_DIR}/NEATGenreads/gen_reads.py; sed -i 's/env source/env python/' ${OPT_DIR}/NEATGenreads/gen_reads.py; [[ -f ${OPT_DIR}/bin/gen_reads.py ]] || ln -s ../NEATGenreads/genReads.py ${OPT_DIR}/bin/gen_reads.py
    EOF
    {:git => "https://github.com/zstephens/neat-genreads.git",
     :commands => commands}
  end

  CMD.tool "gen_reads.py", Rbbt.software.opt.NEATGenreads, "gen_reads.py --help"

end

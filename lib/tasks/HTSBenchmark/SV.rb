require 'HTSBenchmark/SV'

module HTSBenchmark

  input :reference, :binary, "Reference file", nil, :nofile => true
  input :svs, :tsv, "SVs to apply to reference"
  extension 'fa.gz'
  task :SV_reference => :binary do |reference,svs|
    CMD.cmd("bgzip > #{self.tmp_path}", :in => HTSBenchmark.apply_SVs_reference(reference, svs))
    nil
  end

  input :svs, :tsv, "SVs to apply to reference" 
  input :mutations, :array, "Mutations to transpose"
  input :duplicate, :boolean, "Duplicate mutations or choose one", true
  task :SV_mutations => :array do |svs,mutations,duplicate|
    mutations = Open.read(mutations).split("\n") if Misc.is_filename?(mutations)
    mutation_translations = HTSBenchmark.transpose_mutations(svs, mutations)
    if duplicate
      mutation_translations.values.flatten.uniq
    else
      mutation_translations.values.collect{|v| v.shuffle.first }.uniq
    end
  end

end


module HTSBenchmark

  helper :subclones do |mutations,driver,info|
    genotype = [driver]
    info["num"].to_i.times do 
      genotype << mutations.gets.strip
    end

    clones = [genotype]
    info.each do |driver,subclones|
      next if driver == "num"
      clones += subclones(mutations, driver, subclones).collect{|list| genotype + list }
    end

    clones
  end

  dep :simulate_somatic_hg38
  input :evolution, :text, "Evolution"
  task :clonal => :yaml do |evolution|
    mutations = step(:simulate_somatic_hg38).join.path.open 
    evolution = YAML.load(evolution)
    driver = evolution.keys.first
    subclones = evolution[driver] 
    result = subclones mutations, driver, subclones 
    iii result
    result
  end
end

module HTSBenchmark
  dep :fp_positions
  dep :bundle
  dep_task :BAM_fp_images, HTS, :mutation_BAM_img, 
    :positions => :fp_positions, :reference => 'hg38', :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    varcal = dependencies.flatten.select{|dep| dep.task_name.to_s == 'fp_positions' }.first.step(:somatic_variant_calling)

    tumor = varcal.file('BAM.bam')
    normal = varcal.file('BAM_normal.bam')

    options[:tumor] = tumor
    options[:normal] = normal
    {:inputs => options, :jobname => jobname}
  end

  dep :fn_positions
  dep :bundle
  dep_task :BAM_fn_images, HTS, :mutation_BAM_img, 
    :positions => :fn_positions, :reference => 'hg38', :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    varcal = dependencies.flatten.select{|dep| dep.task_name.to_s == 'fn_positions' }.first.step(:somatic_variant_calling)

    tumor = varcal.file('BAM.bam')
    normal = varcal.file('BAM_normal.bam')

    options[:tumor] = tumor
    options[:normal] = normal
    {:inputs => options, :jobname => jobname}
  end


  dep :fp_positions
  dep :bundle
  dep_task :golden_BAM_fp_images, HTS, :mutation_BAM_img, 
    :positions => :fp_positions, :reference => 'hg38', :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    bundle = dependencies.flatten.select{|dep| dep.task_name.to_s == 'bundle' }.first

    normal = bundle.file('bundle/golden_BAM/normal_golden.bam')
    tumor = bundle.file('bundle/golden_BAM/tumor_golden.bam')

    options[:normal] = normal
    options[:tumor] = tumor
    {:inputs => options, :jobname => jobname}
  end

  dep :fn_positions
  dep :bundle
  dep_task :golden_BAM_fn_images, HTS, :mutation_BAM_img, 
    :positions => :fn_positions, :reference => 'hg38', :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    bundle = dependencies.flatten.select{|dep| dep.task_name.to_s == 'bundle' }.first

    normal = bundle.file('bundle/golden_BAM/normal_golden.bam')
    tumor = bundle.file('bundle/golden_BAM/tumor_golden.bam')

    options[:normal] = normal
    options[:tumor] = tumor
    {:inputs => options, :jobname => jobname}
  end

  dep :BAM_fp_images, :variant_caller => 'consensus_somatic_variants', :min_callers => 1
  dep :BAM_fn_images, :variant_caller => 'consensus_somatic_variants', :min_callers => 4
  dep :golden_BAM_fp_images, :variant_caller => 'consensus_somatic_variants', :min_callers => 1
  dep :golden_BAM_fn_images, :variant_caller => 'consensus_somatic_variants', :min_callers => 4
  dep :somatic_variant_calling, :variant_caller => 'combined_caller_vcfs'
  dep Sequence, :expanded_vcf, :vcf_file => :somatic_variant_calling
  dep :vcfeval do |jobname,options|
    jobs = %w(mutect2 somatic_sniper muse strelka).collect{|c| {:inputs => options.merge(:variant_caller => c)} }
    jobs += (1..4).to_a.collect{|min| {:inputs => options.merge(:variant_caller => 'consensus_somatic_variants', :min_callers => min)} }
    jobs
  end
  extension :pdf
  task :IGV_report => :binary do 
    require 'prawn'
    mutation_callers = step(:expanded_vcf).join.path.tsv.index :target => "Filter", :fields => "Genomic Mutation"
    fp_files = Path.setup(step(:BAM_fp_images).files_dir).glob("*.png").collect{|file| File.basename(file)}
    fn_files = Path.setup(step(:BAM_fn_images).files_dir).glob("*.png").collect{|file| File.basename(file)}

    vcfeval_jobs = dependencies.select{|d| d.task_name.to_s == "vcfeval"}

    pages = {}

    page = 1

    fp_info = {}
    fp_files.each do |file|
      name = file.sub('.png', '')
      iii name if name == "1:988599:G"
      callers = mutation_callers[name].nil? ? [] : mutation_callers[name].split(";").collect{|e| e.split("--").first}
      fp_info[name] = callers
    end

    fn_info = {}
    fn_files.each do |file|
      name = file.sub('.png', '')
      callers = mutation_callers[name].nil? ? [] : mutation_callers[name].split(";").collect{|e| e.split("--").first}
      fn_info[name] = callers
    end

    fp_max_callers = fn_max_callers = (fn_info.values + fp_info.values).flatten.uniq.length

    Prawn::Document.generate(self.tmp_path) do |pdf| 
      pdf.font_size 42
      pages["Performance"] = page
      pdf.text 'Performance'
      pdf.start_new_page
      page += 1

      pdf.font_size 20
      vcfeval_jobs.each do |job|
        name = job.recursive_inputs["variant_caller"]
        num = job.recursive_inputs["min_callers"]

        name = name + " #{num} callers" if name == 'consensus_somatic_variants'

        title = "Performance #{name}"
        pages[title] = page
        pdf.text title
        pdf.stroke_horizontal_rule
        pdf.move_down 30
        tsv = job.load
        k = tsv.keys.last
        tsv.fields.each do |f|
          pdf.text "#{f}: #{tsv[k][f]}"
        end
        pdf.start_new_page
        page += 1
      end


      pdf.font_size 42
      pdf.text 'False Positives'
      pages["False Positives"] = page

      pdf.start_new_page
      page += 1

      pdf.font_size 20

      (0..fp_max_callers).to_a.reverse.each do |num|
        muts = fp_info.select{|name,callers| callers.length == num }.length
        next if muts == 0
        pdf.text "False positives called by #{num} callers (#{muts})"
        pages["False positives called by #{num} callers (#{muts})"] = page
        pdf.start_new_page
        page += 1
        fp_info.select{|name,callers| callers.length == num }.sort{|a,b| Misc.genomic_location_cmp_strict(a[0],b[0])}.each do |name,callers|
          pdf.text "False positive: " + name
          pdf.text "Called by: " + callers.sort * ", "
          pdf.image File.join(step(:BAM_fp_images).files_dir, name + '.png'), at: [-20, 600], width: 300
          pdf.image File.join(step(:golden_BAM_fp_images).files_dir, name + '.png'), at: [280, 600], width: 300
          pdf.start_new_page
          page += 1
        end
      end

      pdf.text 'False Negatives'
      pages["False Negatives"] = page

      pdf.start_new_page
      page += 1

      (0..fn_max_callers).to_a.each do |num|
        muts = fn_info.select{|name,callers| callers.length == num }.length
        next if muts == 0
        pdf.text "False negatives missed by #{fp_max_callers - num} callers (#{muts})"
        pages["False negatives missed by #{fp_max_callers - num} callers (#{muts})"] = page
        pdf.start_new_page
        page += 1
        fn_info.select{|name,callers| callers.length == num }.sort{|a,b| Misc.genomic_location_cmp_strict(a[0],b[0])}.each do |name,callers|
          pdf.text "False negative: " + name
          pdf.text "Called by: " + callers.sort * ", "
          pdf.image File.join(step(:BAM_fn_images).files_dir, name + '.png'), at: [-20, 600], width: 300
          pdf.image File.join(step(:golden_BAM_fn_images).files_dir, name + '.png'), at: [280, 600], width: 300
          pdf.start_new_page
          page += 1
        end
      end

      pdf.outline.define do
        pages.each do |name, page|
          pdf.outline.section name, destination: page
        end
      end
    end
  end
end

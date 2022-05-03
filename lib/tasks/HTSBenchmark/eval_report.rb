module HTSBenchmark
  input :mutation_eval_class, :select, "Class of mutation evaluation false positive or false negative (FP or FN)", "FN", :select_options => %w(FN FP), :required => true
  dep :eval_bundle
  dep Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    eval_bundle = dependencies.flatten.select{|d| d.task_name == :eval_bundle }.first
    type = options[:mutation_eval_class].downcase == "fn" ? "fn" : "fp"
    options[:vcf_file]= eval_bundle.file("output/#{type}.vcf.gz")
    {:inputs => options, :jobname => jobname}
  end
  dep_task :eval_bundle_mutation_images, HTS, :mutation_BAM_img, :positions => :genomic_mutations, :reference => :placeholder, :normal => :placeholder, :tumor => :placeholder do |jobname,options,dependencies|
    eval_bundle = dependencies.flatten.flatten.first
    bam = eval_bundle.step(:mutect2_pre).inputs[:tumor]
    bam_normal = eval_bundle.step(:mutect2_pre).inputs[:normal]
    options[:tumor] = bam
    options[:normal] = bam_normal
    options[:reference] = bam.recursive_inputs[:reference]
    {:inputs => options, :jobname => jobname}
  end

  dep :eval_bundle
  dep :eval_bundle_mutation_images, :mutation_eval_class => "FP"
  dep :eval_bundle_mutation_images, :mutation_eval_class => "FN"
  extension :pdf
  task :eval_bundle_report => :binary do 
    require 'prawn'
    fp_imgs, fn_imgs = dependencies.select{|d| d.task_name.to_s == "eval_bundle_mutation_images"}

    fp_files = Path.setup(fp_imgs.files_dir).glob("*.png").collect{|file| file}
    fn_files = Path.setup(fn_imgs.files_dir).glob("*.png").collect{|file| file}

    vcfeval = step(:eval_bundle)

    pages = {}

    page = 1

    Prawn::Document.generate(self.tmp_path) do |pdf| 
      pdf.font_size 42
      pages["Performance"] = page
      pdf.text 'Performance'
      pdf.start_new_page
      page += 1

      pdf.font_size 20

      name = step("run_bundle").inputs["variant_caller"]

      name = name + " #{num} callers" if name == 'consensus_somatic_variants'

      title = "Performance #{name}"
      pages[title] = page
      pdf.text title
      pdf.stroke_horizontal_rule
      pdf.move_down 30
      tsv = vcfeval.load
      k = tsv.keys.last
      tsv.fields.each do |f|
        pdf.text "#{f}: #{tsv[k][f]}"
      end
      pdf.start_new_page
      page += 1


      log :FP, "Preparing false positive report"
      pdf.font_size 42
      pdf.text 'False Positives'
      pdf.font_size 20
      pages["False Positives"] = page
      pdf.start_new_page
      page += 1

      TSV.traverse fp_files, :type => :array, :bar => self.progress_bar("Processing FP images") do |file|
        name = File.basename(file, '.png')
        pdf.text "False positive: " + name
        pdf.image file, at: [-20, 600], width: 300
        pdf.start_new_page
        page += 1
      end

      log :FN, "Preparing false negative report"
      pdf.font_size 42
      pdf.text 'False Negatives'
      pdf.font_size 20
      pages["False Negatives"] = page
      pdf.start_new_page
      page += 1

      TSV.traverse fn_files, :type => :array, :bar => self.progress_bar("Processing FN images") do |file|
        name = File.basename(file, '.png')
        pdf.text "False negative: " + name
        pdf.image file, at: [-20, 600], width: 300
        pdf.start_new_page
        page += 1
      end
      pdf.outline.define do
        pages.each do |name, page|
          pdf.outline.section name, destination: page
        end
      end
    end
  end

  dep :eval_population
  dep_task :eval_population_report, HTSBenchmark, :eval_bundle_report, :evolution => :placeholder, "HTSBenchmark#eval_bundle" => :eval_population, :not_overriden => true

  dep :eval_simulate_population
  dep_task :eval_simulate_population_report, HTSBenchmark, :eval_bundle_report, "HTSBenchmark#eval_bundle" => :eval_simulate_population, :not_overriden => true

end

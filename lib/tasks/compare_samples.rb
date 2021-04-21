module HTSBenchmark

  input :sample_input, :string, "Sample code for input sample"
  input :sample_truth, :string, "Sample code for truth sample"
  dep Sample, :vcf_file do |jobname,options|
    options.values_at(:sample_input, :sample_truth).
      collect{|s| {:jobname => s, :inputs => options} }
  end
  dep_task :cmp_vcfeval, HTS, :vcfeval, :reference => :placeholder do |jobname,options,dependencies|

    sample = options[:sample]

    input, truth = dependencies.flatten

    reference = input.recursive_inputs[:reference]

    options[:input_vcf] = input
    options[:truth_vcf] = truth
    options[:reference] = reference
    {:inputs => options, :jobname => sample}
  end

  dep :cmp_vcfeval
  dep_task :cmp_fp_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fp.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end

  dep :cmp_vcfeval
  dep_task :cmp_fn_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/fn.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end

  dep :cmp_vcfeval
  dep_task :cmp_tp_positions, Sequence, :genomic_mutations, :vcf_file => :placeholder do |jobname,options,dependencies|
    vcfeval = dependencies.flatten.first
    options = options.merge(:vcf_file => vcfeval.file('output/tp.vcf.gz'))
    {:inputs => options, :jobname => jobname}
  end


  dep :cmp_fp_positions
  dep :cmp_fn_positions
  input :max_mut_imgs, :integer, "Max number of mutations to show", nil
  task :cmp_discordant_positions => :array do |max_mut_imgs|
    if max_mut_imgs
      Misc.sort_mutations_strict(dependencies.collect{|d| d.load}.flatten.shuffle[0..max_mut_imgs-1])
    else
      dependencies.collect{|d| d.load}.flatten
    end
  end

  dep :cmp_tp_positions
  input :max_mut_imgs, :integer, "Max number of mutations to show", nil
  task :cmp_concordant_positions => :array do |max_mut_imgs|
    if max_mut_imgs
      Misc.sort_mutations_strict(dependencies.collect{|d| d.load}.flatten.shuffle[0..max_mut_imgs-1])
    else
      dependencies.collect{|d| d.load}.flatten
    end
  end


  dep :cmp_discordant_positions
  dep_task :cmp_BAM_discordant_images, HTS, :mutation_BAM_img, 
    :positions => :cmp_discordant_positions, :reference => :placeholder, :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    input, truth = dependencies.first.rec_dependencies.select{|dep| dep.task_name.to_s == 'mutect2' && dep.workflow.to_s == "Sample" }

    options[:reference] = input.recursive_inputs[:reference]

    [input, truth].collect{|varcal|
      tumor = varcal.step(:BAM)
      normal = varcal.step(:BAM_normal)

      options[:tumor] = tumor
      options[:normal] = normal
      {:inputs => options.dup, :jobname => varcal.clean_name}
    }
  end

  dep :cmp_concordant_positions
  dep_task :cmp_BAM_concordant_images, HTS, :mutation_BAM_img, 
    :positions => :cmp_concordant_positions, :reference => :placeholder, :depth => 80,
    :tumor => :placeholder, :normal => :placeholder do  |jobname,options,dependencies|

    input, truth = dependencies.first.rec_dependencies.select{|dep| dep.task_name.to_s == 'mutect2' && dep.workflow.to_s == "Sample" }

    options[:reference] = input.recursive_inputs[:reference]

    [input, truth].collect{|varcal|
      tumor = varcal.step(:BAM)
      normal = varcal.step(:BAM_normal)

      options[:tumor] = tumor
      options[:normal] = normal
      {:inputs => options.dup, :jobname => varcal.clean_name}
    }
  end


  input :type_of_comparison, :select, "Show discordant or concordant positions", 'discordant', :select_options => %w(discordant concordant)
  dep :cmp_BAM_discordant_images do |jobname, options|
    if options[:type_of_comparison] == "concordant"
      {:inputs => options, :jobname => jobname, :task => :cmp_BAM_concordant_images}
    else
      {:inputs => options, :jobname => jobname}
    end
  end
  extension :pdf
  task :cmp_IGV_report => :binary do 
    require 'prawn'

    jimages = dependencies.first
    positions = jimages.dependencies.first.load

    jimages1, jimages2 = jimages.dependencies.select{|dep| dep.task_name.to_s == "mutation_BAM_img"}

    concordant = self.recursive_inputs[:type_of_comparison] == "concordant"

    if concordant
      tps = jimages.step(:cmp_tp_positions).load
      fps = fns = nil
      tps = tps & positions
    else
      tps = nil
      fps = jimages.step(:cmp_fp_positions).load
      fns = jimages.step(:cmp_fn_positions).load
      fps = fps & positions
      fns = fns & positions
    end


    pages = {}
    page = 1
    Prawn::Document.generate(self.tmp_path) do |pdf| 
      pdf.font_size 42

      if tps
        pages["True positives"] = page
        pdf.text 'True positives'
        pdf.start_new_page
        page += 1

        pdf.font_size 20
        tps.each do |tp|
          pdf.text "True positive: " + tp
          pdf.image File.join(jimages1.files_dir, tp + '.png'), at: [-20, 600], width: 300
          pdf.image File.join(jimages2.files_dir, tp + '.png'), at: [280, 600], width: 300
          pdf.start_new_page
          page += 1
        end
      end


      if fps
        pages["False positives"] = page
        pdf.text 'False positives'
        pdf.start_new_page
        page += 1

        pdf.font_size 20
        fps.each do |fp|
          pdf.text "False positive: " + fp
          pdf.image File.join(jimages1.files_dir, fp + '.png'), at: [-20, 600], width: 300
          pdf.image File.join(jimages2.files_dir, fp + '.png'), at: [280, 600], width: 300
          pdf.start_new_page
          page += 1
        end
      end

      if fns
        pdf.font_size 42
        pages["False negatives"] = page
        pdf.text 'False negatives'
        pdf.start_new_page
        page += 1

        pdf.font_size 20
        fns.each do |fn|
          pdf.text "False negative: " + fn
          pdf.image File.join(jimages1.files_dir, fn + '.png'), at: [-20, 600], width: 300
          pdf.image File.join(jimages2.files_dir, fn + '.png'), at: [280, 600], width: 300
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

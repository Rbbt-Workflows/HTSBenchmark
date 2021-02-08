module HTSBenchmark

  input :sample1, :string, "Sample code"
  input :sample2, :string, "Sample code"
  dep Sample, :mutect2 do |jobname,options|
    options.values_at(:sample1, :sample2).
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

  dep :cmp_fp_positions
  dep :cmp_fn_positions
  task :cmp_discordant_positions => :array do
    dependencies.collect{|d| d.load}.flatten
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
      {:inputs => options, :jobname => varcal.clean_name}
    }
  end

  dep :cmp_BAM_discordant_images
  extension :pdf
  task :cmp_IGV_report => :binary do 
    require 'prawn'

    jimages = step(:cmp_BAM_discordant_images)
    jimages1, jimages2 = jimages.dependencies.select{|dep| dep.task_name.to_s == "mutation_BAM_img"}

    fps = jimages.step(:cmp_fp_positions).load
    fns = jimages.step(:cmp_fn_positions).load

    Prawn::Document.generate(self.tmp_path) do |pdf| 
      pdf.font_size 42
      pages["False positives"] = page
      pdf.text 'False positives'
      pdf.start_new_page
      page += 1

      pdf.font_size 20
      fps.each do |fp|
        pdf.text "False positive: " + fp
        pdf.image File.join(jimages1.files_dir, fp + '.png'), at: [-20, 600], width: 300
        pdf.image File.join(jimages2.files_dir, fp + '.png'), at: [-20, 600], width: 300
        pdf.start_new_page
        page += 1
      end

      pdf.font_size 42
      pages["False negatives"] = page
      pdf.text 'False negatives'
      pdf.start_new_page
      page += 1

      pdf.font_size 20
      fns.each do |fn|
        pdf.text "False negative: " + fp
        pdf.image File.join(jimages1.files_dir, fn + '.png'), at: [-20, 600], width: 300
        pdf.image File.join(jimages2.files_dir, fn + '.png'), at: [-20, 600], width: 300
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


end

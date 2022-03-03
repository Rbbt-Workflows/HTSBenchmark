require File.join(File.expand_path(File.dirname(__FILE__)), '..', 'test_helper.rb')
require 'rbbt-util'
require 'HTSBenchmark/sliceref'
require 'rbbt/workflow'

Workflow.require_workflow "HTS"
class TestSliceRef < Test::Unit::TestCase
  def test_slice_fasta
    reference =<<-EOF
>chr1
12345678901234567890
>chr2
abcdefghijabcdefghij
>chr3
ABCDEFGHIJABCDEFGHIJ
>chr4
12345678901234567890
>chr5
ABCDEFGHIJABCDEFGHIJ
    EOF

    ranges = {}
    ranges["chr1"] = [[2,5], [8,9]]
    ranges["chr2"] = [[2,5], [8,9]]
    ranges["chr3"] = [[2,5], [8,9]]
    ranges["chr4"] = [[2,5], [8,9]]
    ranges["chr5"] = [[2,5], [8,9]]
    TmpFile.with_file(reference) do |ref|
      TmpFile.with_file do |newref|
        HTSBenchmark.slice_fasta(ref, newref, ranges)
        ppp Open.read(newref)
      end
    end
  end

  def test_slice_vcf
    vcf =<<-EOF
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	1	.	A	T	.	PASS	GT	1/0
chr1	4	.	A	T	.	PASS	GT	1/0
chr2	1	.	A	T	.	PASS	GT	1/0
chr2	8	.	A	T	.	PASS	GT	1/0
chr3	2	.	A	T	.	PASS	GT	1/0
chr3	9	.	A	T	.	PASS	GT	1/0
    EOF

    ranges = {}
    ranges["chr1"] = [[2,5], [8,9]]
    ranges["chr2"] = [[2,5], [8,9]]
    ranges["chr3"] = [[2,5], [8,9]]
    ranges["chr4"] = [[2,5], [8,9]]
    ranges["chr5"] = [[2,5], [8,9]]
    Log.with_severity 0 do
    TmpFile.with_file(vcf, false, :extension => 'vcf') do |vcf|
      CMD.cmd("bgzip #{vcf}")
      vcf += '.gz'
      TmpFile.with_file(nil, false, :extension => 'vcf') do |newvcf|
        TmpFile.with_file(nil, false, :extension => 'bed') do |bed|
          Open.write(bed) do |f|
            ranges.each do |chr,list|
              list.each do |s,e|
                f.puts [chr, s, e] * "\t"
              end
            end
          end
          HTSBenchmark.slice_vcf(vcf, newvcf, bed)
          ppp Open.read(newvcf)
        end
      end
    end
    end
  end

  def test_restore_sliced_vcf
    vcf_text =<<-EOF
##fileformat=VCFv4.1
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE
chr1	10	.	A	T	.	PASS	GT	1/0
chr1	20	.	A	T	.	PASS	GT	1/0
chr1	30	.	A	T	.	PASS	GT	1/0
    EOF

    ranges = {}
    ranges["chr1"] = [[5,15], [15,25], [25, 35]]
    ranges["chr2"] = [[2,5], [8,9]]
    ranges["chr3"] = [[2,5], [8,9]]
    ranges["chr4"] = [[2,5], [8,9]]
    ranges["chr5"] = [[2,5], [8,9]]
    Log.with_severity 0 do

    TmpFile.with_file(vcf_text, false, :extension => 'vcf') do |vcf|
      CMD.cmd("bgzip #{vcf}")
      vcf += '.gz'
      TmpFile.with_file(nil, false, :extension => 'vcf') do |newvcf|
        TmpFile.with_file(nil, false, :extension => 'bed') do |bed|
          Open.write(bed) do |f|
            ranges.each do |chr,list|
              list.each do |s,e|
                f.puts [chr, s, e] * "\t"
              end
            end
          end
          HTSBenchmark.slice_vcf(vcf, newvcf, bed)

          CMD.cmd("bgzip #{newvcf}")
          newvcf += '.gz'

          TmpFile.with_file(nil, false, :extension => 'vcf') do |restored_vcf|
            HTSBenchmark.restore_sliced_vcf(newvcf, restored_vcf, bed)
            assert_equal vcf_text, Open.read(restored_vcf)
          end
        end
      end
    end
    end
  end
end


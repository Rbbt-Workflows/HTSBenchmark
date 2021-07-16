require File.join(File.expand_path(File.dirname(__FILE__)), '..', 'test_helper.rb')
require 'rbbt-util'
require 'HTSBenchmark/sliceref'

class TestSliceRef < Test::Unit::TestCase
  def _test_slice_fasta
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
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
chr1	1	.	a	A	.
chr1	4	.	b	B	.
chr2	1	.	c	C	.
chr2	8	.	d	D	.
chr3	2	.	e	E	.
chr3	9	.	f	G	.
    EOF

    ranges = {}
    ranges["chr1"] = [[2,5], [8,9]]
    ranges["chr2"] = [[2,5], [8,9]]
    ranges["chr3"] = [[2,5], [8,9]]
    ranges["chr4"] = [[2,5], [8,9]]
    ranges["chr5"] = [[2,5], [8,9]]
    TmpFile.with_file(vcf) do |vcf|
      TmpFile.with_file do |newvcf|
        HTSBenchmark.slice_vcf(vcf, newvcf, ranges)
        ppp Open.read(newvcf)
      end
    end
  end
end


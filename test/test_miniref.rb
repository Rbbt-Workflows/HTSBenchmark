
require File.join(File.expand_path(File.dirname(__FILE__)), '../test', 'test_helper.rb')
require 'rbbt/workflow'
require 'miniref'

class TestMinify < Test::Unit::TestCase
  def setup
    @log = Log.severity
    Log.severity = 0
  end

  def teardown
    Log.severity = @log if @log
  end

  def test_sizes_from_positions
    Workflow.require_workflow "PCAWG"
    positions = Study.setup("Bladder-TCC").genomic_mutations
    sizes = HTSBenchmark.calculate_sizes(positions)
    iii sizes
  end

  def test_minify_fasta
    fasta = Rbbt.share.organisms.Hsa.b37["b37.fa.gz"].find
    TmpFile.with_file do |tmp|
      HTSBenchmark.minify_fasta(fasta, tmp, "1" => 1000)
      iif Open.read(tmp)
    end
  end

  def test_minify_vcf
    vcf = Rbbt.share.organisms.Hsa.b37.known_sites["1000G_phase1.indels.vcf.gz"].find
    TmpFile.with_file do |tmp|
      HTSBenchmark.minify_vcf(vcf, tmp, "1" => 100000)
      ppp Open.read(tmp)
    end
  end

end


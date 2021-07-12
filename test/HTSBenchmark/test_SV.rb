require File.join(File.expand_path(File.dirname(__FILE__)), '../test_helper.rb')
require 'rbbt-util'
require 'HTSBenchmark/SV'

class TestSV < Test::Unit::TestCase
  def test_collect_fragment
    input =<<-EOF
>chr1
AAACCCTTTGGG
>chr2
AAAAAACCCCC
CTTTTTTGGGGGG
>chr3
TTTTTTGGGGGGA
AAAAACCCCCC
    EOF

    input = StringIO.new input

    ranges = {
      "chr1" => [[3,5], [7,10]],
      "chr2" => [[3,5], [7,10]],
      "chr3" => [[3,5], [7,10]],
    }

    fragments = HTSBenchmark.collect_fragments(input, ranges)
    assert_equal fragments["chr1"][[3,5,nil]], "ACC"
    assert_equal fragments["chr2"][[3,5,nil]], "AAA"
  end

  def test_remove_fragment
    input =<<-EOF
>chr1
AAACCC
TTTGGG
>chr2
AAAAAACC
CCCCTTTTTTGG
GGGG
>chr3
TTTTTTGGGGG
GAAAAAACCCCCC
    EOF

    target =<<-EOF
>chr1
ACTTGGG
>chr2
AAAAAACCCCCCTTTTTTGGGGGG
>chr3
TTTTTTGGGGGGAAAAAACCCCCC
    EOF


    ranges = {
      "chr1" => [[2,4], [6,7]],
    }

    output = HTSBenchmark.remove_fragments(StringIO.new(input), ranges)

    assert_equal target, output.read
  end

  def test_insert_fragment
    input =<<-EOF
>chr1
AAA
CCCTTT
GGG
>chr2
AAAAAACCCCCC
TTTTTTGGGGGG
>chr3
TTTTTTGG
GGGGAAAAAA
CCCCCC
    EOF

    target =<<-EOF
>chr1
ATTTTAACGGCCTTTGGG
>chr2
AAAAAACCCCCCTTTTTTGGGGGG
>chr3
TTTTTTGGGGGGAAAAAACCCCCC
    EOF


    inserts = {
      "chr1" => [[2,'TTTT'], [5, 'GG']]
    }

    output = HTSBenchmark.insert_fragments(StringIO.new(input), inserts)

    assert_equal target, output.read
  end

  def test_shift_ranges_by_deletes
    ranges = {}

    ranges_txt = <<-EOF
chr1:1000-1100
chr1:1300-1400
chr2:1000:A
chr2:1300:A
chr2:2300:A
    EOF

    ranges_txt.split("\n").each{|r|
      chr, range = r.split(":")
      start, eend = range.split("-")
      eend = start if eend.nil?
      ranges[chr] ||= []
      ranges[chr] << [start.to_i, eend.to_i, r]
    }

    deletes = {}

    deletes_txt = <<-EOF
chr1:100-200
chr1:500-700
chr1:2000-3000
chr2:100-200
chr2:500-700
chr2:1100-2000
    EOF

    deletes_txt.split("\n").each{|r|
      chr, range = r.split(":")
      start, eend = range.split("-")
      eend = start if eend.nil?
      deletes[chr] ||= []
      deletes[chr] << [start.to_i, eend.to_i, r]
    }

    new_ranges = HTSBenchmark.shift_ranges_by_deletes ranges, deletes

    id_ranges = {}
    new_ranges.values.flatten(1).each do |values|
      start,eend,id = values
      id_ranges[id] = [start, eend]
    end

    assert_equal [998, 1098], id_ranges["chr1:1300-1400"]
    assert_equal [nil, nil], id_ranges["chr2:1300:A"]

  end

  def test_shift_ranges_by_inserts
    ranges = {}

    ranges_txt = <<-EOF
chr1:1000-1100
chr1:1300-1400
chr2:1000:A
chr2:1300:A
chr2:2300:A
    EOF

    ranges_txt.split("\n").each{|r|
      chr, range = r.split(":")
      start, eend = range.split("-")
      eend = start if eend.nil?
      ranges[chr] ||= []
      ranges[chr] << [start.to_i, eend.to_i, r]
    }

    inserts = {}

    inserts_txt = <<-EOF
chr1:100:100
chr1:1350:100
chr2:100:100
chr2:500:200
    EOF

    inserts_txt.split("\n").each{|r|
      chr, start, length = r.split(":")
      inserts[chr] ||= []
      inserts[chr] << [start.to_i, length.to_i, r]
    }

    new_ranges = HTSBenchmark.shift_ranges_by_inserts ranges, inserts

    id_ranges = {}
    new_ranges.values.flatten(1).each do |values|
      start,eend,id = values
      id_ranges[id] = [start, eend]
    end

    assert_equal [1100, 1200], id_ranges["chr1:1000-1100"]
    assert_equal [nil, nil], id_ranges["chr1:1300-1400"]
    assert_equal [1600, 1600], id_ranges["chr2:1300:A"]

  end

  def test_duplicate_positions
    positions = {}

    positions_txt = <<-EOF
chr1:1000:A
chr1:1300:A
chr1:2300:A
    EOF

    positions_txt.split("\n").each{|r|
      chr, pos = r.split(":")
      positions[chr] ||= []
      positions[chr] << [pos.to_i, pos.to_i, r]
    }

    duplications = {}

    duplications_txt = <<-EOF
chr1:1200-1400:chr1:200
chr1:900-1100:chr1:1200:inverse
chr1:1200-1400:chr2:200
    EOF

    duplications_txt.split("\n").each{|r|
      chr, range, tchr, tpos, inverse = r.split(":")
      inverse = inverse == 'inverse'
      start, eend = range.split("-")
      duplications[chr] ||= []
      duplications[chr] << [start.to_i, eend.to_i, tchr, tpos.to_i, inverse, r]
    }

    new_positions = HTSBenchmark.duplicate_positions positions, duplications

    new_positions.each do |chr,pos,id|
      positions[chr] ||= []
      positions[chr] << [pos, pos, id]
    end

    assert positions["chr2"].include?([300, 300, "chr1:1300:A"])

  end

  def test_apply_SVs
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


    mutations = <<-EOF.split("\n")
chr1:1:chr1-1(1)
chr2:8:chr2-8(h)
chr2:15:chr2-15(e)
    EOF

    svs_txt =<<-EOF
#: :type=:list#:sep=,
#ID,Type,Chromosome,Start,End,Target chromosome,Target start,Target end
DUP.1,INS,chr1,1,6,chr2,8
DUP.2,INS,chr1,9,10,chr2,6
DEL.1,DEL,chr1,1,6
DEL.2,DEL,chr2,15,20
INV.1,INV,chr2,6,9,chr3,2
    EOF

    svs = TSV.open(StringIO.new(svs_txt))

    reference, mutation_translation = TmpFile.with_file(reference) do |reference|
      reference_io  = HTSBenchmark.apply_SVs_reference reference, svs
      mutation_translation = HTSBenchmark.transpose_mutations svs, mutations
      [reference_io.read, mutation_translation]
    end

    assert reference.include?("abcde90fg123456hijabcd")
    assert_equal ["chr2:10:chr1-1(1)"], mutation_translation["chr1:1:chr1-1(1)"]
  end

  def ___test_apply_SVs
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

    svs_txt1 =<<-EOF
#: :type=:list#:sep=,
#ID,Type,Chromosome,Start,End,Target chromosome,Target start,Target end
DUP.1,INS,chr1,1,6,chr2,8
    EOF

    svs_txt2 =<<-EOF
#: :type=:list#:sep=,
#ID,Type,Chromosome,Start,End,Target chromosome,Target start,Target end
DUP.2,INS,chr1,9,10,chr2,6
    EOF

    svs = TSV.open(StringIO.new(svs_txt))

    reference, mutation_translation = TmpFile.with_file(reference) do |reference|
      reference_io  = HTSBenchmark.apply_SVs_reference reference, svs
      mutation_translation = HTSBenchmark.transpose_mutations reference, svs, mutations
      [reference_io.read, mutation_translation]
    end

    assert reference.include?("abcde90fg123456hijabcd")
    assert_equal ["chr2:10:chr1-1(1)"], mutation_translation["chr1:1:chr1-1(1)"]
  end
end


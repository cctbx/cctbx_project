from iotbx import bioinformatics

import unittest

sequence1 = "VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS"
sequence2 = "LVKFYGRKFVLLFMDQKTFDKYESKLETSGYEKFFIFCMASPISPA"
name1 = "gi|159164330|pdb|2E6P|A"
name2 = "Q8E5Q5_STRA3"
description_tmplt = "Chain %s, Solution Structure Of The Ig-Like Domain"

class test_sequence(unittest.TestCase):

  def testEquality(self):

    seq1 = bioinformatics.sequence(
      sequence = sequence1,
      name = name1
      )

    seq2 = bioinformatics.sequence(
      sequence = sequence1,
      name = name2
      )

    fasta1 = bioinformatics.fasta_sequence(
      sequence = sequence1,
      name = name2,
      description = description_tmplt % "A"
      )

    fasta2 = bioinformatics.fasta_sequence(
      sequence = sequence1,
      name = name2,
      description = description_tmplt % "B"
      )

    pir1 = bioinformatics.pir_sequence(
      sequence = sequence1,
      name = name2,
      type = "P1",
      description = description_tmplt % "C"
      )

    self.assertEqual( seq1, seq2 )
    self.assertEqual( seq1, fasta1 )
    self.assertEqual( seq1, fasta2 )
    self.assertEqual( seq1, pir1 )

    self.assertEqual( seq2, fasta1 )
    self.assertEqual( seq2, fasta2 )
    self.assertEqual( seq2, pir1 )

    self.assertEqual( fasta1, fasta2 )
    self.assertEqual( fasta1, pir1 )

    self.assertEqual( fasta2, pir1 )


  def testNonEquality(self):

    seq1 = bioinformatics.sequence(
      sequence = sequence1,
      name = name1
      )

    seq2 = bioinformatics.sequence(
      sequence = sequence2,
      name = name1
      )

    seq3 = bioinformatics.sequence(
      sequence = sequence1 + "A",
      name = name1
      )

    fasta1 = bioinformatics.fasta_sequence(
      sequence = sequence2,
      name = name1,
      description = description_tmplt % "A"
      )

    pir1 = bioinformatics.pir_sequence(
      sequence = sequence2,
      name = name1,
      type = "P1",
      description = description_tmplt % "A"
      )

    self.assert_( seq1 != seq2 )
    self.assert_( seq1 != seq3 )
    self.assert_( seq1 != fasta1 )
    self.assert_( seq1 != pir1 )


class test_fasta_sequence(unittest.TestCase):

  def setUp(self):

    self.short = bioinformatics.fasta_sequence(
      sequence = sequence1,
      name = name1,
      description = description_tmplt % "A"
      )
    self.long = bioinformatics.fasta_sequence(
      sequence = sequence1 * 4,
      name = name1,
      description = description_tmplt % "A"
      )


  def testFormat(self):

    self.assertEqual(
      self.short.format( 70 ),
      ">gi|159164330|pdb|2E6P|A Chain A, Solution Structure Of The Ig-Like Do\n"
      + "VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS"
      )

    self.assertEqual(
      self.short.format( 20 ),
      ">gi|159164330|pdb|2E\n"
      + "VVKMDGRKHRLILPEAKVQD\n"
      + "SGEFECRTEGVSAFFGVTVQ\n"
      + "DPSGPS"
      )


  def testStr(self):

    self.assertEqual(
      str( self.short ),
      ">gi|159164330|pdb|2E6P|A Chain A, Solution Structure Of The Ig-Like Do\n"
      + "VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS"
      )

    self.assertEqual(
      str( self.long ),
      ">gi|159164330|pdb|2E6P|A Chain A, Solution Structure Of The Ig-Like Do\n"
      + "VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPSVVKMDGRKHRLILPEAKVQDSGEF\n"
      + "ECRTEGVSAFFGVTVQDPSGPSVVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPSVV\n"
      + "KMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS"
      )


class test_pir_sequence(unittest.TestCase):

  def setUp(self):

    self.short = bioinformatics.pir_sequence(
      sequence = sequence1,
      name = name1,
      type = "P1",
      description = description_tmplt % "A"
      )
    self.long = bioinformatics.pir_sequence(
      sequence = sequence1 * 4,
      name = name1,
      type = "P1",
      description = description_tmplt % "A"
      )


  def testFormat(self):

    self.assertEqual(
      self.short.format( 70 ),
      ">P1;gi|159164330|pdb|2E6P|A\n"
      + "Chain A, Solution Structure Of The Ig-Like Domain\n"
      + "  VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS*"
      )

    self.assertEqual(
      self.short.format( 20 ),
      ">P1;gi|159164330|pdb\n"
      + "Chain A, Solution St\n"
      + "  VVKMDGRKHRLILPEAKV\n"
      + "  QDSGEFECRTEGVSAFFG\n"
                + "  VTVQDPSGPS*"
      )


  def testStr(self):

    self.assertEqual(
      str( self.short ),
      ">P1;gi|159164330|pdb|2E6P|A\n"
      + "Chain A, Solution Structure Of The Ig-Like Domain\n"
      + "  VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS*"
      )

    self.assertEqual(
      str( self.long ),
      ">P1;gi|159164330|pdb|2E6P|A\n"
      + "Chain A, Solution Structure Of The Ig-Like Domain\n"
      + "  VVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPSVVKMDGRKHRLILPEAKVQDSG\n"
      + "  EFECRTEGVSAFFGVTVQDPSGPSVVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSG\n"
      + "  PSVVKMDGRKHRLILPEAKVQDSGEFECRTEGVSAFFGVTVQDPSGPS*"
      )


class test_midline(unittest.TestCase):

  def setUp(self):

    self.midline = bioinformatics.midline( identical = "#" )


  def testConservationCode(self):

    self.assertEqual(
      self.midline.conservation_code( [ "A", "A", "A", "A" ] ),
      self.midline.identical
      )
    self.assertEqual(
      self.midline.conservation_code( [ "A", "A", "A", "B" ] ),
      self.midline.differ
      )


  def testMidline(self):

    self.assertEqual(
      self.midline.compare( [ "ABCD-EFGH", "ABD--EKGH" ] ),
      "##   # ##"
      )
    self.assertEqual(
      self.midline.compare( [] ),
      ""
      )

ali_1hml = "-KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYDTQAIVENN--ESTEYGLFQISNKLWCKSSQVPQSR"
ali_1hfy = "-EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNN--DSTEYGLFQINNKIWCKDDQNPHSR"
ali_1ghl = "GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFNTGATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK"
ali_1lz3 = "-KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFDTHATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK"
ali_empt = "---------------------------------------------------------------------------"


class test_alignment(unittest.TestCase):

  def setUp(self):

    self.alignment1 = bioinformatics.alignment(
      alignments = [ ali_1hml, ali_1hfy ],
      names = [ "1hml", "1hfy" ]
      )
    self.alignment2 = bioinformatics.alignment(
      alignments = [ ali_1hml, ali_1hfy, ali_1ghl, ali_1lz3 ],
      names = [ "1hml", "1hfya", "1ghla", "1lz3" ]
      )
    self.alignment3 = bioinformatics.alignment(
      alignments = [],
      names = []
      )
    self.alignment4 = bioinformatics.alignment(
      alignments = [ ali_1hfy, ali_1ghl, ali_empt ],
      names = [ "1hfya", "1ghla", "empty" ]
      )


  def test_error(self):

    self.assertRaises(
      ValueError,
      bioinformatics.alignment,
      [ ali_1hml, ali_1hml + "A" ],
      [ "1hml", "1hfy" ]
      )
    self.assertRaises(
      ValueError,
      bioinformatics.alignment,
      [ ali_1hml, ali_1hml ],
      [ "1hml", "1hfy", "1ghla" ]
      )


  def test_identity_count(self):

    self.assertEqual( self.alignment1.identity_count(), 49 )
    self.assertEqual( self.alignment2.identity_count(), 21 )
    self.assertEqual( self.alignment3.identity_count(), 0 )
    self.assertEqual( self.alignment4.identity_count(), 0 )


  def test_identity_fraction(self):

    self.assertAlmostEqual( self.alignment1.identity_fraction(), 0.700, 3 )
    self.assertAlmostEqual( self.alignment2.identity_fraction(), 0.300, 3 )
    self.assertAlmostEqual( self.alignment3.identity_fraction(), 1.000, 3 )
    self.assertAlmostEqual( self.alignment4.identity_fraction(), 1.000, 3 )


  def test_multiplicity(self):

    self.assertEqual( self.alignment1.multiplicity(), 2 )
    self.assertEqual( self.alignment2.multiplicity(), 4 )
    self.assertEqual( self.alignment3.multiplicity(), 0 )
    self.assertEqual( self.alignment4.multiplicity(), 3 )


  def test_length(self):

    self.assertEqual( self.alignment1.length(), 75 )
    self.assertEqual( self.alignment2.length(), 75 )
    self.assertEqual( self.alignment3.length(), 0 )
    self.assertEqual( self.alignment4.length(), 75 )


  def test_copy(self):

    c = self.alignment1.copy(
      alignments = [ "AB", "BA" ],
      names = [ "XY", "YX" ],
      gap = "#"
      )

    self.assert_( c is not self.alignment1 )
    self.assertEqual( c.alignments, [ "AB", "BA" ] )
    self.assertEqual( c.names, [ "XY", "YX" ] )
    self.assertEqual( c.gap, "#" )


  def test__str__(self):
    self.assertEqual(
      str( self.alignment1 ),
      ">1hml\n"
      + "-KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYDTQAIVENN--ESTEYGLFQISNKLWCKSSQ\n"
      + "VPQSR\n\n"
      + ">1hfy\n"
      + "-EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNN--DSTEYGLFQINNKIWCKDDQ\n"
      + "NPHSR"
      )

  def test_extend(self):

    sequences = [
      "ABKQFTKCELSQLLKDIDGYGGIALPELICTMFHTSGYDTQAIVENNESTEYGLFQISNKLWCKSSQVPQSRNOPQR",
      "CDEEQLTKCEVFQKLKDLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNNDSTEYGLFQINNKIWCKDDQNPHSRSTUV",
      "FGHIGKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFNTGATNRNTDGSTDYGILQINSRWWCNDGRTPGSKXYZ",
      "JKILMKVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFDTHATNRNTDGSTDYGILQINSRWWCNDGRTPGSKAB",
      ]
    ali = self.alignment2.copy()
    ali.extend( sequences = sequences )
    self.assertEqual( ali.names, self.alignment2.names )
    self.assertEqual( ali.gap, self.alignment2.gap )
    self.assertEqual(
      ali.alignments,
      [
        "AB------------" + self.alignment2.alignments[0] + "NOPQR---------",
        "--CDE---------" + self.alignment2.alignments[1] + "-----STUV-----",
        "-----FGHI-----" + self.alignment2.alignments[2] + "---------XYZ--",
        "---------JKILM" + self.alignment2.alignments[3] + "------------AB",
        ]
      )


class test_fasta_alignment(unittest.TestCase):

  def setUp(self):

    self.alignment = bioinformatics.fasta_alignment(
      alignments = [ ali_1hml, ali_1hfy, ali_1ghl, ali_1lz3 ],
      names = [ "1hml", "1hfya", "1ghla", "1lz3" ],
      descriptions = [
        "alpha lactalbumin:Homo sapiens",
        "alpha-lactalbumin:Capra hircus",
        "lysozyme:Phasianus colchicus",
        "lysozyme:Meleagris gallopavo",
        ]
      )

  def test_error(self):

    self.assertRaises(
      ValueError,
      bioinformatics.fasta_alignment,
      [ ali_1hml, ali_1hml ],
      [ "1hml", "1hfy" ],
      [
        "alpha lactalbumin:Homo sapiens",
        "alpha-lactalbumin:Capra hircus",
        "lysozyme:Phasianus colchicus",
        ]
      )


  def test_format(self):

    self.assertEqual(
      self.alignment.format( 40 ),
      ">1hml alpha lactalbumin:Homo sapiens\n"
      + "-KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYD\n"
      + "TQAIVENN--ESTEYGLFQISNKLWCKSSQVPQSR\n\n"
      + ">1hfya alpha-lactalbumin:Capra hircus\n"
      + "-EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYD\n"
      + "TQAIVQNN--DSTEYGLFQINNKIWCKDDQNPHSR\n\n"
      + ">1ghla lysozyme:Phasianus colchicus\n"
      + "GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFN\n"
      + "TGATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK\n\n"
      + ">1lz3 lysozyme:Meleagris gallopavo\n"
      + "-KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFD\n"
      + "THATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK"
      )


  def test__str__(self):
    self.assertEqual(
      str( self.alignment ),
      ">1hml alpha lactalbumin:Homo sapiens\n"
      + "-KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYDTQAIVENN--ESTEYGLFQISNKLWCKSSQ\n"
      + "VPQSR\n\n"
      + ">1hfya alpha-lactalbumin:Capra hircus\n"
      + "-EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNN--DSTEYGLFQINNKIWCKDDQ\n"
      + "NPHSR\n\n"
      + ">1ghla lysozyme:Phasianus colchicus\n"
      + "GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFNTGATNRNT-DGSTDYGILQINSRWWCNDGR\n"
      + "TPGSK\n\n"
      + ">1lz3 lysozyme:Meleagris gallopavo\n"
      + "-KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFDTHATNRNT-DGSTDYGILQINSRWWCNDGR\n"
      + "TPGSK"
      )


  def test_copy(self):

    c = self.alignment.copy( descriptions = [ "A", "B", "C", "D" ] )
    self.assertEqual( c.descriptions, [ "A", "B", "C", "D" ] )
    self.assertEqual( c.alignments, self.alignment.alignments )
    self.assertEqual( c.names, self.alignment.names )
    self.assertEqual( c.gap, self.alignment.gap )


class test_pir_alignment(unittest.TestCase):

  def setUp(self):

    self.alignment = bioinformatics.pir_alignment(
      alignments = [ ali_1hml, ali_1hfy, ali_1ghl, ali_1lz3 ],
      names = [ "1hml", "1hfya", "1ghla", "1lz3" ],
      types = [ "P1" ] * 4,
      descriptions = [
        "alpha lactalbumin:Homo sapiens",
        "alpha-lactalbumin:Capra hircus",
        "lysozyme:Phasianus colchicus",
        "lysozyme:Meleagris gallopavo",
        ]
      )

  def test_error(self):

    self.assertRaises(
      ValueError,
      bioinformatics.pir_alignment,
      [ ali_1hml, ali_1hml ],
      [ "1hml", "1hfy" ],
      [ "P1", "P1", "P1" ],
      [
        "alpha lactalbumin:Homo sapiens",
        "alpha-lactalbumin:Capra hircus",
        ]
      )
    self.assertRaises(
      ValueError,
      bioinformatics.pir_alignment,
      [ ali_1hml, ali_1hml ],
      [ "1hml", "1hfy" ],
      [ "P1", "P1" ],
      [
        "alpha lactalbumin:Homo sapiens",
        "alpha-lactalbumin:Capra hircus",
        "lysozyme:Phasianus colchicus",
        ]
      )


  def test_format(self):

    self.assertEqual(
      self.alignment.format( 40 ),
      ">P1;1hml\nalpha lactalbumin:Homo sapiens\n"
      + "  -KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSG\n"
      + "  YDTQAIVENN--ESTEYGLFQISNKLWCKSSQVPQSR*\n\n"
      + ">P1;1hfya\nalpha-lactalbumin:Capra hircus\n"
      + "  -EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSG\n"
      + "  YDTQAIVQNN--DSTEYGLFQINNKIWCKDDQNPHSR*\n\n"
      + ">P1;1ghla\nlysozyme:Phasianus colchicus\n"
      + "  GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESN\n"
      + "  FNTGATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK*\n\n"
      + ">P1;1lz3\nlysozyme:Meleagris gallopavo\n"
      + "  -KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESN\n"
      + "  FDTHATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK*"
      )


  def test__str__(self):
    self.assertEqual(
      str( self.alignment ),
      ">P1;1hml\nalpha lactalbumin:Homo sapiens\n"
      + "  -KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYDTQAIVENN--ESTEYGLFQISNKLWCKS\n"
      + "  SQVPQSR*\n\n"
      + ">P1;1hfya\nalpha-lactalbumin:Capra hircus\n"
      + "  -EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNN--DSTEYGLFQINNKIWCKD\n"
      + "  DQNPHSR*\n\n"
      + ">P1;1ghla\nlysozyme:Phasianus colchicus\n"
      + "  GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFNTGATNRNT-DGSTDYGILQINSRWWCND\n"
      + "  GRTPGSK*\n\n"
      + ">P1;1lz3\nlysozyme:Meleagris gallopavo\n"
      + "  -KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFDTHATNRNT-DGSTDYGILQINSRWWCND\n"
      + "  GRTPGSK*"
      )

  def test_copy(self):

    c = self.alignment.copy(
      descriptions = [ "A", "B", "C", "D" ],
      types = [ "P1", "P2", "P3", "P4" ]
      )
    self.assertEqual( c.descriptions, [ "A", "B", "C", "D" ] )
    self.assertEqual( c.types, [ "P1", "P2", "P3", "P4" ] )
    self.assertEqual( c.alignments, self.alignment.alignments )
    self.assertEqual( c.names, self.alignment.names )
    self.assertEqual( c.gap, self.alignment.gap )


class test_clustal_alignment(unittest.TestCase):

  def setUp(self):

    self.alignment = bioinformatics.clustal_alignment(
      alignments = [ ali_1hml, ali_1hfy, ali_1ghl, ali_1lz3 ],
      names = [ "1hml", "1hfya", "1ghla", "1lz3" ],
      program = "CLUSTAL X 1.35"
      )

  def test_format(self):

    self.assertEqual(
      self.alignment.format( 40, 10 ),
      "CLUSTAL X 1.35 multiple sequence alignment\n\n"
      + "1hml       -KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYD 37\n"
      + "1hfya      -EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYD 37\n"
      + "1ghla      GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFN 40\n"
      + "1lz3       -KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFD 39\n"
      + "                 **     *      * *  *    *     *   \n\n"
      + "1hml       TQAIVENN--ESTEYGLFQISNKLWCKSSQVPQSR 70\n"
      + "1hfya      TQAIVQNN--DSTEYGLFQINNKIWCKDDQNPHSR 70\n"
      + "1ghla      TGATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK 74\n"
      + "1lz3       THATNRNT-DGSTDYGILQINSRWWCNDGRTPGSK 73\n"
      + "           * *   *    ** **  **    **     * * "
      )


  def test__str__(self):

    self.assertEqual(
      str( self.alignment ),
      "CLUSTAL X 1.35 multiple sequence alignment\n\n"
      + "1hml            -KQFTKCELSQLLK--DIDGYGGIALPELICTMFHTSGYDTQAIVENN--ESTEYGLFQI 55\n"
      + "1hfya           -EQLTKCEVFQKLK--DLKDYGGVSLPEWVCTAFHTSGYDTQAIVQNN--DSTEYGLFQI 55\n"
      + "1ghla           GKVYGRCELAAAMKRMGLDNYRGYSLGNWVCAAKFESNFNTGATNRNT-DGSTDYGILQI 59\n"
      + "1lz3            -KVYGRCELAAAMKRLGLDNYRGYSLGNWVCAAKFESNFDTHATNRNT-DGSTDYGILQI 58\n"
      + "                      **     *      * *  *    *     *   * *   *    ** **  **\n\n"
      + "1hml            SNKLWCKSSQVPQSR 70\n"
      + "1hfya           NNKIWCKDDQNPHSR 70\n"
      + "1ghla           NSRWWCNDGRTPGSK 74\n"
      + "1lz3            NSRWWCNDGRTPGSK 73\n"
      + "                    **     * * "
      )

  def test_copy(self):

    c = self.alignment.copy( program = "Foo" )
    self.assertEqual( c.program, "Foo" )
    self.assertEqual( c.alignments, self.alignment.alignments )
    self.assertEqual( c.names, self.alignment.names )
    self.assertEqual( c.gap, self.alignment.gap )


seq = """
Error
>chain_A
 TPDCVTGKVE YTKYNDDDTF TVKVGDKELF TNRWNLQSLL LSAQITGMTV TIKTNACHNG
 GGFSEVIFR
> chain A
 TPDCVTGKVE YTKYNDDDTF TVKVGDKELF TNRWNLQSLL LSAQITGMTV TIKTNACHNG
 GGFSEVIFR

> chain_A
 TPDCVTGKVE YTKYNDDDTF TVKVGDKELF TNRWNLQSLL LSAQITGMTV TIKTNACHNG
 GGFSEVIFR*
>chain_A"""

seq_sequence = "TPDCVTGKVEYTKYNDDDTFTVKVGDKELFTNRWNLQSLLLSAQITGMTVTIKTNACHNGGGFSEVIFR"

fasta = """
Error
>FOSB_MOUSE Protein fosB. 338 bp
     MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
     ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
     GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT
     DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD
     LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY
     TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL

> FOSB_MOUSE Protein fosB. 338 bp
     MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
     ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
     GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT
     DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD
     LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY
     TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL
> 1YJP:A|PDBID|CHAIN|SEQUENCE

     MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA
     ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS
     GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT
     DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD
     LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY
     TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL


"""

fasta_sequence = (
    "MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA"
    + "ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS"
    + "GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT"
    + "DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD"
    + "LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY"
    + "TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL"
    )

pir = """
Error
>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

>P1;CRAB_ANAPL

  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK

> P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

>P11;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

>P1;CRAB_ ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

"""

tolerant_pir = """
Error
>P1;CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

Error2

>CRAB_ANAPL
ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN).
  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK*

>CRAB_ANAPL

  MDITIHNPLI RRPLFSWLAP SRIFDQIFGE HLQESELLPA SPSLSPFLMR
  SPIFRMPSWL ETGLSEMRLE KDKFSVNLDV KHFSPEELKV KVLGDMVEIH
  GKHEERQDEH GFIAREFNRK YRIPADVDPL TITSSLSLDG VLTVSAPRKQ
  SDVPERSIPI TREEKPAIAG AQRK
"""

tf = """
> anb.pdb
AAKDVKFGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTT
TATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKV
GKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLL
IIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINK
DTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAA
VEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEY
GNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLP

AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASK
ANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVG
KLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLE
AVAKAGKPLLIIAEDVEGEALQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVA
QIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKL
ADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQ
YAASVAGLMITTECMVTDLP

AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVAST
ATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVG
KEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLI
IAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKD
TTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAV
EEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGKGGDGNYGYNAATEEYGNMIDMGILDP
TKVTRSALQYAASVAGLMITTECMVTDLP

"""

lineseparated = """
AAKDVKFGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTT
TATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKV
GKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLL
IIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINK
DTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAA
VEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEY
GNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLP

AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASK
ANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVG
KLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLE
AVAKAGKPLLIIAEDVEGEALQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVA
QIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKL
ADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQ
YAASVAGLMITTECMVTDLP

"""

dbfetch = """
>PDB:1GAG_A
VFPSSVFVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESAS
LRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAEN
NPGRPPPTLQEMIQMAAEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDI
XETDXXRKGGKGLLPVRWMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQ
VLKFVMDGGYLDQPDNCPERVTDLMRMCWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFF
HSEENK

>PDB:1GAG_B
PATGDFMNMSPVG
"""

pir_sequence = (
    "MDITIHNPLIRRPLFSWLAPSRIFDQIFGEHLQESELLPASPSLSPFLMR"
    + "SPIFRMPSWLETGLSEMRLEKDKFSVNLDVKHFSPEELKVKVLGDMVEIH"
    + "GKHEERQDEHGFIAREFNRKYRIPADVDPLTITSSLSLDGVLTVSAPRKQ"
    + "SDVPERSIPITREEKPAIAGAQRK"
    )

clustal1 = """CLUSTAL 2.0.10 multiple sequence alignment


Horse           VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKA--- 57
chain_A         VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED 60
                *** .:   *  .*:** .... :* : * *:* ..* *   * :*.  : .*::**

Horse           ---HGKKVGDALTLAVGHLDDLPGALSDLSNLHAHKLRVDPVNFKLLSHCLLSTLAVHLP 114
chain_A         LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHP 120
                   ** .*  **   : : ..  . *. *:: ** * ::    ::::*..:: .*  : *

Horse           NDFTPAVHASLDKFLSSVSTVLTSKYR------ 141
chain_A         GDFGADAQGAMNKALELFRKDIAAKYKELGYQG 153
                .** . .:.:::* *. . . :::**:
"""

clustal2 = """CLUSTAL X (1.81) multiple sequence alignment

1vkk VVCEVDPELKETLRKFRFR---KETNNAAIIMKVD--KDRQMVVLEDELQ-NISPEELKL
1ahq -GIAVSDDCVQKFNELKLGHQH-----RYVTFKMNASN--TEVVVEHVGGPNATYEDFKS

1vkk ELPERQPRFVVYSYKYVH--DDGRVSYPLCFIFSSPVGCKPEQQMMYAGSKNRLVQTAE-
1ahq QLPERDCRYAIFDYEFQVDG---GQRNKITFILWAPDSAPIKSKMMYTSTKDSIKKKLVG

1vkk L-TKVFEIRTTDD-LTETWLKEKLAFFR
1ahq IQ-VEVQATD-AAEISEDAVSERAKK--
"""

pir_ali1 = """
>P1;Horse

VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKA---
---HGKKVGDALTLAVGHLDDLPGALSDLSNLHAHKLRVDPVNFKLLSHCLLSTLAVHLP
NDFTPAVHASLDKFLSSVSTVLTSKYR------*

>P3;chain_A

VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED
LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHP
GDFGADAQGAMNKALELFRKDIAAKYKELGYQG*
"""

fasta_ali1 = """
>2QZU:A Putative sulfatase from B. fragilis
QPTPNLVFIXADQYRGDAIGCIGKEPVKTPHLDKLASEGINFTNAISSYP
VSSPARGXLXTGXYPIGSKVT-GNCNSETAPYGVELSQNARCWSDVLKDQ
GYNXGYIGKWHLDAPYKPYVDTYNNRGKVAWNEWCP---PERRHGFDHWI
AYGTYDY-----------------------------H-LKPXYWNTTAPR
DSFYYVNQ-----------WGPEYEASKAIEYINGQKDQKQPFALVVSXN
PPHTGYE--------------LVPDRYKEIYKDLDVEALCKGRPDIP---
-------------------AKG----------------------TEXGDY
FRNNIRNYYACITGVDENVGRIIEALKQNNLFDNTIVVFTSDHGICXG--
---------------AHENAGKD-IFYEESXRIPXILSWPDQIKPRKSDP
--LXIAFA-DLYPTLLSXXGFSKEIPETVQTFDLSNEVLTGKNKKD-LVQ
PYYFVKF--------DNHA------------------TGYRGLRTDRYTY
AVHAT-DGK------------------------IDNVILFDRTNDPHEXN
NIASQ--QLKLTHTFNRQLKTWLEKTNDPF
>3B5Q:B Putative sulfatase from B. thetaiotaomicron
-EKPNFLIIQCDHLTQRVVGAYGQTQGCTLPIDEVASRGVIFSNAYVGCP
LSQPSRAALWSGXXPHQTNVR-SNSS---EPVNTRLPENVPTLGSLFSES
GYEAVHFGKTHDX---------------------------GSLRGFKHKE
P-------------------------------------------------
---VAKPFTDPEFPVNNDSFLDVGTCEDAVAYLSNP--PKEPFICIADFQ
NPHNICGFIGENAGVHTDRPI----------SGPLPEL----PDNFDVED
WSNIPTPVQYICCSHRRXT---QAA-------------------HWNEEN
YRHYIAAFQHYTKXVSKQVDSVLKALYSTPAGRNTIVVIXADHGDGXA--
---------------SHRXVTKHISFYDEXTNVPFIFAG-PGIKQQKKPV
DHLLTQPTLDLLPTLCDLAGIA--VPAEKAGISLAPTLRGEKQKKSHPYV
VSEWHSEYEYVT-------------------------TPGRXVRGPRYKY
THYLE----------------------------GNGEELYDXKKDPGERK
NLAKDPKYSKILAEHRALLDDYITRSKDDY
"""

ali_ali1 = """
>2QZU:A
QPTPNLVFIXADQYRGDAIGCIGKEPVKTPHLDKLASEGINFTNAISSYP
VSSPARGXLXTGXYPIGSKVT-GNCNSETAPYGVELSQNARCWSDVLKDQ
GYNXGYIGKWHLDAPYKPYVDTYNNRGKVAWNEWCP---PERRHGFDHWI
AYGTYDY-----------------------------H-LKPXYWNTTAPR
DSFYYVNQ-----------WGPEYEASKAIEYINGQKDQKQPFALVVSXN
PPHTGYE--------------LVPDRYKEIYKDLDVEALCKGRPDIP---
-------------------AKG----------------------TEXGDY
FRNNIRNYYACITGVDENVGRIIEALKQNNLFDNTIVVFTSDHGICXG--
---------------AHENAGKD-IFYEESXRIPXILSWPDQIKPRKSDP
--LXIAFA-DLYPTLLSXXGFSKEIPETVQTFDLSNEVLTGKNKKD-LVQ
PYYFVKF--------DNHA------------------TGYRGLRTDRYTY
AVHAT-DGK------------------------IDNVILFDRTNDPHEXN
NIASQ--QLKLTHTFNRQLKTWLEKTNDPF
>3B5Q:B
-EKPNFLIIQCDHLTQRVVGAYGQTQGCTLPIDEVASRGVIFSNAYVGCP
LSQPSRAALWSGXXPHQTNVR-SNSS---EPVNTRLPENVPTLGSLFSES
GYEAVHFGKTHDX---------------------------GSLRGFKHKE
P-------------------------------------------------
---VAKPFTDPEFPVNNDSFLDVGTCEDAVAYLSNP--PKEPFICIADFQ
NPHNICGFIGENAGVHTDRPI----------SGPLPEL----PDNFDVED
WSNIPTPVQYICCSHRRXT---QAA-------------------HWNEEN
YRHYIAAFQHYTKXVSKQVDSVLKALYSTPAGRNTIVVIXADHGDGXA--
---------------SHRXVTKHISFYDEXTNVPFIFAG-PGIKQQKKPV
DHLLTQPTLDLLPTLCDLAGIA--VPAEKAGISLAPTLRGEKQKKSHPYV
VSEWHSEYEYVT-------------------------TPGRXVRGPRYKY
THYLE----------------------------GNGEELYDXKKDPGERK
NLAKDPKYSKILAEHRALLDDYITRSKDDY
"""

class test_sequence_parse(unittest.TestCase):

  def testSequence(self):

    ( sequences, unknowns ) = bioinformatics.seq_sequence_parse( seq )

    self.assertEqual( len( sequences ), 3 )
    self.assertEqual( unknowns, [ "\nError\n", ">chain_A" ] )

    self.assertEqual(
      sequences[0].name,
      "chain_A"
      )
    self.assertEqual(
      sequences[0].sequence,
      seq_sequence
      )

    self.assertEqual(
      sequences[1].name,
      "chain A"
      )
    self.assertEqual(
      sequences[1].sequence,
       seq_sequence
       )

    self.assertEqual(
      sequences[2].name,
      "chain_A"
      )
    self.assertEqual(
      sequences[2].sequence,
      seq_sequence
      )

    ( fastas, unknowns ) = bioinformatics.seq_sequence_parse( fasta )

    self.assertEqual( len( fastas ), 3 )
    self.assertEqual( unknowns, [ "\nError\n" ] )

    self.assertEqual(
      fastas[0].name,
      "FOSB_MOUSE Protein fosB. 338 bp"
      )
    self.assertEqual(
      fastas[0].sequence,
      fasta_sequence
      )

    self.assertEqual(
      fastas[1].name,
      "FOSB_MOUSE Protein fosB. 338 bp"
      )
    self.assertEqual(
      fastas[1].sequence,
      fasta_sequence
      )

    self.assertEqual(
      fastas[2].name,
      "1YJP:A|PDBID|CHAIN|SEQUENCE"
      )
    self.assertEqual(
      fastas[2].sequence,
      fasta_sequence
      )

  def testFasta(self):

    ( fastas, unknowns ) = bioinformatics.fasta_sequence_parse( fasta )

    self.assertEqual( len( fastas ), 3 )
    self.assertEqual(
      unknowns,
      [
        "\nError\n",
        ] )

    self.assertEqual(
      fastas[0].name,
      "FOSB_MOUSE"
      )
    self.assertEqual(
      fastas[0].sequence,
      fasta_sequence
      )
    self.assertEqual(
      fastas[0].description,
      "Protein fosB. 338 bp"
      )

    self.assertEqual(
      fastas[1].name,
      "FOSB_MOUSE"
      )
    self.assertEqual(
      fastas[1].sequence,
      fasta_sequence
      )
    self.assertEqual(
      fastas[1].description,
      "Protein fosB. 338 bp"
      )

    self.assertEqual(
      fastas[2].name,
      "1YJP"
      )
    self.assertEqual(
      fastas[2].description,
      "A|PDBID|CHAIN|SEQUENCE"
      )
    self.assertEqual(
      fastas[2].sequence,
      fasta_sequence
      )


  def testPir(self):

    ( pirs, unknowns ) = bioinformatics.pir_sequence_parse( pir )

    self.assertEqual( len( pirs ), 2 )
    self.assertEqual( unknowns, [ "\nError\n", pir[490:] ] )

    self.assertEqual(
      pirs[0].name,
      "CRAB_ANAPL"
      )
    self.assertEqual(
      pirs[0].type,
      "P1"
      )
    self.assertEqual(
      pirs[0].description,
      "ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN)."
      )
    self.assertEqual(
      pirs[0].sequence,
      pir_sequence
      )

    self.assertEqual(
      pirs[1].name,
      "CRAB_ANAPL"
      )
    self.assertEqual(
      pirs[1].type,
      "P1"
      )
    self.assertEqual(
      pirs[1].description,
      ""
      )
    self.assertEqual(
      pirs[1].sequence,
      pir_sequence
      )


  def testTolerantPir(self):

    ( pirs, unknowns ) = bioinformatics.tolerant_pir_sequence_parse( tolerant_pir )

    self.assertEqual( len( pirs ), 3 )
    self.assertEqual( unknowns, [ "\nError\n", "Error2\n\n" ] )

    self.assertEqual(
      pirs[0].name,
      "P1;CRAB_ANAPL"
      )
    self.assertEqual(
      pirs[0].type,
      "P1"
      )
    self.assertEqual(
      pirs[0].description,
      "ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN)."
      )
    self.assertEqual(
      pirs[0].sequence,
      pir_sequence
      )

    self.assertEqual(
      pirs[1].name,
      "CRAB_ANAPL"
      )
    self.assertEqual(
      pirs[1].type,
      "P1"
      )
    self.assertEqual(
      pirs[1].description,
      "ALPHA CRYSTALLIN B CHAIN (ALPHA(B)-CRYSTALLIN)."
      )
    self.assertEqual(
      pirs[1].sequence,
      pir_sequence
      )
    self.assertEqual(
      pirs[2].name,
      "CRAB_ANAPL"
      )
    self.assertEqual(
      pirs[2].type,
      "P1"
      )
    self.assertEqual(
      pirs[2].description,
      ""
      )
    self.assertEqual(
      pirs[2].sequence,
      pir_sequence
      )


  def test_tf_parse(self):

    ( seqs, unknowns ) = bioinformatics.tf_sequence_parse( text = tf )
    self.assertEqual( unknowns, [] )
    self.assertEqual( len( seqs ), 3 )
    self.assertEqual( seqs[0].name, "anb.pdb" )
    self.assertEqual(
      seqs[0].sequence,
      "AAKDVKFGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTT"
      + "TATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKV"
      + "GKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLL"
      + "IIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINK"
      + "DTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAA"
      + "VEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEY"
      + "GNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLP"
      )
    self.assertEqual( seqs[1].name, "anb.pdb" )
    self.assertEqual(
      seqs[1].sequence,
      "AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASK"
      + "ANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVG"
      + "KLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLE"
      + "AVAKAGKPLLIIAEDVEGEALQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVA"
      + "QIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKL"
      + "ADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQ"
      + "YAASVAGLMITTECMVTDLP"
      )
    self.assertEqual( seqs[2].name, "anb.pdb" )
    self.assertEqual(
      seqs[2].sequence,
      "AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVAST"
      + "ATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVG"
      + "KEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLI"
      + "IAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKD"
      + "TTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAV"
      + "EEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGKGGDGNYGYNAATEEYGNMIDMGILDP"
      + "TKVTRSALQYAASVAGLMITTECMVTDLP"
      )


  def test_ls_parse(self):

    ( seqs, unknowns ) = bioinformatics.lineseparated_sequence_parse.parse(
      text = lineseparated,
      name = "foo"
      )
    self.assertEqual( unknowns, [] )
    self.assertEqual( len( seqs ), 2 )
    self.assertEqual( seqs[0].name, "foo" )
    self.assertEqual(
      seqs[0].sequence,
      "AAKDVKFGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTT"
      + "TATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKV"
      + "GKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLL"
      + "IIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINK"
      + "DTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAA"
      + "VEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEY"
      + "GNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLP"
      )
    self.assertEqual( seqs[1].name, "foo" )
    self.assertEqual(
      seqs[1].sequence,
      "AAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASK"
      + "ANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVG"
      + "KLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLE"
      + "AVAKAGKPLLIIAEDVEGEALQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVA"
      + "QIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKL"
      + "ADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQ"
      + "YAASVAGLMITTECMVTDLP"
      )


  def test_dbfetch_parse(self):

      ( seqs, unknowns ) = bioinformatics.dbfetch_sequence_parse(
        text = dbfetch
        )
      self.assertEqual( unknowns, [] )
      self.assertEqual( len( seqs ), 2 )
      self.assertEqual( seqs[0].name, "1GAG" )
      self.assertEqual( seqs[0].chain, "A" )
      self.assertEqual(
        seqs[0].sequence,
        "VFPSSVFVPDEWEVSREKITLLRELGQGSFGMVYEGNARDIIKGEAETRVAVKTVNESAS"
        + "LRERIEFLNEASVMKGFTCHHVVRLLGVVSKGQPTLVVMELMAHGDLKSYLRSLRPEAEN"
        + "NPGRPPPTLQEMIQMAAEIADGMAYLNAKKFVHRDLAARNCMVAHDFTVKIGDFGMTRDI"
        + "XETDXXRKGGKGLLPVRWMAPESLKDGVFTTSSDMWSFGVVLWEITSLAEQPYQGLSNEQ"
        + "VLKFVMDGGYLDQPDNCPERVTDLMRMCWQFNPKMRPTFLEIVNLLKDDLHPSFPEVSFF"
        + "HSEENK"
        )
      self.assertEqual( seqs[1].name, "1GAG" )
      self.assertEqual( seqs[1].chain, "B" )
      self.assertEqual( seqs[1].sequence, "PATGDFMNMSPVG" )


  def test_filename_selection(self):

    self.assertEqual(
      bioinformatics.sequence_parser_for( "dummy.pir" ),
      bioinformatics.pir_sequence_parse
      )
    self.assertEqual( bioinformatics.sequence_parser_for( "seq" ), None )


  def test_known_formats(self):
    self.assertEqual(
      sorted( bioinformatics.known_sequence_formats() ),
      ['.dat', '.fa', '.faa', '.fasta', '.pir', '.seq']
      )


class test_alignment_parse(unittest.TestCase):

  def testClustal(self):

    ( ali, unknowns ) = bioinformatics.clustal_alignment_parse( clustal1 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.program, "CLUSTAL 2.0.10" )
    self.assertEqual( ali.names, [ "Horse", "chain_A" ] )
    self.assertEqual(
      ali.alignments,
      [
        "VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKA---"
        + "---HGKKVGDALTLAVGHLDDLPGALSDLSNLHAHKLRVDPVNFKLLSHCLLSTLAVHLP"
        + "NDFTPAVHASLDKFLSSVSTVLTSKYR------",
        "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED"
        + "LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHP"
        + "GDFGADAQGAMNKALELFRKDIAAKYKELGYQG"
        ]
      )

    ( ali, unknowns ) = bioinformatics.clustal_alignment_parse( clustal2 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.program, "CLUSTAL X (1.81)" )
    self.assertEqual( ali.names, [ "1vkk", "1ahq" ] )
    self.assertEqual(
      ali.alignments,
      [
        "VVCEVDPELKETLRKFRFR---KETNNAAIIMKVD--KDRQMVVLEDELQ-NISPEELKL"
        + "ELPERQPRFVVYSYKYVH--DDGRVSYPLCFIFSSPVGCKPEQQMMYAGSKNRLVQTAE-"
        + "L-TKVFEIRTTDD-LTETWLKEKLAFFR",
        "-GIAVSDDCVQKFNELKLGHQH-----RYVTFKMNASN--TEVVVEHVGGPNATYEDFKS"
        + "QLPERDCRYAIFDYEFQVDG---GQRNKITFILWAPDSAPIKSKMMYTSTKDSIKKKLVG"
        + "IQ-VEVQATD-AAEISEDAVSERAKK--"
        ]
      )

    ( ali, unknowns ) = bioinformatics.clustal_alignment_parse( "\n" + clustal2 )

    # Bad format
    self.assertEqual( unknowns, "\n" + clustal2 )
    self.assertEqual( ali, None )

    # Empty alignment
    ( ali, unknowns ) = bioinformatics.clustal_alignment_parse(
      "CLUSTAL 2.0.10 multiple sequence alignment\n\n"
      )

    self.assertEqual( unknowns, "" )
    self.assertEqual( ali.program, "CLUSTAL 2.0.10" )
    self.assertEqual( ali.names, [] )
    self.assertEqual( ali.alignments, [] )


  def testPir(self):

    ( ali, unknowns ) = bioinformatics.pir_alignment_parse( pir_ali1 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.names, [ "Horse", "chain_A" ] )
    self.assertEqual( ali.types, [ "P1", "P3" ] )
    self.assertEqual( ali.descriptions, [ "", "" ] )
    self.assertEqual(
      ali.alignments,
      [
        "VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFPTTKTYFPHFDLSHGSAQVKA---"
        + "---HGKKVGDALTLAVGHLDDLPGALSDLSNLHAHKLRVDPVNFKLLSHCLLSTLAVHLP"
        + "NDFTPAVHASLDKFLSSVSTVLTSKYR------",
        "VLSEGEWQLVLHVWAKVEADVAGHGQDILIRLFKSHPETLEKFDRFKHLKTEAEMKASED"
        + "LKKHGVTVLTALGAILKKKGHHEAELKPLAQSHATKHKIPIKYLEFISEAIIHVLHSRHP"
        + "GDFGADAQGAMNKALELFRKDIAAKYKELGYQG"
        ]
      )

    # Bad format
    ( ali, unknowns ) = bioinformatics.pir_alignment_parse( "\n" + clustal2 )

    self.assertEqual( unknowns, "\n" + clustal2 )
    self.assertEqual( ali, None )

    # Empty alignment
    ( ali, unknowns ) = bioinformatics.pir_alignment_parse( "" )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.names, [] )
    self.assertEqual( ali.types, [] )
    self.assertEqual( ali.descriptions, [] )
    self.assertEqual( ali.alignments, [] )


  def testFasta(self):

    ( ali, unknowns ) = bioinformatics.fasta_alignment_parse( fasta_ali1 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.names, [ "2QZU:A", "3B5Q:B" ] )
    self.assertEqual(
      ali.descriptions,
      [
        "Putative sulfatase from B. fragilis",
        "Putative sulfatase from B. thetaiotaomicron",
        ]
      )
    self.assertEqual(
      ali.alignments,
      [
        "QPTPNLVFIXADQYRGDAIGCIGKEPVKTPHLDKLASEGINFTNAISSYP"
        + "VSSPARGXLXTGXYPIGSKVT-GNCNSETAPYGVELSQNARCWSDVLKDQ"
        + "GYNXGYIGKWHLDAPYKPYVDTYNNRGKVAWNEWCP---PERRHGFDHWI"
        + "AYGTYDY-----------------------------H-LKPXYWNTTAPR"
        + "DSFYYVNQ-----------WGPEYEASKAIEYINGQKDQKQPFALVVSXN"
        + "PPHTGYE--------------LVPDRYKEIYKDLDVEALCKGRPDIP---"
        + "-------------------AKG----------------------TEXGDY"
        + "FRNNIRNYYACITGVDENVGRIIEALKQNNLFDNTIVVFTSDHGICXG--"
        + "---------------AHENAGKD-IFYEESXRIPXILSWPDQIKPRKSDP"
        + "--LXIAFA-DLYPTLLSXXGFSKEIPETVQTFDLSNEVLTGKNKKD-LVQ"
        + "PYYFVKF--------DNHA------------------TGYRGLRTDRYTY"
        + "AVHAT-DGK------------------------IDNVILFDRTNDPHEXN"
        + "NIASQ--QLKLTHTFNRQLKTWLEKTNDPF",
        "-EKPNFLIIQCDHLTQRVVGAYGQTQGCTLPIDEVASRGVIFSNAYVGCP"
        + "LSQPSRAALWSGXXPHQTNVR-SNSS---EPVNTRLPENVPTLGSLFSES"
        + "GYEAVHFGKTHDX---------------------------GSLRGFKHKE"
        + "P-------------------------------------------------"
        + "---VAKPFTDPEFPVNNDSFLDVGTCEDAVAYLSNP--PKEPFICIADFQ"
        + "NPHNICGFIGENAGVHTDRPI----------SGPLPEL----PDNFDVED"
        + "WSNIPTPVQYICCSHRRXT---QAA-------------------HWNEEN"
        + "YRHYIAAFQHYTKXVSKQVDSVLKALYSTPAGRNTIVVIXADHGDGXA--"
        + "---------------SHRXVTKHISFYDEXTNVPFIFAG-PGIKQQKKPV"
        + "DHLLTQPTLDLLPTLCDLAGIA--VPAEKAGISLAPTLRGEKQKKSHPYV"
        + "VSEWHSEYEYVT-------------------------TPGRXVRGPRYKY"
        + "THYLE----------------------------GNGEELYDXKKDPGERK"
        + "NLAKDPKYSKILAEHRALLDDYITRSKDDY",
        ]
      )


  def testAli(self):

    ( ali, unknowns ) = bioinformatics.ali_alignment_parse( ali_ali1 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.names, [ "2QZU:A", "3B5Q:B" ] )
    self.assertEqual(
      ali.alignments,
      [
        "QPTPNLVFIXADQYRGDAIGCIGKEPVKTPHLDKLASEGINFTNAISSYP"
        + "VSSPARGXLXTGXYPIGSKVT-GNCNSETAPYGVELSQNARCWSDVLKDQ"
        + "GYNXGYIGKWHLDAPYKPYVDTYNNRGKVAWNEWCP---PERRHGFDHWI"
        + "AYGTYDY-----------------------------H-LKPXYWNTTAPR"
        + "DSFYYVNQ-----------WGPEYEASKAIEYINGQKDQKQPFALVVSXN"
        + "PPHTGYE--------------LVPDRYKEIYKDLDVEALCKGRPDIP---"
        + "-------------------AKG----------------------TEXGDY"
        + "FRNNIRNYYACITGVDENVGRIIEALKQNNLFDNTIVVFTSDHGICXG--"
        + "---------------AHENAGKD-IFYEESXRIPXILSWPDQIKPRKSDP"
        + "--LXIAFA-DLYPTLLSXXGFSKEIPETVQTFDLSNEVLTGKNKKD-LVQ"
        + "PYYFVKF--------DNHA------------------TGYRGLRTDRYTY"
        + "AVHAT-DGK------------------------IDNVILFDRTNDPHEXN"
        + "NIASQ--QLKLTHTFNRQLKTWLEKTNDPF",
        "-EKPNFLIIQCDHLTQRVVGAYGQTQGCTLPIDEVASRGVIFSNAYVGCP"
        + "LSQPSRAALWSGXXPHQTNVR-SNSS---EPVNTRLPENVPTLGSLFSES"
        + "GYEAVHFGKTHDX---------------------------GSLRGFKHKE"
        + "P-------------------------------------------------"
        + "---VAKPFTDPEFPVNNDSFLDVGTCEDAVAYLSNP--PKEPFICIADFQ"
        + "NPHNICGFIGENAGVHTDRPI----------SGPLPEL----PDNFDVED"
        + "WSNIPTPVQYICCSHRRXT---QAA-------------------HWNEEN"
        + "YRHYIAAFQHYTKXVSKQVDSVLKALYSTPAGRNTIVVIXADHGDGXA--"
        + "---------------SHRXVTKHISFYDEXTNVPFIFAG-PGIKQQKKPV"
        + "DHLLTQPTLDLLPTLCDLAGIA--VPAEKAGISLAPTLRGEKQKKSHPYV"
        + "VSEWHSEYEYVT-------------------------TPGRXVRGPRYKY"
        + "THYLE----------------------------GNGEELYDXKKDPGERK"
        + "NLAKDPKYSKILAEHRALLDDYITRSKDDY",
        ]
      )

  def test_filename_selection(self):

    self.assertEqual(
      bioinformatics.alignment_parser_for( "dummy.pir" ),
      bioinformatics.pir_alignment_parse
      )
    self.assertEqual( bioinformatics.alignment_parser_for( "clustal" ), None )


  def test_known_formats(self):

    self.assertEqual(
      sorted( bioinformatics.known_alignment_formats() ),
      [ ".ali", ".aln", ".clustal", ".fa", ".fasta", ".hhr", ".pir" ]
      )


hhsearchout = \
"""
Query         Thu_Sep_09_20:47:04_+0200_2010
Match_columns 205
No_of_seqs    280 out of 1436
Neff          7.1
Searched_HMMs 22773
Date          Thu Sep  9 20:51:32 2010
Command       /cluster/toolkit/production/bioprogs/hhpred/hhsearch -cpu 4 -v 1 -i /cluster/toolkit/production/tmp/production/426773/3338418.hhm -d /cluster/toolkit/production/databases/hhpred/new_dbs/pdb70_9Sep10/db/pdb.hhm -o /cluster/toolkit/production/tmp/production/426773/3338418.hhr -p 20 -P 20 -Z 100 -B 100 -seq 1 -aliw 80 -local -ssm 2 -norealign -sc 1 -dbstrlen 10000 -cs /cluster/toolkit/production/bioprogs/csblast/data/clusters.prf

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 3cng_A Nudix hydrolase; struct 100.0 4.5E-28   2E-32  196.5  15.7  139   59-204    27-172 (189)
  2 3i7u_A AP4A hydrolase; nudix p  99.9 1.4E-27 6.2E-32  182.2  13.1  127   69-196     4-130 (134)

No 1
>3cng_A Nudix hydrolase; structural genomics, APC7497, PSI-2, protein structure initiative, midwest center for structural genomics; 2.00A {Nitrosomonas europaea atcc 19718}
Probab=99.96  E-value=4.5e-28  Score=196.52  Aligned_cols=139  Identities=21%  Similarity=0.317  Sum_probs=0.0

Q ss_pred             hccCCC---CCceEEEEEEEEEECCEEEEEEcCC----CEEECCceeeCCCCCHHHHHHHHHHHHhCCccccceEEEEEe
Q Thu_Sep_09_20:   59 FCNETG---YQTPKLDTRAAIFQEDKILLVQEND----GLWSLPGGWCDVDQSVKDNVVKEVKEEAGLDVEAQRVVAILD  131 (205)
Q Consensus        59 ~~~~~g---y~Tpkv~v~a~v~~d~kiLLv~~~~----g~W~lPGG~ve~gEs~~eaa~REv~EETGl~v~~~~ll~v~~  131 (205)
                       |+.||   |+.|++.|.++|.++++|||++|+.    |.|+||||++|.|||+++||+||++||||+++....+++++.
T Consensus        27 ~C~~C~~~~y~~P~v~v~~ii~~~~~vLLv~r~~~~~~g~W~lPGG~ve~GEs~e~aa~REv~EEtGl~v~~~~l~~~~~  106 (189)
T 3cng_A           27 ICPKCHTIHYQNPKVIVGCIPEWENKVLLCKRAIAPYRGKWTLPAGFMENNETLVQGAARETLEEANARVEIRELYAVYS  106 (189)
T ss_dssp             EETTTTEEECCCCEEEEEEEEEETTEEEEEEESSSSSTTCEECSEEECCTTCCHHHHHHHHHHHHHCCCEEEEEEEEEEE
T ss_pred             eCCCCCCcccCCCceEEEEEeecCceEEEEeccCCCCCCCEeCCcccCcCCCCHHHHHHHHHHhhhceeeeeeEEEEEee


Q ss_pred             eccccCCCCceEEEEEEEEEEecCCccCCCCCeEeEEEEcHHHCccccccCCCHHHHHHHHHHHhCCCCCCCc
Q Thu_Sep_09_20:  132 KHKNNPAKSAHRVTKVFILCRLLGGEFQPNSETVASGFFSLDDLPPLYLGKNTAEQLALCLEASRSEHWETRF  204 (205)
Q Consensus       132 ~~~~~~~~~~~~~~~~~f~~~~~~~~~~~~~E~~e~~Wf~~deLp~Ls~~r~~~~~i~~~f~~~r~~~~~t~f  204 (205)
                      ....       +...++|.|...++.+.++.|+.+++||++++||...+.-....+.-+.|...+..+....+
T Consensus       107 ~~~~-------~~~~~~f~~~~~~~~~~~~~E~~e~~wf~~~elp~~~la~~~~~~~l~~~~~~~~~g~~~~~  172 (189)
T 3cng_A          107 LPHI-------SQVYMLFRAKLLDLDFFPGIESLEVRLFGEQEIPWNDIAFRVIHDPLKRYMEERHHGQPAFH  172 (189)
T ss_dssp             EGGG-------TEEEEEEEEEECCSCCCCCTTEEEEEEECTTTCCGGGBSCHHHHHHHHHHHHHHHHSSCCCE
T ss_pred             cccc-------ceeEEEEEEEeccCcccCcccceeEEEEcHHHCCchhcCcHHHHHHHHHHHHHhhcCCCccc


No 2
>3i7u_A AP4A hydrolase; nudix protein, diadenosine polyphosphate, structural genomics, NPPSFA; HET: PGE; 1.80A {Aquifex aeolicus} PDB: 2pq1_A* 3i7u_A* 3i7v_A*
Probab=99.95  E-value=1.4e-27  Score=182.17  Aligned_cols=127  Identities=20%  Similarity=0.336  Sum_probs=0.0

Q ss_pred             EEEEEEEEEECCEEEEEEcCCCEEECCceeeCCCCCHHHHHHHHHHHHhCCccccceEEEEEeeccccCCCCceEEEEEE
Q Thu_Sep_09_20:   69 KLDTRAAIFQEDKILLVQENDGLWSLPGGWCDVDQSVKDNVVKEVKEEAGLDVEAQRVVAILDKHKNNPAKSAHRVTKVF  148 (205)
Q Consensus        69 kv~v~a~v~~d~kiLLv~~~~g~W~lPGG~ve~gEs~~eaa~REv~EETGl~v~~~~ll~v~~~~~~~~~~~~~~~~~~~  148 (205)
                      ++.++|+|+++|+|||+++.+|.|++|||++|.|||+.+||+||++||||+++....+++..+.. +.......+...++
T Consensus         4 ~~aAg~vv~~~~~vLlv~r~~~~w~~PgG~ve~gEt~~~aa~RE~~EEtGl~~~~~~~~~~~~~~-~~~~~~~~~~~~~~   82 (134)
T 3i7u_A            4 EFSAGGVLFKDGEVLLIKTPSNVWSFPKGNIEPGEKPEETAVREVWEETGVKGEILDYIGEIHYW-YTLKGERIFKTVKY   82 (134)
T ss_dssp             EEEEEEEEEETTEEEEEECTTSCEECCEEECCTTCCHHHHHHHHHHHHHSEEEEEEEEEEEEEEE-EEETTEEEEEEEEE
T ss_pred             EEEEEEEEEECCEEEEEEeCCCcEECceeEeCCCCCHHHHHHhhhhheeceeEEEeeeeeeeeee-ccCCCceEEEEEEE


Q ss_pred             EEEEecCCccCCCCCeEeEEEEcHHHCccccccCCCHHHHHHHHHHHh
Q Thu_Sep_09_20:  149 ILCRLLGGEFQPNSETVASGFFSLDDLPPLYLGKNTAEQLALCLEASR  196 (205)
Q Consensus       149 f~~~~~~~~~~~~~E~~e~~Wf~~deLp~Ls~~r~~~~~i~~~f~~~r  196 (205)
                      |.+...++++.++.|+.+++|+++++++++....+....|..+++..+
T Consensus        83 f~~~~~~~~~~~~~E~~~~~W~~~~e~~~~l~~~~~r~il~~~~~l~~  130 (134)
T 3i7u_A           83 YLMKYKEGEPRPSWEVKDAKFFPIKEAKKLLKYKGDKEIFEKALKLKE  130 (134)
T ss_dssp             EEEEEEEECCCCCTTSSEEEEEEHHHHHHHCCSHHHHHHHHHHHHHHH
T ss_pred             EEEeccCCcccCChhheEEEEEeHHHHHhhcCChHHHHHHHHHHHHHh


Done!
"""

hhalignout = """
Query         1XVQ.pdb
Match_columns 155
No_of_seqs    111 out of 3723
Neff          9.0
Searched_HMMs 0
Date          Tue Jul 13 19:50:02 2010
Command       /work/dimaio/sequence/hhsearch/hhalign -i 1XVQ.hhm -t 1FOH.hhm -glob -o 1FOH_1XVQ.hhr

 No Hit                             Prob E-value P-value  Score    SS Cols Query HMM  Template HMM
  1 1FOH_aln.pdb                     1.0       1       1   56.1  11.5  141    1-155   460-631 (649)

No 1
>1FOH_aln.pdb
Probab=1.00    E-value=1  Score=56.09  Aligned_columns=141  Identities=11%

Q ss_pred             CCCCCCCCCCCCCCEEEEC-CCCCEEEHHH-HC--CCEEEEEEECCCCCCCHHHHHHHHHHHHHHC--------------
Q ss_conf             9988554688178517897-9997764889-69--9848999832437711126669999998420--------------
Q 1XVQ.pdb          1 TVGELPAVGSPAPAFTLTG-GDLGVISSDQ-FR--GKSVLLNIFPSVDTPVCATSVRTFDERAAAS--------------   62 (155)
Q Consensus         1 ~~~~~p~vG~~aPdF~l~~-~~G~~v~ls~-~~--Gk~vvl~f~~~~~~p~c~~~~~~l~~~~~~~--------------   62 (155)
                      -...-+.+|..+|+.-+.. .+|....+.+ +.  |++.++.|......+.-...+..+.+...................
T Consensus       460 ~~~~~~~~G~r~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  539 (649)
T 1FOH_aln.pdb    460 ELAKNCVVGTRFKSQPVVRHSEGLWMHFGDRLVTDGRFRIIVFAGKATDATQMSRIKKFSAYLDSENSVISLYTPKVSDR  539 (649)
T ss_pred             CCCCCCCCCCCCCCCCEEECCCCCCCCHHHHHCCCCCEEEEEECCCCCCCHHHHHHHHHHHHHCCCCCEEEECCCCCCCC
T ss_conf             33688888864688503532788831125541278837999943888761135678888875204563011022466766


Q ss_pred             --CEEECCCCCCCHH-----HHHHHHHHHCCCC------EEECCCCCHHHHHHCCCCCCCCCCCCCEEEEEEEECCCCEE
Q ss_conf             --1000023467889-----9999999828975------11014433488997077344565557101279998699839
Q 1XVQ.pdb         63 --GATVLCVSKDLPF-----AQKRFCGAEGTEN------VMPASAFRDSFGEDYGVTIADGPMAGLLARAIVVIGADGNV  129 (155)
Q Consensus        63 --g~~v~~i~~~~~~-----~~~~~~~~~~~~~------~~~~~d~~~~~~~~~g~~~~~~~~~~~~~p~~fiID~~G~I  129 (155)
                        .+.++.|......     .............      .....|........||+..+        ..+.+||=|||-|
T Consensus       540 ~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~g~~~~--------~~~~vlVRPD~~V  611 (649)
T 1FOH_aln.pdb    540 NSRIDVITIHSCHRDDIEMHDFPAPALHPKWQYDFIYADCDSWHHPHPKSYQAWGVDET--------KGAVVVVRPDGYT  611 (649)
T ss_pred             CCCEEEEEEECCCCCCCCHHHCHHHHCCCCCCCEEEEECCCCCCCCCHHHHHHHCCCHH--------CCEEEEECCCCEE
T ss_conf             56168999854776544233320200246742016650444444667789998098711--------4379998699558


Q ss_pred             EEEEECCCCCCCCCHHHHHHHHHHHC
Q ss_conf             99998178888889999999998519
Q 1XVQ.pdb        130 AYTELVPEIAQEPNYEAALAALGATS  155 (155)
Q Consensus       130 ~~~~~~~~~~~~~~~~~il~~lkal~  155 (155)
                      .++.  +    .-+.+++.+.+...=
T Consensus       612 a~~~--~----~~~~~~l~~~~~~~l  631 (649)
T 1FOH_aln.pdb    612 SLVT--D----LEGTAEIDRYFSGIL  631 (649)
T ss_pred             EEEE--C----HHHHHHHHHHHHHHC
T ss_conf             9997--3----013589999999742


Done!
"""


class test_hhsearch_parser(unittest.TestCase):

  def setUp(self):

    self.hss = bioinformatics.hhsearch_parser( output = hhsearchout )


  def test_process_header(self):

    self.assertEqual( self.hss.query, "Thu_Sep_09_20:47:04_+0200_2010" )
    self.assertEqual( self.hss.match_columns, 205 )
    self.assertEqual( self.hss.no_of_sequences, ( 280, 1436 ) )
    self.assertEqual( self.hss.neff, "7.1" )
    self.assertEqual( self.hss.searched_hmms, 22773 )
    self.assertEqual( self.hss.date, "Thu Sep  9 20:51:32 2010" )
    self.assertEqual(
      self.hss.command,
      "/cluster/toolkit/production/bioprogs/hhpred/hhsearch -cpu 4 -v 1 -i /cluster/toolkit/production/tmp/production/426773/3338418.hhm -d /cluster/toolkit/production/databases/hhpred/new_dbs/pdb70_9Sep10/db/pdb.hhm -o /cluster/toolkit/production/tmp/production/426773/3338418.hhr -p 20 -P 20 -Z 100 -B 100 -seq 1 -aliw 80 -local -ssm 2 -norealign -sc 1 -dbstrlen 10000 -cs /cluster/toolkit/production/bioprogs/csblast/data/clusters.prf"
      )


  def test_process_hits(self):

    self.assertEqual( self.hss.indices, [ 1, 2 ] )
    self.assertEqual( self.hss.pdbs, [ "3CNG", "3I7U" ] )
    self.assertEqual( self.hss.chains, [ "A", "A" ] )
    self.assertEqual(
      self.hss.annotations,
      [
        "Nudix hydrolase; structural genomics, APC7497, PSI-2, protein structure initiative, midwest center for structural genomics; 2.00A {Nitrosomonas europaea atcc 19718}",
        "AP4A hydrolase; nudix protein, diadenosine polyphosphate, structural genomics, NPPSFA; HET: PGE; 1.80A {Aquifex aeolicus} PDB: 2pq1_A* 3i7u_A* 3i7v_A*",
        ]
        )

    self.assertIterablesAlmostEqual( self.hss.probabs, [ 99.96, 99.95 ], 2 )
    self.assertIterablesAlmostEqual( self.hss.e_values, [ 4.5e-28, 1.4e-27 ], 28 )
    self.assertIterablesAlmostEqual( self.hss.scores, [ 196.52, 182.17 ], 2 )
    self.assertEqual( self.hss.aligned_cols, [ 139, 127 ] )
    self.assertIterablesAlmostEqual( self.hss.identities, [ 21, 20 ], 1 )
    self.assertIterablesAlmostEqual( self.hss.similarities, [ 0.317, 0.336 ], 3 )
    self.assertIterablesAlmostEqual( self.hss.sum_probs, [ 0, 0 ], 1 )
    self.assertEqual( self.hss.query_starts, [ 59, 69 ] )
    self.assertEqual( self.hss.query_ends, [ 204, 196 ] )
    self.assertEqual( self.hss.query_others, [ 205, 205 ] )
    self.assertEqual(
      self.hss.query_alignments,
      [
        "FCNETG---YQTPKLDTRAAIFQEDKILLVQEND----GLWSLPGGWCDVDQSVKDNVVKEVKEEAGLDVEAQRVVAILD"
          + "KHKNNPAKSAHRVTKVFILCRLLGGEFQPNSETVASGFFSLDDLPPLYLGKNTAEQLALCLEASRSEHWETRF",
        "KLDTRAAIFQEDKILLVQENDGLWSLPGGWCDVDQSVKDNVVKEVKEEAGLDVEAQRVVAILDKHKNNPAKSAHRVTKVF"
          + "ILCRLLGGEFQPNSETVASGFFSLDDLPPLYLGKNTAEQLALCLEASR",
        ]
      )
    self.assertEqual(
      self.hss.query_consensi,
      [
        "~~~~~g---y~Tpkv~v~a~v~~d~kiLLv~~~~----g~W~lPGG~ve~gEs~~eaa~REv~EETGl~v~~~~ll~v~~"
          + "~~~~~~~~~~~~~~~~~f~~~~~~~~~~~~~E~~e~~Wf~~deLp~Ls~~r~~~~~i~~~f~~~r~~~~~t~f",
        "kv~v~a~v~~d~kiLLv~~~~g~W~lPGG~ve~gEs~~eaa~REv~EETGl~v~~~~ll~v~~~~~~~~~~~~~~~~~~~"
          + "f~~~~~~~~~~~~~E~~e~~Wf~~deLp~Ls~~r~~~~~i~~~f~~~r",
        ]
      )
    self.assertEqual(
      self.hss.query_ss_preds,
      [
        "hccCCC---CCceEEEEEEEEEECCEEEEEEcCC----CEEECCceeeCCCCCHHHHHHHHHHHHhCCccccceEEEEEe"
          + "eccccCCCCceEEEEEEEEEEecCCccCCCCCeEeEEEEcHHHCccccccCCCHHHHHHHHHHHhCCCCCCCc",
        "EEEEEEEEEECCEEEEEEcCCCEEECCceeeCCCCCHHHHHHHHHHHHhCCccccceEEEEEeeccccCCCCceEEEEEE"
          + "EEEEecCCccCCCCCeEeEEEEcHHHCccccccCCCHHHHHHHHHHHh",
        ]
      )
    self.assertEqual(
      self.hss.midlines,
      [
        " |+.||   |+.|++.|.++|.++++|||++|+.    |.|+||||++|.|||+++||+||++||||+++....+++++."
          + "....       +...++|.|...++.+.++.|+.+++||++++||...+.-....+.-+.|...+..+....+",
        "++.++|+|+++|+|||+++.+|.|++|||++|.|||+.+||+||++||||+++....+++..+.. +.......+...++"
          + "|.+...++++.++.|+.+++|+++++++++....+....|..+++..+",
        ]
      )

    self.assertEqual( self.hss.hit_starts, [ 27, 4 ] )
    self.assertEqual( self.hss.hit_ends, [ 172, 130 ] )
    self.assertEqual( self.hss.hit_others, [ 189, 134 ] )
    self.assertEqual(
      self.hss.hit_ss_preds,
      [
        "eCCCCCCcccCCCceEEEEEeecCceEEEEeccCCCCCCCEeCCcccCcCCCCHHHHHHHHHHhhhceeeeeeEEEEEee"
          + "cccc-------ceeEEEEEEEeccCcccCcccceeEEEEcHHHCCchhcCcHHHHHHHHHHHHHhhcCCCccc",
        "EEEEEEEEEECCEEEEEEeCCCcEECceeEeCCCCCHHHHHHhhhhheeceeEEEeeeeeeeeee-ccCCCceEEEEEEE"
          + "EEEeccCCcccCChhheEEEEEeHHHHHhhcCChHHHHHHHHHHHHHh",
        ]
      )
    self.assertEqual(
      self.hss.hit_ss_dssps,
      [
        "EETTTTEEECCCCEEEEEEEEEETTEEEEEEESSSSSTTCEECSEEECCTTCCHHHHHHHHHHHHHCCCEEEEEEEEEEE"
          + "EGGG-------TEEEEEEEEEECCSCCCCCTTEEEEEEECTTTCCGGGBSCHHHHHHHHHHHHHHHHSSCCCE",
        "EEEEEEEEEETTEEEEEECTTSCEECCEEECCTTCCHHHHHHHHHHHHHSEEEEEEEEEEEEEEE-EEETTEEEEEEEEE"
          + "EEEEEEEECCCCCTTSSEEEEEEHHHHHHHCCSHHHHHHHHHHHHHHH",
        ]
      )
    self.assertEqual(
      self.hss.hit_consensi,
      [
        "~C~~C~~~~y~~P~v~v~~ii~~~~~vLLv~r~~~~~~g~W~lPGG~ve~GEs~e~aa~REv~EEtGl~v~~~~l~~~~~"
          + "~~~~-------~~~~~~f~~~~~~~~~~~~~E~~e~~wf~~~elp~~~la~~~~~~~l~~~~~~~~~g~~~~~",
        "~~aAg~vv~~~~~vLlv~r~~~~w~~PgG~ve~gEt~~~aa~RE~~EEtGl~~~~~~~~~~~~~~-~~~~~~~~~~~~~~"
          + "f~~~~~~~~~~~~~E~~~~~W~~~~e~~~~l~~~~~r~il~~~~~l~~",
        ]
      )
    self.assertEqual(
      self.hss.hit_alignments,
      [
        "ICPKCHTIHYQNPKVIVGCIPEWENKVLLCKRAIAPYRGKWTLPAGFMENNETLVQGAARETLEEANARVEIRELYAVYS"
          + "LPHI-------SQVYMLFRAKLLDLDFFPGIESLEVRLFGEQEIPWNDIAFRVIHDPLKRYMEERHHGQPAFH",
        "EFSAGGVLFKDGEVLLIKTPSNVWSFPKGNIEPGEKPEETAVREVWEETGVKGEILDYIGEIHYW-YTLKGERIFKTVKY"
          + "YLMKYKEGEPRPSWEVKDAKFFPIKEAKKLLKYKGDKEIFEKALKLKE",
        ]
      )


  def assertIterablesAlmostEqual(self, i1, i2, digits):

    self.assertEqual( len( i1 ), len( i2 ) )

    for ( o, e ) in zip( i1, i2 ):
      self.assertAlmostEqual( o, e, digits )


class test_hhalign_parser(unittest.TestCase):

  def setUp(self):

    self.hss = bioinformatics.hhalign_parser( output = hhalignout )


  def test_process_header(self):

    self.assertEqual( self.hss.query, "1XVQ.pdb" )
    self.assertEqual( self.hss.match_columns, 155 )
    self.assertEqual( self.hss.no_of_sequences, ( 111, 3723 ) )
    self.assertEqual( self.hss.neff, "9.0" )
    self.assertEqual( self.hss.searched_hmms, 0 )
    self.assertEqual( self.hss.date, "Tue Jul 13 19:50:02 2010" )
    self.assertEqual(
      self.hss.command,
      "/work/dimaio/sequence/hhsearch/hhalign -i 1XVQ.hhm -t 1FOH.hhm -glob -o 1FOH_1XVQ.hhr"
      )


  def test_process_hits(self):

    self.assertEqual( self.hss.indices, [ 1 ] )
    self.assertEqual( self.hss.annotations, [ "1FOH_aln.pdb" ] )

    self.assertIterablesAlmostEqual( self.hss.probabs, [ 1.0 ], 2 )
    self.assertIterablesAlmostEqual( self.hss.e_values, [ 1.0 ], 2 )
    self.assertIterablesAlmostEqual( self.hss.scores, [ 56.09 ], 2 )
    self.assertEqual( self.hss.aligned_cols, [ 141 ] )
    self.assertIterablesAlmostEqual( self.hss.identities, [ 11 ], 1 )
    self.assertEqual( self.hss.query_starts, [ 1 ] )
    self.assertEqual( self.hss.query_ends, [ 155 ] )
    self.assertEqual( self.hss.query_others, [ 155 ] )
    self.assertEqual(
      self.hss.query_alignments,
      [
        "TVGELPAVGSPAPAFTLTG-GDLGVISSDQ-FR--GKSVLLNIFPSVDTPVCATSVRTFDERAAAS--------------"
          + "--GATVLCVSKDLPF-----AQKRFCGAEGTEN------VMPASAFRDSFGEDYGVTIADGPMAGLLARAIVVIGADGNV"
          + "AYTELVPEIAQEPNYEAALAALGATS",
        ]
      )
    self.assertEqual(
      self.hss.query_consensi,
      [
        "~~~~~p~vG~~aPdF~l~~-~~G~~v~ls~-~~--Gk~vvl~f~~~~~~p~c~~~~~~l~~~~~~~--------------"
          + "--g~~v~~i~~~~~~-----~~~~~~~~~~~~~------~~~~~d~~~~~~~~~g~~~~~~~~~~~~~p~~fiID~~G~I"
          + "~~~~~~~~~~~~~~~~~il~~lkal~",
        ]
      )
    self.assertEqual(
      self.hss.query_ss_preds,
      [
        "CCCCCCCCCCCCCCEEEEC-CCCCEEEHHH-HC--CCEEEEEEECCCCCCCHHHHHHHHHHHHHHC--------------"
          + "--CEEECCCCCCCHH-----HHHHHHHHHCCCC------EEECCCCCHHHHHHCCCCCCCCCCCCCEEEEEEEECCCCEE"
          + "EEEEECCCCCCCCCHHHHHHHHHHHC",
        ]
      )
    self.assertEqual(
      self.hss.query_ss_confs,
      [
        "9988554688178517897-9997764889-69--9848999832437711126669999998420--------------"
          + "--1000023467889-----9999999828975------11014433488997077344565557101279998699839"
          + "99998178888889999999998519",
        ]
      )

    self.assertEqual(
      self.hss.midlines,
      [
        "-...-+.+|..+|+.-+.. .+|....+.+ +.  |++.++.|......+.-...+..+.+..................."
          + "  .+.++.|......     .............      .....|........||+..+        ..+.+||=|||-|"
          + ".++.  +    .-+.+++.+.+...=",
        ]
      )

    self.assertEqual( self.hss.hit_starts, [ 460 ] )
    self.assertEqual( self.hss.hit_ends, [ 631 ] )
    self.assertEqual( self.hss.hit_others, [ 649 ] )
    self.assertEqual(
      self.hss.hit_ss_preds,
      [
        "CCCCCCCCCCCCCCCCEEECCCCCCCCHHHHHCCCCCEEEEEECCCCCCCHHHHHHHHHHHHHCCCCCEEEECCCCCCCC"
          + "CCCEEEEEEECCCCCCCCHHHCHHHHCCCCCCCEEEEECCCCCCCCCHHHHHHHCCCHH--------CCEEEEECCCCEE"
          + "EEEE--C----HHHHHHHHHHHHHHC",
        ]
      )
    self.assertEqual(
      self.hss.hit_ss_confs,
      [
        "33688888864688503532788831125541278837999943888761135678888875204563011022466766"
          + "56168999854776544233320200246742016650444444667789998098711--------4379998699558"
          + "9997--3----013589999999742",
        ]
      )
    self.assertEqual(
      self.hss.hit_consensi,
      [
        "~~~~~~~~G~r~p~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
          + "~~~~~~~~v~~~~~~~~~~~~~~~~~~~~~~~~~~v~~~~~~~~~~~~~~~~~~g~~~~--------~~~~vlVRPD~~V"
          + "a~~~--~----~~~~~~l~~~~~~~l",
        ]
      )
    self.assertEqual(
      self.hss.hit_alignments,
      [
        "ELAKNCVVGTRFKSQPVVRHSEGLWMHFGDRLVTDGRFRIIVFAGKATDATQMSRIKKFSAYLDSENSVISLYTPKVSDR"
          + "NSRIDVITIHSCHRDDIEMHDFPAPALHPKWQYDFIYADCDSWHHPHPKSYQAWGVDET--------KGAVVVVRPDGYT"
          + "SLVT--D----LEGTAEIDRYFSGIL",
        ]
      )


  def assertIterablesAlmostEqual(self, i1, i2, digits):

    self.assertEqual( len( i1 ), len( i2 ) )

    for ( o, e ) in zip( i1, i2 ):
      self.assertAlmostEqual( o, e, digits )


suite_sequence = unittest.TestLoader().loadTestsFromTestCase(
  test_sequence
  )
suite_fasta_sequence = unittest.TestLoader().loadTestsFromTestCase(
  test_fasta_sequence
  )
suite_pir_sequence = unittest.TestLoader().loadTestsFromTestCase(
  test_pir_sequence
  )
suite_midline = unittest.TestLoader().loadTestsFromTestCase(
  test_midline
  )
suite_alignment = unittest.TestLoader().loadTestsFromTestCase(
  test_alignment
  )
suite_fasta_alignment = unittest.TestLoader().loadTestsFromTestCase(
  test_fasta_alignment
  )
suite_pir_alignment = unittest.TestLoader().loadTestsFromTestCase(
  test_pir_alignment
  )
suite_clustal_alignment = unittest.TestLoader().loadTestsFromTestCase(
  test_clustal_alignment
  )
suite_sequence_parse = unittest.TestLoader().loadTestsFromTestCase(
    test_sequence_parse
    )
suite_alignment_parse = unittest.TestLoader().loadTestsFromTestCase(
    test_alignment_parse
    )
suite_hhsearch_parser = unittest.TestLoader().loadTestsFromTestCase(
    test_hhsearch_parser
    )
suite_hhalign_parser = unittest.TestLoader().loadTestsFromTestCase(
    test_hhalign_parser
    )

alltests = unittest.TestSuite(
  [
    suite_sequence,
    suite_fasta_sequence,
    suite_pir_sequence,
    suite_midline,
    suite_alignment,
    suite_fasta_alignment,
    suite_pir_alignment,
    suite_clustal_alignment,
    suite_sequence_parse,
    suite_alignment_parse,
    suite_hhsearch_parser,
    suite_hhalign_parser,
    ]
  )

if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

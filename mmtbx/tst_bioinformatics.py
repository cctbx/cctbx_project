from mmtbx import bioinformatics

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
    

  def testError(self):

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
    

  def testIdentityCount(self):

    self.assertEqual( self.alignment1.identity_count(), 49 )
    self.assertEqual( self.alignment2.identity_count(), 21 )
    self.assertEqual( self.alignment3.identity_count(), 0 )
    self.assertEqual( self.alignment4.identity_count(), 0 )
    
    
  def testIdentityFraction(self):
      
    self.assertAlmostEqual( self.alignment1.identity_fraction(), 0.700, 3 )
    self.assertAlmostEqual( self.alignment2.identity_fraction(), 0.300, 3 )
    self.assertAlmostEqual( self.alignment3.identity_fraction(), 1.000, 3 )
    self.assertAlmostEqual( self.alignment4.identity_fraction(), 1.000, 3 )


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

  def testError(self):

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


  def testFormat(self):

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


  def testStr(self):
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

  def testError(self):

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


  def testFormat(self):

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


  def testStr(self):
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


class test_clustal_alignment(unittest.TestCase):

  def setUp(self):

    self.alignment = bioinformatics.clustal_alignment(
      alignments = [ ali_1hml, ali_1hfy, ali_1ghl, ali_1lz3 ],
      names = [ "1hml", "1hfya", "1ghla", "1lz3" ],
      program = "X",
      version = "1.35"
      )

  def testFormat(self):

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


  def testStr(self):

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

    self.assertEqual( len( fastas ), 2 )
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


  def testFasta(self):

    ( fastas, unknowns ) = bioinformatics.fasta_sequence_parse( fasta )

    self.assertEqual( len( fastas ), 1 )
    self.assertEqual(
      unknowns,
      [
        "\nError\n",
        "> FOSB_MOUSE Protein fosB. 338 bp\n"
        + "     MFQAFPGDYDSGSRCSSSPSAESQYLSSVDSFGSPPTAAASQECAGLGEMPGSFVPTVTA\n"
        + "     ITTSQDLQWLVQPTLISSMAQSQGQPLASQPPAVDPYDMPGTSYSTPGLSAYSTGGASGS\n"
        + "     GGPSTSTTTSGPVSARPARARPRRPREETLTPEEEEKRRVRRERNKLAAAKCRNRRRELT\n"
        + "     DRLQAETDQLEEEKAELESEIAELQKEKERLEFVLVAHKPGCKIPYEEGPGPGPLAEVRD\n"
        + "     LPGSTSAKEDGFGWLLPPPPPPPLPFQSSRDAPPNLTASLFTHSEVQVLGDPFPVVSPSY\n"
        + "     TSSFVLTCPEVSAFAGAQRTSGSEQPSDPLNSPSLLAL\n\n"
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


  def test_filename_selection(self):

    self.assertEqual(
      bioinformatics.sequence_parser_for( "dummy.pir" ),
      bioinformatics.pir_sequence_parse
      )
    self.assertEqual( bioinformatics.sequence_parser_for( "seq" ), None )


  def test_known_formats(self):

    self.assertEqual(
      sorted( bioinformatics.known_sequence_formats() ),
      [ ".fasta", ".pir", ".seq" ]
      )


class test_alignment_parse(unittest.TestCase):

  def testClustal(self):

    ( ali, unknowns ) = bioinformatics.clustal_alignment_parse( clustal1 )

    self.assertEqual( unknowns, "" )

    self.assertEqual( ali.program, "" )
    self.assertEqual( ali.version, "2.0.10" )
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

    self.assertEqual( ali.program, "X" )
    self.assertEqual( ali.version, "1.81" )
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
    self.assertEqual( ali.program, "" )
    self.assertEqual( ali.version, "2.0.10" )
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
      [ ".ali", ".aln", ".clustal", ".fasta", ".pir" ]
      )


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
    ]
  )

if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

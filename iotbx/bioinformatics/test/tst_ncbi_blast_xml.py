from __future__ import absolute_import, division, print_function

from iotbx.bioinformatics import ncbi_blast_xml

import unittest


class TestParse(unittest.TestCase):

  HSEARCH = None

  def setUp(self):

    if self.__class__.HSEARCH is None:
      import libtbx.load_env
      root = libtbx.env.under_dist( "iotbx", "bioinformatics" )

      import os.path
      import gzip
      self.__class__.HSEARCH = ncbi_blast_xml.parse(
        source = gzip.open( os.path.join( root, "test", "ncbi_blast.xml.gz" ) )
        )


  def testHsps(self):

    hsp = self.HSEARCH.root.iterations[0].hits[4].hsps[0]

    self.assertEqual( hsp.num, 1 )
    self.assertAlmostEqual( hsp.bit_score, 159.073, 3 )
    self.assertEqual( hsp.score, 401 )
    self.assertAlmostEqual( hsp.evalue, 6.95488e-40, 8 )
    self.assertEqual( hsp.query.start, 8 )
    self.assertEqual( hsp.query.end, 192 )
    self.assertEqual( hsp.query.frame, 0 )
    self.assertEqual( hsp.hit.start, 2 )
    self.assertEqual( hsp.hit.end, 193 )
    self.assertEqual( hsp.hit.frame, 0 )
    self.assertEqual( hsp.identity, 78 )
    self.assertEqual( hsp.positive, 127 )
    self.assertEqual( hsp.gaps, 7 )
    self.assertEqual( hsp.length, 192 )
    self.assertEqual(
      hsp.query.seq,
      "KSKIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRAEVSS-GSARGKMLSEIMEKGQLVPLET"
      + "VLDMLR---DAMVAKVDTSKGFLIDGYPREVKQGEEFERKI---GQPTLLLYVDAGPETMTKRL"
      + "LKRGETSGRVDDNEETIKKRLETYYKATEPVIAFYEKRGIVRKVNAEGSVDDVFSQVCTHLD"
      )
    self.assertEqual(
      hsp.hit.seq,
      "KPLVVFVLGGPGAGKGTQCARIVEKYGYTHLSAGELLRDERKNPDSQYGELIEKYIKEGKIVPVEI"
      + "TISLLKREMDQTMAANAQKNKFLIDGFPRNQDNLQGWNKTMDGKADVSFVLFFDCNNEICIERC"
      + "LERGKSSGRSDDNRESLEKRIQTYLQSTKPIIDLYEEMGKVKKIDASKSVDEVFDEVVQIFD"
      )
    self.assertEqual(
      hsp.midline,
      "K  ++FV+GGPG+GKGTQC +IV+KYGYTHLS G+LLR E  +  S  G+++ + +++G++VP+E "
      + " + +L+   D  +A       FLIDG+PR     + + + +      + +L+ D   E   +R "
      + "L+RG++SGR DDN E+++KR++TY ++T+P+I  YE+ G V+K++A  SVD+VF +V    D"
      )


  def testHits(self):

    hit = self.HSEARCH.root.iterations[0].hits[4]

    self.assertEqual( hit.num, 5 )
    self.assertEqual( hit.identifier, "gi|50513735|pdb|1TEV|A" )
    self.assertEqual(
      hit.annotation,
      "Chain A, Crystal Structure Of The Human UmpCMP KINASE IN OPEN Conformation"
      )
    self.assertEqual( hit.accession, "1TEV_A" )
    self.assertEqual( hit.length, 196 )
    self.assertEqual( len( hit.hsps ), 1 )


  def testIterations(self):

    iteration = self.HSEARCH.root.iterations[0]

    self.assertEqual( iteration.num, 1 )
    self.assertEqual( iteration.query_id, "77636" )
    self.assertEqual( iteration.query_def, "3ADK:A|PDBID|CHAIN|SEQUENCE" )
    self.assertEqual( iteration.query_len, 195 )
    self.assertEqual( len( iteration.hits ), 68 )
    self.assertEqual( iteration.statistics.db_num, 40680 )
    self.assertEqual( iteration.statistics.db_len, 9289891 )
    self.assertEqual( iteration.statistics.hsp_len, 0 )
    self.assertEqual( iteration.statistics.eff_space, 0 )
    self.assertEqual( iteration.statistics.kappa, 0.041 )
    self.assertEqual( iteration.statistics.lambdav, 0.267 )
    self.assertEqual( iteration.statistics.entropy, 0.14 )


  def testRoot(self):

    output = self.HSEARCH.root

    self.assertEqual( output.program, "blastp" )
    self.assertEqual( output.version, "BLASTP 2.2.20+" )
    self.assertEqual(
      output.reference,
      'Alejandro A. Sch&auml;ffer, L. Aravind, Thomas L. Madden, '
      + 'Sergei Shavirin, John L. Spouge, Yuri I. Wolf, Eugene V. Koonin, '
      + 'and Stephen F. Altschul (2001), '
      +'"Improving the accuracy of PSI-BLAST protein database searches with '
      + 'composition-based statistics and other refinements", '
      + 'Nucleic Acids Res. 29:2994-3005.'
      )
    self.assertEqual( output.db, "pdb" )
    self.assertEqual( output.query_id, "77636" )
    self.assertEqual( output.query_def, "3ADK:A|PDBID|CHAIN|SEQUENCE" )
    self.assertEqual( output.query_len, 195 )
    self.assertEqual( len( output.iterations ), 1 )
    self.assertEqual( output.parameters.matrix, "BLOSUM62" )
    self.assertEqual( output.parameters.expect, 10 )
    self.assertEqual( output.parameters.gap_open, 11 )
    self.assertEqual( output.parameters.gap_extend, 1 )
    self.assertEqual( output.parameters.filter, "F" )


suite_parse = unittest.TestLoader().loadTestsFromTestCase( TestParse )

alltests = unittest.TestSuite(
  [
    suite_parse,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )


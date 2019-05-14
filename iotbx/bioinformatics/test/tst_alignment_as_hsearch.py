from __future__ import absolute_import, division, print_function

from iotbx.bioinformatics import alignment_as_hsearch

import unittest


class TestParse(unittest.TestCase):

  HSEARCH = None

  def setUp(self):

    if self.__class__.HSEARCH is None:
      import libtbx.load_env
      root = libtbx.env.under_dist( "iotbx", "bioinformatics" )

      import os.path
      self.__class__.HSEARCH = alignment_as_hsearch.parse(
        source = open( os.path.join( root, "test", "alignment.pir" ) )
        )

  def test_len(self):

    self.assertEqual( len( self.HSEARCH ), 2 )


  def test_hits(self):

    hits = list( self.HSEARCH.hits() )

    self.assertEqual( len( hits ), 2 )

    self.assertEqual( hits[0].identifier, "3ADK" )
    self.assertEqual( hits[0].chain, "A" )
    self.assertEqual( hits[0].annotation, "First hit" )
    self.assertEqual( hits[0].alignment.multiplicity(), 2 )
    self.assertEqual( hits[0].alignment.names[0], "target" )
    self.assertEqual( hits[0].alignment.names[1], "3ADK_A" )
    self.assertEqual(
      hits[0].alignment.alignments[0],
      "MEEKLKKSKIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRAEVSSGSARGKMLSEIMEKGQLVPLETVLDMLRDAM"
      + "VAKVDTSKGFLIDGYPREVKQGEEFERKIGQPTLLLYVDAGPETMTKRLLKRGETSGRVDDNEETIKKRLETYYKATEPV"
      + "IAFYEKRGIVRKVNAEGSVDDVFSQVCTHLDTLK",
      )
    self.assertEqual(
      hits[0].alignment.alignments[1],
      "MEEKLKKSKIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRAEVSSGSARGKMLSEIMEKGQLVPLETVLDMLRDAM"
      + "VAKVDTSKGFLIDGYPREVKQGEEFERKIGQPTLLLYVDAGPETMTKRLLKRGETSGRVDDNEETIKKRLETYYKATEPV"
      + "IAFYEKRGIVRKVNAEGSVDDVFSQVCTHLDTLK",
      )

    self.assertEqual( hits[1].identifier, "1Z83" )
    self.assertEqual( hits[1].chain, "A" )
    self.assertEqual( hits[1].annotation, "Second hit" )
    self.assertEqual( hits[1].alignment.multiplicity(), 2 )
    self.assertEqual( hits[1].alignment.names[0], "target" )
    self.assertEqual( hits[1].alignment.names[1], "1Z83_A" )
    self.assertEqual(
      hits[1].alignment.alignments[0],
      "MEEKLKKSKIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRAEVSSGSARGKMLSEIMEKGQLVPLETVLDMLRDAM"
      + "VAKVDTSKGFLIDGYPREVKQGEEFERKIGQPTLLLYVDAGPETMTKRLLKRGETSGRVDDNEETIKKRLETYYKATEPV"
      + "IAFYEKRGIVRKVNAEGSVDDVFSQVCTHLDTLK",
      )
    self.assertEqual(
      hits[1].alignment.alignments[1],
      "MEEKLKKTNIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRSEVSSGSARGKKLSEIMEKGQLVPLETVLDMLRDAM"
      + "VAKVNTSKGFLIDGYPREVQQGEEFERRIGQPTLLLYVDAGPETMTQRLLKRGETSGRVDDNEETIKKRLETYYKATEPV"
      + "IAFYEKRGIVRKVNAEGSVDSVFSQVCTHLDAL-"
      )


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


from __future__ import absolute_import, division, print_function

from iotbx.bioinformatics import ebi_wu_blast_xml

import unittest


class TestParse(unittest.TestCase):

  HSEARCH = None

  def setUp(self):

    if self.__class__.HSEARCH is None:
      import libtbx.load_env
      root = libtbx.env.under_dist( "iotbx", "bioinformatics" )

      import os.path
      import gzip
      self.__class__.HSEARCH = ebi_wu_blast_xml.parse(
        source = gzip.open( os.path.join( root, "test", "ebi_blast.xml.gz" ) )
        )


  def testAlignment(self):

    ali = self.HSEARCH.root.hits[6].alignments[0]

    self.assertEqual( ali.number, 1 )
    self.assertEqual( ali.score, 605 )
    self.assertAlmostEqual( ali.bits, 218.0, 1 )
    self.assertAlmostEqual( ali.expectation, 2.0e-59, 10 )
    self.assertAlmostEqual( ali.probability, 2.0e-59, 10 )
    self.assertEqual( ali.identity, 58 )
    self.assertEqual( ali.positives, 79 )
    self.assertEqual( ali.query.start, 4 )
    self.assertEqual( ali.query.end, 194 )
    self.assertEqual(
      ali.query.seq,
      "EKLKKSKIIFVVGGPGSGKGTQCEKIVQKYGYTHLSTGDLLRAEVSSGSARGKMLSEIMEKGQLVP"
      + "LETVLDMLRDAMVAKVDTSKGFLIDGYPREVKQGEEFERKIGQPTLLLYVDAGPETMTKRLLKR"
      + "GETSGRVDDNEETIKKRLETYYKATEPVIAFYEKRGIVRKVNAEGSVDDVFSQVCTHLDTL"
      )
    self.assertEqual( ali.match.start, 7 )
    self.assertEqual( ali.match.end, 197 )
    self.assertEqual(
      ali.match.seq,
      "EDLRKCKIIFIIGGPGSGKGTQCEKLVEKYGFTHLSTGELLREELASESERSKLIRDIMERGDLVP"
      + "SGIVLELLKEAMVASLGDTRGFLIDGYPREVKQGEEFGRRIGDPQLVICMDCSADTMTNRLLQM"
      + "SRSSLPVDDTTKTIAKRLEAYYRASIPVIAYYETKTQLHKINAEGTPEDVFLQLCTAIDSI"
      )
    self.assertEqual(
      ali.pattern,
      "E L+K KIIF++GGPGSGKGTQCEK+V+KYG+THLSTG+LLR E++S S R K++ +IME+G LVP"
      + "   VL++L++AMVA +  ++GFLIDGYPREVKQGEEF R+IG P L++ +D   +TMT RLL+ "
      + "  +S  VDD  +TI KRLE YY+A+ PVIA+YE +  + K+NAEG+ +DVF Q+CT +D++"
      )


  def testHit(self):

    hit = self.HSEARCH.root.hits[0]

    self.assertEqual( hit.number, 1 )
    self.assertEqual( hit.database, "pdb" )
    self.assertEqual( hit.identifier, "3ADK_A" )
    self.assertEqual( hit.length, 195 )
    self.assertEqual( hit.description, "mol:protein length:195  ADENYLATE KINASE" )
    self.assertEqual( len( hit.alignments ), 1 )


  def testRoot(self):

    header = self.HSEARCH.root.header

    self.assertEqual( header.program.name, "WU-blastp" )
    self.assertEqual( header.program.version, "2.0MP-WashU [04-May-2006]" )
    self.assertEqual( header.program.citation, "PMID:12824421" )
    self.assertEqual(
      header.command,
      "/ebi/extserv/bin/wu-blast/blastp \"pdb\" /ebi/extserv/blast-work/interactive/blast-20090331-16013"
      + "83158.input E=10 B=50 V=100 -mformat=1 -matrix BLOSUM62 -sump  -filter seg -cpus 8 -sort_by_pvalue -putenv='"
      + "WUBLASTMAT=/ebi/extserv/bin/wu-blast/matrix' -putenv=\"WUBLASTDB=$IDATA_CURRENT/blastdb\" -putenv="
      + "'WUBLASTFILTER=/ebi/extserv/bin/wu-blast/filter' "
      )
    self.assertEqual( header.time.start, "2009-03-31T16:01:44+01:00" )
    self.assertEqual( header.time.end, "2009-03-31T16:01:45+01:00" )
    self.assertEqual( header.time.search, "PT01S" )
    self.assertEqual( len( header.parameters.sequences ), 1 )

    seq = header.parameters.sequences[0]
    self.assertEqual( seq.number, 1 )
    self.assertEqual( seq.name, "Sequence" )
    self.assertEqual( seq.type, "p" )
    self.assertEqual( seq.length, 195 )
    self.assertEqual( len( header.parameters.databases ), 1 )

    db = header.parameters.databases[0]
    self.assertEqual( db.number, 1 )
    self.assertEqual( db.name, "pdb" )
    self.assertEqual( db.type, "p" )
    self.assertEqual( db.created, "2009-03-27T00:00:07+01:00" )

    self.assertEqual( header.parameters.scores, 100.0 )
    self.assertEqual( header.parameters.alignments, 50 )
    self.assertEqual( header.parameters.matrix, "BLOSUM62" )
    self.assertEqual( header.parameters.expectation_upper, 10.0 )
    self.assertEqual( header.parameters.statistics, "sump" )
    self.assertEqual( header.parameters.filter, "seg" )

    hits = self.HSEARCH.root.hits
    self.assertEqual( len( hits ), 50 )


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


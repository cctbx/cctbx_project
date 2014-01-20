from __future__ import division

from mmtbx.geometry import topology

import unittest


class TestSidechainMatch(unittest.TestCase):

  def test_asn_leu(self):

    l_ca = topology.Atom( label = "CA", coordinates = ( -1.0085, -0.590773,  0.814318 ) )
    l_cb = topology.Atom( label = "C",  coordinates = (  0.0275, -0.557773, -0.314682 ) )
    l_cg = topology.Atom( label = "C",  coordinates = (  1.2335,  0.374227, -0.138682 ) )
    l_cd1 = topology.Atom( label = "C", coordinates = (  2.3065,  0.046227, -1.16768  ) )
    l_cd2 = topology.Atom( label = "C", coordinates = (  0.8395,  1.84323,  -0.230682 ) )
    leu = topology.Molecule()
    leu.add( atom = l_ca )
    leu.add( atom = l_cb )
    leu.add( atom = l_cg )
    leu.add( atom = l_cd1 )
    leu.add( atom = l_cd2 )

    a_ca = topology.Atom( label = "CA", coordinates = ( -1.03327, -0.544348,  0.860946 ) )
    a_cb = topology.Atom( label = "C",  coordinates = (  0.10486, -0.548357, -0.164901 ) )
    a_cg = topology.Atom( label = "C",  coordinates = (  0.990984, 0.682823, -0.070521 ) )
    a_od1 = topology.Atom( label = "C", coordinates = (  1.39496,  1.24684,  -1.08724  ) )
    a_nd2 = topology.Atom( label = "C", coordinates = (  1.29745,  1.10599,   1.15228  ) )
    asn = topology.Molecule()
    asn.add( atom = a_ca )
    asn.add( atom = a_cb )
    asn.add( atom = a_cg )
    asn.add( atom = a_od1 )
    asn.add( atom = a_nd2 )

    res = topology.sidechain_match( molecule1 = leu, molecule2 = asn, tolerance = 0.1 )
    self.assertEqual( len( res ), 3 )
    self.assertTrue( ( l_ca, a_ca ) in res )
    self.assertTrue( ( l_cb, a_cb ) in res )
    self.assertTrue( ( l_cg, a_cg ) in res )
    self.assertTrue( ( l_cd1, a_od1 ) not in res )


suite_sidechain_match= unittest.TestLoader().loadTestsFromTestCase(
  TestSidechainMatch
  )


alltests = unittest.TestSuite(
  [
    suite_sidechain_match,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )


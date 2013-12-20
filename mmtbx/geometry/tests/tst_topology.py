from __future__ import division

from mmtbx.geometry import topology

import unittest

class collector(object):

  def __init__(self):

    self.collected = []


  def __call__(self, data):

    self.collected.append( data )
    return True


class TestGraph(unittest.TestCase):

  def setUp(self):

    self.graph = topology.graph()


  def test_manipulation(self):

    self.assertEqual( self.graph.add_vertex( label = "CA" ), 0 )
    self.assertEqual( self.graph.add_vertex( label = "C" ), 1 )
    self.assertEqual( self.graph.add_vertex( label = "C" ), 2 )

    res = self.graph.add_edge( vertex1 = 0, vertex2 = 1, weight = 1.5 )
    self.assertEqual( len( res ), 2 )
    self.assertTrue( res[1] )

    res = self.graph.add_edge( vertex1 = 0, vertex2 = 1, weight = 1.5 )
    self.assertEqual( len( res ), 2 )
    self.assertFalse( res[1] )

    res = self.graph.add_edge( vertex1 = 0, vertex2 = 2, weight = 1.5 )
    self.assertEqual( len( res ), 2 )
    self.assertTrue( res[1] )


  def test_matching(self):

    leu = topology.graph()
    v0 = leu.add_vertex( "CA" )
    v1 = leu.add_vertex( "C" )  # CB
    v2 = leu.add_vertex( "C" )  # CG
    v3 = leu.add_vertex( "C" )  # CD2
    v4 = leu.add_vertex( "C" )  # CD1

    leu.add_edge(v0, v1, 1.53 )
    leu.add_edge(v0, v2, 2.62 )
    leu.add_edge(v0, v3, 3.23 )
    leu.add_edge(v0, v4, 3.91 )

    leu.add_edge(v1, v2, 1.53 )
    leu.add_edge(v1, v3, 2.54 )
    leu.add_edge(v1, v4, 2.51 )

    leu.add_edge(v2, v3, 1.52 )
    leu.add_edge(v2, v4, 1.52 )

    leu.add_edge(v3, v4, 2.50 )

    asn = topology.graph()
    w0 = asn.add_vertex( "CA" )
    w1 = asn.add_vertex( "C" )  # CB
    w2 = asn.add_vertex( "C" )  # CG
    w3 = asn.add_vertex( "C" )  # ND2
    w4 = asn.add_vertex( "C" )  # OD1

    asn.add_edge(w0, w1, 1.53 )
    asn.add_edge(w0, w2, 2.54 )
    asn.add_edge(w0, w3, 2.87 )
    asn.add_edge(w0, w4, 3.59 )

    asn.add_edge(w1, w2, 1.52 )
    asn.add_edge(w1, w3, 2.43 )
    asn.add_edge(w1, w4, 2.40 )

    asn.add_edge(w2, w3, 1.33 )
    asn.add_edge(w2, w4, 1.23 )

    asn.add_edge(w3, w4, 2.25 )

    callback = collector()

    import operator

    topology.mcgregor_common_subgraphs_unique(
      graph1 = leu,
      graph2 = asn,
      vertex_equality = operator.eq,
      edge_equality = lambda l, r: abs( l - r ) <= 0.1,
      callback = callback,
      )
    self.assertEqual( len( callback.collected ), 13 )
    self.assertEqual( max( len( m ) for m in callback.collected ), 3 )


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


suite_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestGraph
  )
suite_sidechain_match= unittest.TestLoader().loadTestsFromTestCase(
  TestSidechainMatch
  )


alltests = unittest.TestSuite(
  [
    suite_graph,
    suite_sidechain_match,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )


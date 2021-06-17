from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import graph_structure_comparison

import unittest

class collector(object):

  def __init__(self):

    self.collected = []


  def __call__(self, data):

    self.collected.append( data )
    return True


class TestMcGregorCommonSubgraphsUnique(unittest.TestCase):

  def build_leu(self, leu):

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

    return ( v0, v1, v2 )


  def build_asn(self, asn):

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

    return ( w0, w1, w2 )


  def manipulate(self, leu, asn, matchings):

    ( leu_ca, leu_cb, leu_cg ) = self.build_leu( leu = leu )
    ( asn_ca, asn_cb, asn_cg ) = self.build_asn( asn = asn )

    callback = collector()

    import operator

    graph_structure_comparison.mcgregor_common_subgraphs_unique(
      graph1 = leu,
      graph2 = asn,
      vertex_equality = operator.eq,
      edge_equality = lambda l, r: abs( l - r ) <= 0.1,
      callback = callback,
      )
    self.assertEqual( len( callback.collected ), matchings )

    largest = max( callback.collected, key = len )
    self.assertEqual( len( largest ), 3 )
    self.assertEqual(
      sorted( [ ( leu.vertex_label( p[0] ), asn.vertex_label( p[1] ) ) for p in largest ] ),
      [ ( "C", "C" ), ( "C", "C" ), ( "CA", "CA" ) ],
      )
    self.assertEqual(
      set( largest ),
      set( [ ( leu_ca, asn_ca ), ( leu_cb, asn_cb ), ( leu_cg, asn_cg ) ] ),
      )



  def test_adjacency_list_undirected_vector_set(self):

    try:
      leu = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "set",
        )
      asn = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulate( leu = leu, asn = asn, matchings = 13 )


  def test_adjacency_list_undirected_list_set(self):

    try:
      leu = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "set",
        )
      asn = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulate( leu = leu, asn = asn, matchings = 13 )


  def test_adjacency_list_undirected_vector_vector(self):

    try:
      leu = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "vector",
        )
      asn = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulate( leu = leu, asn = asn, matchings = 13 )


  def test_adjacency_list_undirected_list_vector(self):

    try:
      leu = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "vector",
        )
      asn = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulate( leu = leu, asn = asn, matchings = 13 )


  def test_adjacency_list_directed_vector_vector(self):

    try:
      leu = graph.adjacency_list(
        graph_type = "directed",
        vertex_type = "vector",
        edge_type = "vector",
        )
      asn = graph.adjacency_list(
        graph_type = "directed",
        vertex_type = "vector",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulate( leu = leu, asn = asn, matchings = 8 )


suite_mcgregor_common_subgraphs_unique = unittest.TestLoader().loadTestsFromTestCase(
  TestMcGregorCommonSubgraphsUnique
  )


alltests = unittest.TestSuite(
  [
    suite_mcgregor_common_subgraphs_unique,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

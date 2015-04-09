from __future__ import division

from boost_adaptbx import graph
from graph import minimum_cut_enumerate

import unittest

class TestStirlingAdaptiveEnumeration(unittest.TestCase):

  def manipulation(self, g):

    vd0 = g.add_vertex( label = 0 )
    vd1 = g.add_vertex( label = 1 )
    vd2 = g.add_vertex( label = 2 )
    vd3 = g.add_vertex( label = 3 )
    vd4 = g.add_vertex( label = 4 )
    vd5 = g.add_vertex( label = 5 )
    vd6 = g.add_vertex( label = 6 )
    vd7 = g.add_vertex( label = 7 )

    g.add_edge( vertex1 = vd3, vertex2 = vd4, weight = 4 )
    g.add_edge( vertex1 = vd3, vertex2 = vd6, weight = 3 )
    g.add_edge( vertex1 = vd3, vertex2 = vd5, weight = 1 )
    g.add_edge( vertex1 = vd0, vertex2 = vd4, weight = 3 )
    g.add_edge( vertex1 = vd0, vertex2 = vd1, weight = 1 )
    g.add_edge( vertex1 = vd0, vertex2 = vd6, weight = 2 )
    g.add_edge( vertex1 = vd0, vertex2 = vd7, weight = 6 )
    g.add_edge( vertex1 = vd0, vertex2 = vd5, weight = 1 )
    g.add_edge( vertex1 = vd0, vertex2 = vd2, weight = 8 )
    g.add_edge( vertex1 = vd4, vertex2 = vd1, weight = 1 )
    g.add_edge( vertex1 = vd1, vertex2 = vd6, weight = 1 )
    g.add_edge( vertex1 = vd1, vertex2 = vd5, weight = 80 )
    g.add_edge( vertex1 = vd6, vertex2 = vd7, weight = 2 )
    g.add_edge( vertex1 = vd7, vertex2 = vd5, weight = 1 )
    g.add_edge( vertex1 = vd5, vertex2 = vd2, weight = 1 )

    stiter = minimum_cut_enumerate.stirling_adaptive_tree_enumeration(
      graph = g,
      reference = vd0,
      sw_path_vertex_selector = minimum_cut_enumerate.arbitrary_selection_from_set,
      bk_path_vertex_selector = minimum_cut_enumerate.BKEqualSizeSelector(),
      )
    cuts_with_weight = {}

    for ( weight, source, sink ) in stiter:
      solu = (
        frozenset( g.vertex_label( vertex = v ) for v in source ),
        frozenset( g.vertex_label( vertex = v ) for v in sink ),
        )
      cuts_with_weight.setdefault( int( weight ), set() ).add( solu )

    count_with_weight = {
      7: 1,
      8: 4,
      9: 2,
      10: 2,
      11: 2,
      13: 5,
      14: 3,
      15: 2,
      16: 3,
      17: 7,
      18: 4,
      19: 3,
      20: 6,
      21: 3,
      22: 2,
      23: 2,
      24: 2,
      25: 2,
      26: 4,
      27: 2,
      29: 1,
      30: 1,
      83: 1,
      84: 1,
      89: 4,
      90: 2,
      91: 4,
      92: 6,
      94: 2,
      95: 4,
      96: 1,
      97: 4,
      98: 7,
      99: 5,
      100: 5,
      101: 1,
      102: 3,
      103: 3,
      104: 3,
      105: 1,
      106: 1,
      107: 3,
      109: 2,
      110: 1,
      }

    self.assertEqual( len( cuts_with_weight ), len( count_with_weight ) )

    for ( weight, cuts ) in cuts_with_weight.iteritems():
      self.assertEqual( len( cuts ), count_with_weight[ weight ] )


  def test_adjacency_list_undirected_vector_set(self):

    try:
      g = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( g )


  def test_adjacency_list_undirected_list_set(self):

    try:
      g = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( g )


  def test_adjacency_list_undirected_vector_vector(self):

    try:
      g = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( g )


  def test_adjacency_list_undirected_list_vector(self):

    try:
      g = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( g )


suite_stirling_adaptive_enumeration = unittest.TestLoader().loadTestsFromTestCase(
  TestStirlingAdaptiveEnumeration
  )

alltests = unittest.TestSuite(
  [
    suite_stirling_adaptive_enumeration,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca

import unittest


class TestConnectedComponents(unittest.TestCase):

  def manipulation(self, g):

    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()

    g.add_edge( vertex1 = vd1, vertex2 = vd2 )
    components = cca.connected_components( graph = g )

    self.assertEqual( len( components ), 2 )
    self.assertEqual(
      set( [ frozenset( c ) for c in components ] ),
      set( [ frozenset( [ vd1, vd2 ] ), frozenset( [ vd3 ] ) ] ),
      )


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

suite_connected_components = unittest.TestLoader().loadTestsFromTestCase(
  TestConnectedComponents
  )


alltests = unittest.TestSuite(
  [
    suite_connected_components,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

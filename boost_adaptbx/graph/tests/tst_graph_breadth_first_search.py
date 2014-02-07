from __future__ import division

from boost_adaptbx import graph
from graph import breadth_first_search as bfs

import unittest


class TestBreadthFirstSearch(unittest.TestCase):

  def manipulation(self, g):

    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()
    vd4 = g.add_vertex()

    g.add_edge( vertex1 = vd1, vertex2 = vd2 )
    g.add_edge( vertex1 = vd1, vertex2 = vd3 )
    g.add_edge( vertex1 = vd2, vertex2 = vd4 )
    visitor = bfs.distance_recording_visitor( start_vertex = vd1 )
    bfs.breadth_first_search( graph = g, vertex = vd1, visitor = visitor )

    self.assertEqual(
      visitor.distance_for,
      { vd1: 0, vd2: 1, vd3: 1, vd4: 2 }
      )
    visitor2 = bfs.distance_recording_visitor( start_vertex = vd2 )
    bfs.breadth_first_search( graph = g, vertex = vd2, visitor = visitor2 )

    self.assertEqual(
      visitor2.distance_for,
      { vd1: 1, vd2: 0, vd3: 2, vd4: 1 }
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

suite_breadth_first_search = unittest.TestLoader().loadTestsFromTestCase(
  TestBreadthFirstSearch
  )


alltests = unittest.TestSuite(
  [
    suite_breadth_first_search,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

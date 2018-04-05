from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import breadth_first_search as bfs

import unittest


class TestVertexRecordingVisitor(unittest.TestCase):

  def manipulation(self, g):

    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()
    vd4 = g.add_vertex()

    g.add_edge( vertex1 = vd1, vertex2 = vd2 )
    visitor = bfs.vertex_recording_visitor( start_vertex = vd1 )
    bfs.breadth_first_search( graph = g, vertex = vd1, visitor = visitor )
    self.assertEqual( visitor.visited_vertices, [ vd1, vd2 ] )

    visitor = bfs.vertex_recording_visitor( start_vertex = vd2 )
    bfs.breadth_first_search( graph = g, vertex = vd2, visitor = visitor )
    self.assertEqual( visitor.visited_vertices, [ vd2, vd1 ] )

    visitor = bfs.vertex_recording_visitor( start_vertex = vd3 )
    bfs.breadth_first_search( graph = g, vertex = vd3, visitor = visitor )
    self.assertEqual( visitor.visited_vertices, [ vd3 ] )

    g.add_edge( vertex1 = vd2, vertex2 = vd3 )
    g.add_edge( vertex1 = vd2, vertex2 = vd4 )
    visitor = bfs.vertex_recording_visitor( start_vertex = vd1 )
    bfs.breadth_first_search( graph = g, vertex = vd1, visitor = visitor )
    self.assertEqual( visitor.visited_vertices[:2], [ vd1, vd2 ] )
    self.assertEqual(
      set( visitor.visited_vertices ),
      set( [ vd1, vd2, vd3, vd4 ] ),
      )

    visitor = bfs.vertex_recording_visitor( start_vertex = vd2 )
    bfs.breadth_first_search( graph = g, vertex = vd2, visitor = visitor )
    self.assertEqual(
      set( visitor.visited_vertices ),
      set( [ vd1, vd2, vd3, vd4 ] ),
      )
    self.assertEqual( visitor.visited_vertices[0], vd2 )

    g.add_edge( vertex1 = vd1, vertex2 = vd3 )
    visitor = bfs.vertex_recording_visitor( start_vertex = vd1 )
    bfs.breadth_first_search( graph = g, vertex = vd1, visitor = visitor )
    self.assertEqual(
      set( visitor.visited_vertices ),
      set( [ vd1, vd2, vd3, vd4 ] ),
      )
    self.assertEqual( len( visitor.visited_vertices ), 4 )


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


class TestDistanceRecordingVisitor(unittest.TestCase):

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


suite_vertex_recording_visitor = unittest.TestLoader().loadTestsFromTestCase(
  TestVertexRecordingVisitor
  )
suite_distance_recording_visitor = unittest.TestLoader().loadTestsFromTestCase(
  TestDistanceRecordingVisitor
  )


alltests = unittest.TestSuite(
  [
    suite_vertex_recording_visitor,
    suite_distance_recording_visitor,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

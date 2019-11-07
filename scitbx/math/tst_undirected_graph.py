from __future__ import absolute_import, division, print_function

from scitbx.math import undirected_graph

import unittest

class TestGraph(unittest.TestCase):

  def setUp(self):

    self.edges = [
      undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ),
      undirected_graph.Edge( vertex1 = 2, vertex2 = 3 ),
      undirected_graph.Edge( vertex1 = 2, vertex2 = 4 ),
      undirected_graph.Edge( vertex1 = 3, vertex2 = 4 ),
      ]
    self.add_vertices = [ 5, 1, 6 ]

    self.vertices = set( [ 1, 2, 3, 4, 5, 6 ] )

    self.edges_from = {
      1: set( [ self.edges[0] ] ),
      2: set( [ self.edges[0], self.edges[1], self.edges[2] ] ),
      3: set( [ self.edges[1], self.edges[3] ] ),
      4: set( [ self.edges[2], self.edges[3] ] ),
      5: set(),
      6: set(),
      }


  def build_and_test_graph(self, graph):

    self.assertEqual( set( graph.vertices ), set() )
    self.assertEqual( set( graph.edges ), set() )

    for e in self.edges:
      graph.add_edge( edge = e )

    for v in self.vertices:
      graph.add_vertex( vertex = v )

    self.assertEqual( set( graph.edges ), set( self.edges ) )
    self.assertEqual( set( graph.vertices ), self.vertices )

    for ( vertex, edges ) in self.edges_from.items():
      self.assertEqual( set( graph.edges_from( vertex = vertex ) ), edges )


  def test_graph(self):

    self.build_and_test_graph( graph = undirected_graph.Graph() )


  def test_vertex_indexed_graph(self):

    self.build_and_test_graph( graph = undirected_graph.VertexIndexedGraph() )


class TestPreorderDepthFirstIteration(unittest.TestCase):

  def create_graph_from_edges(self, edges):

    graph = undirected_graph.Graph()

    for e in edges:
      graph.add_edge( edge = e )

    return graph


  def test_1(self):

    graph = self.create_graph_from_edges(
      edges = [
        undirected_graph.Edge( vertex1 = 1, vertex2 = 1 ),
        ],
      )
    self.assertEqual(
      list( undirected_graph.preorder_depth_first_iteration(vertex = 1, graph = graph ) ),
      [ 1 ],
      )


  def test_2(self):

    graph = self.create_graph_from_edges(
      edges = [
        undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ),
        ],
        )
    self.assertEqual(
      list( undirected_graph.preorder_depth_first_iteration(vertex = 1, graph = graph ) ),
      [ 1, 2 ],
      )
    self.assertEqual(
      list( undirected_graph.preorder_depth_first_iteration(vertex = 2, graph = graph ) ),
      [ 2, 1 ],
      )


  def test_3(self):

    graph = self.create_graph_from_edges(
      edges = [
        undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ),
        undirected_graph.Edge( vertex1 = 2, vertex2 = 3 ),
        ],
      )
    self.assertEqual(
      list( undirected_graph.preorder_depth_first_iteration(vertex = 1, graph = graph ) ),
      [ 1, 2, 3 ],
      )
    self.assertEqual(
      list( undirected_graph.preorder_depth_first_iteration(vertex = 3, graph = graph ) ),
      [ 3, 2, 1 ],
      )
    middle = list(
      undirected_graph.preorder_depth_first_iteration(vertex = 2, graph = graph )
      )
    self.assertTrue( ( middle == [ 2, 1, 3 ] ) or ( middle == [ 2, 3, 1 ] ) )


  def test_4(self):

    graph = self.create_graph_from_edges(
      edges = [
        undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ),
        undirected_graph.Edge( vertex1 = 2, vertex2 = 3 ),
        undirected_graph.Edge( vertex1 = 2, vertex2 = 4 ),
        undirected_graph.Edge( vertex1 = 3, vertex2 = 4 ),
        ],
      )

    res = list(
      undirected_graph.preorder_depth_first_iteration(vertex = 1, graph = graph )
      )
    self.assertTrue( ( res == [ 1, 2, 3, 4 ] ) or ( res == [ 1, 2, 4, 3 ] ) )

    res = list(
      undirected_graph.preorder_depth_first_iteration(vertex = 2, graph = graph )
      )
    self.assertTrue(
      ( res == [ 2, 1, 3, 4 ] )
      or ( res == [ 2, 1, 4, 3 ] )
      or ( res == [ 2, 3, 4, 1 ] )
      or ( res == [ 2, 4, 3, 1 ] )
      )

    res = list(
      undirected_graph.preorder_depth_first_iteration(vertex = 3, graph = graph )
      )
    self.assertTrue(
      ( res == [ 3, 2, 1, 4 ] )
      or ( res == [ 3, 2, 4, 1 ] )
      or ( res == [ 3, 4, 2, 1 ] )
      )
    res = list(
      undirected_graph.preorder_depth_first_iteration(vertex = 4, graph = graph )
      )
    self.assertTrue(
      ( res == [ 4, 3, 2, 1 ] )
      or ( res == [ 4, 2, 1, 3 ] )
      or ( res == [ 4, 2, 3, 1 ] )
      )


class TestConnectedComponents(unittest.TestCase):

  def test_1(self):

    graph = undirected_graph.Graph()
    graph.add_edge( edge = undirected_graph.Edge( vertex1 = 1, vertex2 = 1 ) )

    self.assertEqual(
      list( undirected_graph.connected_components( graph = graph ) ),
      [ [ 1 ] ],
      )

  def test_2(self):

    graph = undirected_graph.Graph()
    graph.add_edge( edge = undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ) )

    self.assertEqual(
      list(
        sorted( c ) for c in undirected_graph.connected_components( graph = graph )
        ),
      [ [ 1, 2 ] ],
      )

    graph = undirected_graph.Graph()
    graph.add_vertex( 1 )
    graph.add_vertex( 2 )

    self.assertEqual(
      sorted( undirected_graph.connected_components( graph = graph ) ),
      [ [ 1 ], [ 2 ] ],
      )

  def test_3(self):

    graph = undirected_graph.Graph()
    graph.add_edge( undirected_graph.Edge( vertex1 = 1, vertex2 = 2 ) )
    graph.add_edge( undirected_graph.Edge( vertex1 = 2, vertex2 = 3 ) )
    graph.add_edge( undirected_graph.Edge( vertex1 = 4, vertex2 = 5 ) )
    graph.add_edge( undirected_graph.Edge( vertex1 = 5, vertex2 = 6 ) )
    graph.add_edge( undirected_graph.Edge( vertex1 = 5, vertex2 = 7 ) )
    graph.add_edge( undirected_graph.Edge( vertex1 = 6, vertex2 = 7 ) )

    self.assertEqual(
      sorted( sorted( c ) for c in undirected_graph.connected_components( graph = graph ) ),
      [ [ 1, 2, 3 ], [ 4, 5, 6, 7 ] ],
      )


suite_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestGraph
  )
suite_preorder_depth_first = unittest.TestLoader().loadTestsFromTestCase(
  TestPreorderDepthFirstIteration
  )
suite_connected_components = unittest.TestLoader().loadTestsFromTestCase(
  TestConnectedComponents
  )

alltests = unittest.TestSuite(
  [
    suite_graph,
    suite_preorder_depth_first,
    suite_connected_components,
    ]
  )

if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )


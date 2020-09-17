from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph

import unittest


class TestGraph(unittest.TestCase):

  def manipulation(self, g, unique_edge, directed):

    vd1 = g.add_vertex( label = "CA" )
    vd2 = g.add_vertex( label = "C" )
    vd3 = g.add_vertex( label = "C" )

    res = g.add_edge( vertex1 = vd1, vertex2 = vd2, weight = 1 )
    self.assertEqual( len( res ), 2 )
    self.assertTrue( res[1] )
    ed1 = res[0]

    res = g.add_edge( vertex1 = vd2, vertex2 = vd3, weight = 2 )
    self.assertEqual( len( res ), 2 )
    self.assertTrue( res[1] )
    ed2 = res[0]

    if unique_edge:
      res = g.add_edge( vertex1 = vd1, vertex2 = vd2, weight = 1 )
      self.assertEqual( len( res ), 2 )
      self.assertFalse( res[1] )

    self.assertEqual( set( g.vertices() ), set( [ vd1, vd2, vd3 ] ) )
    self.assertEqual( set( g.edges() ), set( [ ed1, ed2 ] ) )

    self.assertEqual( set( g.adjacent_vertices( vertex = vd1 ) ), set( [ vd2 ] ) )
    self.assertEqual(
      set( g.adjacent_vertices( vertex = vd2 ) ),
      set( [ vd3 ] ) if directed else set( [ vd1, vd3 ] ),
      )
    self.assertEqual(
      set( g.adjacent_vertices( vertex = vd3 ) ),
      set() if directed else set( [ vd2 ] ),
      )

    self.assertEqual( set( g.out_edges( vertex = vd1 ) ), set( [ ed1 ] ) )
    self.assertEqual(
      set( g.out_edges( vertex = vd2 ) ),
      set( [ ed2 ] ) if directed else set( [ ed1, ed2 ] ),
      )
    self.assertEqual(
      set( g.out_edges( vertex = vd3 ) ),
      set() if directed else set( [ ed2 ] ),
      )

    self.assertEqual( g.vertex_label( vertex = vd1 ), "CA" )
    self.assertEqual( g.vertex_label( vertex = vd2 ), "C" )
    self.assertEqual( g.vertex_label( vertex = vd3 ), "C" )
    self.assertEqual( g.edge_weight( edge = ed1 ), 1 )
    self.assertEqual( g.edge_weight( edge = ed2 ), 2 )

    self.assertEqual( g.source( edge = ed1 ), vd1 )
    self.assertEqual( g.source( edge = ed2 ), vd2 )
    self.assertEqual( g.target( edge = ed1 ), vd2 )
    self.assertEqual( g.target( edge = ed2 ), vd3 )

    g.set_vertex_label( vertex = vd1, label = "CD" )
    self.assertEqual( g.vertex_label( vertex = vd1 ), "CD" )
    g.set_vertex_label( vertex = vd1, label = None )
    self.assertEqual( g.vertex_label( vertex = vd1 ), None )

    g.set_edge_weight( edge = ed1, weight = "FOO" )
    self.assertEqual( g.edge_weight( edge = ed1 ), "FOO" )

    g.remove_edge( edge = ed2 )
    self.assertEqual( len( list( g.edges() ) ), 1 )
    self.assertEqual( g.edge_weight( edge = ed1 ), "FOO" )

    g.remove_vertex( vertex = vd3 )
    self.assertEqual( set( g.vertices() ), set( [ vd1, vd2 ] ) )


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
      self.manipulation( g, unique_edge = True, directed = False )


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
      self.manipulation( g, unique_edge = True, directed = False )


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
      self.manipulation( g, unique_edge = False, directed = False )


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
      self.manipulation( g, unique_edge = False, directed = False )


  def test_adjacency_list_directed_vector_vector(self):

    try:
      g = graph.adjacency_list(
        graph_type = "directed",
        vertex_type = "vector",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( g, unique_edge = False, directed = True )


suite_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestGraph
  )


alltests = unittest.TestSuite(
  [
    suite_graph,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

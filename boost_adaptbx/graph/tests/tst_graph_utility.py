from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import utility

import unittest

class TestCopyGraph(unittest.TestCase):

  def manipulation(self, g):

    vd0 = g.add_vertex( label = object() )
    vd1 = g.add_vertex( label = object() )
    vd2 = g.add_vertex( label = object() )

    g.add_edge( vertex1 = vd0, vertex2 = vd1, weight = object() )
    g.add_edge( vertex1 = vd0, vertex2 = vd2, weight = object() )
    g.add_edge( vertex1 = vd1, vertex2 = vd2, weight = object() )

    copy = utility.copy_graph( graph = g )
    self.assertEqual( copy.num_vertices(), 3 )
    self.assertEqual( copy.num_edges(), 3 )
    self.assertEqual(
      set( g.vertex_label( vertex = v ) for v in g.vertices() ),
      set( copy.vertex_label( vertex = v ) for v in copy.vertices() ),
      )
    self.assertEqual(
      set( g.edge_weight( edge = e ) for e in g.edges() ),
      set( copy.edge_weight( edge = e ) for e in copy.edges() ),
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


class TestCopyGraphAndMapVertices(unittest.TestCase):

  def manipulation(self, g):

    vd0 = g.add_vertex( label = object() )
    vd1 = g.add_vertex( label = object() )
    vd2 = g.add_vertex( label = object() )

    g.add_edge( vertex1 = vd0, vertex2 = vd1, weight = object() )
    g.add_edge( vertex1 = vd0, vertex2 = vd2, weight = object() )
    g.add_edge( vertex1 = vd1, vertex2 = vd2, weight = object() )

    ( copy, mapping ) = utility.copy_graph_and_map_vertices( graph = g )
    self.assertEqual( copy.num_vertices(), 3 )
    self.assertEqual( copy.num_edges(), 3 )

    self.assertEqual( len( mapping ), 3 )

    for vertex in g.vertices():
      self.assertTrue( vertex in mapping )
      self.assertEqual(
        g.vertex_label( vertex = vertex ),
        copy.vertex_label( vertex = mapping[ vertex ] ),
        )

    weight_of_edge_between = {}

    for edge in copy.edges():
      sv = copy.source( edge = edge )
      dv = copy.target( edge = edge )
      weight_of_edge_between[ frozenset( ( sv, dv ) ) ] = copy.edge_weight( edge = edge )

    for edge in g.edges():
      sv = copy.source( edge = edge )
      dv = copy.target( edge = edge )
      key = frozenset( ( mapping[ sv ], mapping[ dv ] ) )
      self.assertTrue( key in weight_of_edge_between )
      self.assertEqual( weight_of_edge_between[ key ], g.edge_weight( edge = edge ) )


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


class TestContractUndirectedGraph(unittest.TestCase):

  def setUp(self):

    import operator
    self.contract_undirected_graph = utility.undirected_graph_vertex_contraction(
      joint_vertex_label = operator.or_,
      joint_edge_weight = operator.add
      )


  def manipulation(self, g):

    vd0 = g.add_vertex( label = frozenset( [ 0 ] ) )
    vd1 = g.add_vertex( label = frozenset( [ 1 ] ) )
    vd2 = g.add_vertex( label = frozenset( [ 2 ] ) )
    vd3 = g.add_vertex( label = frozenset( [ 3 ] ) )
    vd4 = g.add_vertex( label = frozenset( [ 4 ] ) )

    g.add_edge( vertex1 = vd0, vertex2 = vd1, weight = 1 )
    g.add_edge( vertex1 = vd0, vertex2 = vd2, weight = 3 )
    g.add_edge( vertex1 = vd1, vertex2 = vd2, weight = 4 )
    g.add_edge( vertex1 = vd0, vertex2 = vd3, weight = 9 )
    g.add_edge( vertex1 = vd1, vertex2 = vd4, weight = 8 )
    g.add_edge( vertex1 = vd2, vertex2 = vd4, weight = 6 )

    self.assertEqual(
      self.contract_undirected_graph.joint_vertex_label(
        frozenset( [ 0 ] ),
        frozenset( [ 1 ] ),
        ),
      frozenset( [ 0, 1 ] ),
      )
    self.assertEqual( self.contract_undirected_graph.joint_edge_weight( 3, 4 ), 7 )

    self.contract_undirected_graph( graph = g, v1 = vd0, v2 = vd1 )
    vl01 = frozenset( [ 0, 1 ] )
    vl2 = frozenset( [ 2 ] )
    vl3 = frozenset( [ 3 ] )
    vl4 = frozenset( [ 4 ] )

    self.assertEqual( g.num_vertices(), 4 )
    vertices = list( g.vertices() )
    self.assertEqual( len( vertices ), 4 )
    vertex_labels = set( g.vertex_label( vertex = v ) for v in vertices )
    self.assertEqual( len( vertex_labels ), 4 )
    self.assertEqual( vertex_labels, set( [ vl01, vl2, vl3, vl4 ] ) )

    self.assertEqual( g.num_edges(), 4 )
    edges = list( g.edges() )
    self.assertEqual( len( edges ), 4 )

    weight_between = {
      frozenset( [ vl01, vl2 ] ): 7,
      frozenset( [ vl01, vl3 ] ): 9,
      frozenset( [ vl01, vl4 ] ): 8,
      frozenset( [ vl2, vl4 ] ): 6,
      }

    for ed in edges:
      key = frozenset(
        [
          g.vertex_label( vertex = g.source( edge = ed ) ),
          g.vertex_label( vertex = g.target( edge = ed ) ),
          ]
        )
      self.assertTrue( key in weight_between )
      self.assertEqual( g.edge_weight( edge = ed ), weight_between[ key ] )


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


suite_copy_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestCopyGraph
  )
suite_copy_graph_and_map_vertices = unittest.TestLoader().loadTestsFromTestCase(
  TestCopyGraphAndMapVertices
  )
suite_contract_undirected_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestContractUndirectedGraph
  )


alltests = unittest.TestSuite(
  [
    suite_copy_graph,
    suite_copy_graph_and_map_vertices,
    suite_contract_undirected_graph,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

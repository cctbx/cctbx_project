from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import maximum_clique

import unittest

class accumulator(object):

  def __init__(self):

    self.largest = 0
    self.sent = []

  def __call__(self, result):

    size = len( result )

    if self.largest < size:
      self.largest = size
      self.sent = [ result ]

    elif self.largest == size:
      self.sent.append( result )


class TestRascal(unittest.TestCase):

  def manipulation(self, g):

    n = {}
    n[ "1,1'" ] = g.add_vertex() # 0
    n[ "1,5'" ] = g.add_vertex() # 1
    n[ "2,2'" ] = g.add_vertex() # 2
    n[ "2,3'" ] = g.add_vertex() # 3
    n[ "3,2'" ] = g.add_vertex() # 4
    n[ "3,3'" ] = g.add_vertex() # 5
    n[ "4,4'" ] = g.add_vertex() # 6
    n[ "5,1'" ] = g.add_vertex() # 7
    n[ "5,5'" ] = g.add_vertex() # 8

    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,2'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,1'" ] )
    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,5'" ] )

    a = accumulator()
    maximum_clique.rascal( graph = g, callable = a )
    self.assertEqual(
      set(
        [
          frozenset( [ n[ "1,1'" ], n[ "3,3'" ] , n[ "4,4'" ], n[ "5,5'" ] ] ),
          frozenset( [ n[ "2,2'" ], n[ "3,3'" ] , n[ "4,4'" ], n[ "5,5'" ] ] ),
          ]
        ),
      set( frozenset( m ) for m in a.sent ),
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


class TestGreedy(unittest.TestCase):

  def manipulation(self, g):

    n = {}
    n[ "1,1'" ] = g.add_vertex() # 0
    n[ "1,5'" ] = g.add_vertex() # 1
    n[ "2,2'" ] = g.add_vertex() # 2
    n[ "2,3'" ] = g.add_vertex() # 3
    n[ "3,2'" ] = g.add_vertex() # 4
    n[ "3,3'" ] = g.add_vertex() # 5
    n[ "4,4'" ] = g.add_vertex() # 6
    n[ "5,1'" ] = g.add_vertex() # 7
    n[ "5,5'" ] = g.add_vertex() # 8

    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,2'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,1'" ] )
    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,5'" ] )

    result = maximum_clique.greedy( graph = g )
    self.assertEqual(
      set(
        [
          frozenset( [ n[ "1,1'" ], n[ "3,3'" ] , n[ "4,4'" ], n[ "5,5'" ] ] ),
          frozenset( [ n[ "2,2'" ], n[ "3,3'" ] , n[ "4,4'" ], n[ "5,5'" ] ] ),
          ]
        ),
      set( frozenset( m ) for m in result ),
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


class TestCompatibilityGraph(unittest.TestCase):

  def manipulation(self, left, right):

    u1 = left.add_vertex( label = "O" )
    u2 = left.add_vertex( label = "C" )
    u3 = left.add_vertex( label = "C" )
    u4 = left.add_vertex( label = "N" )
    u5 = left.add_vertex( label = "O" )

    left.add_edge( vertex1 = u1, vertex2 = u2, weight = 1 )
    left.add_edge( vertex1 = u2, vertex2 = u3, weight = 1 )
    left.add_edge( vertex1 = u3, vertex2 = u4, weight = 1 )
    left.add_edge( vertex1 = u3, vertex2 = u5, weight = 2 )

    v1 = right.add_vertex( label = "O" )
    v2 = right.add_vertex( label = "C" )
    v3 = right.add_vertex( label = "C" )
    v4 = right.add_vertex( label = "N" )
    v5 = right.add_vertex( label = "O" )

    right.add_edge( vertex1 = v1, vertex2 = v2, weight = 2 )
    right.add_edge( vertex1 = v2, vertex2 = v3, weight = 1 )
    right.add_edge( vertex1 = v3, vertex2 = v4, weight = 1 )
    right.add_edge( vertex1 = v3, vertex2 = v5, weight = 2 )

    result = maximum_clique.compatibility_graph(
      first = left,
      second = right,
      )

    vertices = list( result.vertices() )
    edges = list( result.edges() )

    self.assertEqual( len( vertices ), 9 )
    self.assertEqual( len( edges ), 16 )

    self.assertEqual(
      set( result.vertex_label( vertex = v ) for v in result.vertices() ),
      set(
        [
          ( u1, v1 ), ( u1, v5 ), ( u5, v1 ), ( u5, v5 ),
          ( u2, v2 ), ( u2, v3 ), ( u3, v2 ), ( u3, v3 ),
          ( u4, v4 ),
          ]
        )
      )
    edges_between = set(
      [
        frozenset(
          [
            result.vertex_label( vertex = result.source( edge = e ) ),
            result.vertex_label( vertex = result.target( edge = e ) ),
            ]
          )
        for e in result.edges()
        ]
      )
    self.assertEqual(
      edges_between,
      set(
        [
          frozenset( [ ( u1, v1 ), ( u3, v3 ) ] ),
          frozenset( [ ( u1, v1 ), ( u4, v4 ) ] ),
          frozenset( [ ( u1, v1 ), ( u5, v5 ) ] ),

          frozenset( [ ( u1, v5 ), ( u3, v2 ) ] ),
          frozenset( [ ( u1, v5 ), ( u4, v4 ) ] ),
          frozenset( [ ( u1, v5 ), ( u5, v1 ) ] ),

          frozenset( [ ( u2, v2 ), ( u3, v3 ) ] ),
          frozenset( [ ( u2, v2 ), ( u4, v4 ) ] ),
          frozenset( [ ( u2, v2 ), ( u5, v5 ) ] ),

          frozenset( [ ( u2, v3 ), ( u3, v2 ) ] ),
          frozenset( [ ( u2, v3 ), ( u5, v1 ) ] ),

          frozenset( [ ( u3, v2 ), ( u5, v1 ) ] ),

          frozenset( [ ( u3, v3 ), ( u4, v4 ) ] ),
          frozenset( [ ( u3, v3 ), ( u5, v5 ) ] ),

          frozenset( [ ( u4, v4 ), ( u5, v1 ) ] ),
          frozenset( [ ( u4, v4 ), ( u5, v5 ) ] ),
          ]
        )
      )


  def test_adjacency_list_undirected_vector_set(self):

    try:
      g1 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "set",
        )
      g2 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( left = g1, right = g2 )


  def test_adjacency_list_undirected_list_set(self):

    try:
      g1 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "set",
        )
      g2 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "set",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( left = g1, right = g2 )


  def test_adjacency_list_undirected_vector_vector(self):

    try:
      g1 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "vector",
        )
      g2 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "vector",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( left = g1, right = g2 )


  def test_adjacency_list_undirected_list_vector(self):

    try:
      g1 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "vector",
        )
      g2 = graph.adjacency_list(
        graph_type = "undirected",
        vertex_type = "list",
        edge_type = "vector",
        )

    except NotImplementedError:
      pass

    else:
      self.manipulation( left = g1, right = g2 )


class TestSelectedSubgraph(unittest.TestCase):

  def manipulation(self, g):

    n = {}
    n[ "1,1'" ] = g.add_vertex( label = "1,1'" ) # 0
    n[ "1,5'" ] = g.add_vertex( label = "1,5'" ) # 1
    n[ "2,2'" ] = g.add_vertex( label = "2,2'" ) # 2
    n[ "2,3'" ] = g.add_vertex( label = "2,3'" ) # 3
    n[ "3,2'" ] = g.add_vertex( label = "3,2'" ) # 4
    n[ "3,3'" ] = g.add_vertex( label = "3,3'" ) # 5
    n[ "4,4'" ] = g.add_vertex( label = "4,4'" ) # 6
    n[ "5,1'" ] = g.add_vertex( label = "5,1'" ) # 7
    n[ "5,5'" ] = g.add_vertex( label = "5,5'" ) # 8

    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,1'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "1,5'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "3,3'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "2,2'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "3,2'" ] )
    g.add_edge( vertex1 = n[ "2,3'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,2'" ], vertex2 = n[ "5,1'" ] )

    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "4,4'" ] )
    g.add_edge( vertex1 = n[ "3,3'" ], vertex2 = n[ "5,5'" ] )

    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,1'" ] )
    g.add_edge( vertex1 = n[ "4,4'" ], vertex2 = n[ "5,5'" ] )

    subgraph = maximum_clique.selected_subgraph(
      graph = g,
      vertices = g.adjacent_vertices( vertex = n[ "2,2'"] )
      )
    vertices = list( subgraph.vertices() )
    self.assertEqual( len( vertices ), 3 )
    self.assertEqual(
      set( subgraph.vertex_label( vertex = v ) for v in vertices ),
      set( [ "3,3'", "4,4'", "5,5'" ] )
      )

    edges = list( subgraph.edges() )
    self.assertEqual( len( edges ), 3 )
    self.assertEqual(
      set(
        frozenset(
          [
            subgraph.vertex_label( vertex = subgraph.source( edge = e ) ),
            subgraph.vertex_label( vertex = subgraph.target( edge = e ) ),
            ]
          )
        for e in edges
        ),
      set(
        [
          frozenset( [ "3,3'", "4,4'" ] ),
          frozenset( [ "3,3'", "5,5'" ] ),
          frozenset( [ "4,4'", "5,5'" ] ),
          ]
        )
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


suite_rascal = unittest.TestLoader().loadTestsFromTestCase(
  TestRascal
  )
suite_greedy = unittest.TestLoader().loadTestsFromTestCase(
  TestGreedy
  )
suite_compatibility_graph = unittest.TestLoader().loadTestsFromTestCase(
  TestCompatibilityGraph
  )
suite_selected_subgraph = unittest.TestLoader().loadTestsFromTestCase(
  TestSelectedSubgraph
  )


alltests = unittest.TestSuite(
  [
    suite_rascal,
    suite_greedy,
    suite_compatibility_graph,
    suite_selected_subgraph,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

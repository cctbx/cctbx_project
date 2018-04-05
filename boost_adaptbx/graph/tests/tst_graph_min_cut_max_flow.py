from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import min_cut_max_flow

import unittest

class TestStoerWagnerMinCut(unittest.TestCase):

  def manipulation(self, g):

    vd0 = g.add_vertex()
    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()
    vd4 = g.add_vertex()
    vd5 = g.add_vertex()
    vd6 = g.add_vertex()
    vd7 = g.add_vertex()

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

    result = min_cut_max_flow.stoer_wagner_min_cut( graph = g )
    self.assertAlmostEqual( result.weight, 7, 7 )
    self.assertEqual( g.num_vertices(), len( result.parities ) )
    ( source, sink ) = result.cutsets

    if len( sink ) < len( source ):
      ( source, sink ) = ( sink, source )

    self.assertEqual( len( source ), 2 )
    self.assertEqual( len( sink ), 6 )
    self.assertTrue( source, set( ( vd1, vd5 ) ) )
    self.assertTrue( sink, set( ( vd0, vd2, vd3, vd4, vd6, vd7 ) ) )


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


class TestBoykovKolmogorovMaxFlow(unittest.TestCase):

  def create_edge(self, graph, left, right, weight, reverse_edge_map):

    edf = graph.add_edge( vertex1 = left, vertex2 = right, weight = weight )[0]
    edb = graph.add_edge( vertex1 = right, vertex2 = left, weight = weight )[0]
    reverse_edge_map[ edf ] = edb
    reverse_edge_map[ edb ] = edf


  def manipulation(self, g):

    vd0 = g.add_vertex()
    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()
    vd4 = g.add_vertex()
    vd5 = g.add_vertex()
    vd6 = g.add_vertex()
    vd7 = g.add_vertex()

    revmap = {}

    self.create_edge( graph = g, left = vd3, right = vd4, weight = 4, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd3, right = vd6, weight = 3, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd3, right = vd5, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd4, weight = 3, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd1, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd6, weight = 2, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd7, weight = 6, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd5, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd0, right = vd2, weight = 8, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd4, right = vd1, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd1, right = vd6, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd1, right = vd5, weight = 80, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd6, right = vd7, weight = 2, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd7, right = vd5, weight = 1, reverse_edge_map = revmap )
    self.create_edge( graph = g, left = vd5, right = vd2, weight = 1, reverse_edge_map = revmap )

    result = min_cut_max_flow.boykov_kolmogorov_max_flow(
      graph = g,
      reverse_edge_map = revmap,
      source = vd0,
      sink = vd1,
      )
    self.assertAlmostEqual( result.weight, 7, 7 )

    ( source, sink ) = result.cutsets

    self.assertTrue( vd0 in source )
    self.assertTrue( vd1 in sink )
    self.assertEqual( len( source ), 6 )
    self.assertEqual( source, set( ( vd0, vd2, vd3, vd4, vd6, vd7 ) ) )
    self.assertEqual( len( sink ), 2 )
    self.assertEqual( sink, set( ( vd1, vd5 ) ) )


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
      self.manipulation( g )


class TestBoykovKolmogorovMinSTCut(unittest.TestCase):

  def manipulation(self, g):

    vd0 = g.add_vertex()
    vd1 = g.add_vertex()
    vd2 = g.add_vertex()
    vd3 = g.add_vertex()
    vd4 = g.add_vertex()
    vd5 = g.add_vertex()
    vd6 = g.add_vertex()
    vd7 = g.add_vertex()

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

    result = min_cut_max_flow.boykov_kolmogorov_min_st_cut( graph = g, source = vd0, sink = vd1 )

    self.assertAlmostEqual( result.weight, 7, 7 )
    ( source, sink ) = result.cutsets

    self.assertTrue( vd0 in source )
    self.assertTrue( vd1 in sink )
    self.assertEqual( len( source ), 6 )
    self.assertEqual( len( sink ), 2 )
    self.assertTrue( source, set( ( vd0, vd2, vd3, vd4, vd6, vd7 ) ) )
    self.assertTrue( sink, set( ( vd1, vd5 ) ) )


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

suite_stoer_wagner_min_cut = unittest.TestLoader().loadTestsFromTestCase(
  TestStoerWagnerMinCut
  )
suite_boykov_kolmogorov_max_flow = unittest.TestLoader().loadTestsFromTestCase(
  TestBoykovKolmogorovMaxFlow
  )
suite_boykov_kolmogorov_min_st_cut = unittest.TestLoader().loadTestsFromTestCase(
  TestBoykovKolmogorovMinSTCut
  )


alltests = unittest.TestSuite(
  [
    suite_stoer_wagner_min_cut,
    suite_boykov_kolmogorov_max_flow,
    suite_boykov_kolmogorov_min_st_cut,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

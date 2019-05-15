from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import clustering_algorithm

import unittest
from six.moves import range


class TestBetweennessCentralityClustering(unittest.TestCase):

  def graph_build(self, g):

    vds = [ g.add_vertex() for i in range( 6 ) ]
    return (
      vds,
      (
        g.add_edge( vertex1 = vds[0], vertex2 = vds[1], weight = 1 )[0],
        g.add_edge( vertex1 = vds[0], vertex2 = vds[2], weight = 1 )[0],
        g.add_edge( vertex1 = vds[1], vertex2 = vds[2], weight = 1 )[0],
        g.add_edge( vertex1 = vds[2], vertex2 = vds[3], weight = 1 )[0],
        g.add_edge( vertex1 = vds[4], vertex2 = vds[5], weight = 1 )[0],
        g.add_edge( vertex1 = vds[4], vertex2 = vds[3], weight = 1 )[0],
        g.add_edge( vertex1 = vds[3], vertex2 = vds[5], weight = 1 )[0],
        )
      )


  def undirected_params(self):

    return ( 9, [ 1, 1, 1, 1, 1, 1 ], [ [ 0, 1, 2 ], [ 3, 4, 5 ] ] )


  def directed_params(self):

    return ( 6, [ 1, 1, 1, 1, 1, 1 ], None )


  def build_and_test(self, g, params):

    ( threshold, exp_ecs, exp_comps ) = params
    ( vds, eds ) = self.graph_build( g )
    ecmap = clustering_algorithm.betweenness_centrality_clustering(
      graph = g,
      threshold = threshold,
      )

    self.assertTrue( eds[3] not in ecmap )
    self.assertEqual( len( ecmap ), len( exp_ecs ) )
    self.assertEqual( exp_ecs, [ ecmap[ ed ] for ed in eds if ed != eds[3] ] )

    if exp_comps is not None:
      from boost_adaptbx.graph import connected_component_algorithm as cca
      comps = cca.connected_components( graph = g )
      self.assertEqual(
        set( frozenset( c ) for c in comps ),
        set( frozenset( vds[i] for i in c ) for c in exp_comps )
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
      self.build_and_test( g, self.undirected_params() )


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
      self.build_and_test( g, self.undirected_params() )


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
      self.build_and_test( g, self.undirected_params() )


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
      self.build_and_test( g, self.undirected_params() )


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
      self.build_and_test( g, self.directed_params() )


suite_betweenness_centrality_clustering = unittest.TestLoader().loadTestsFromTestCase(
  TestBetweennessCentralityClustering
  )


alltests = unittest.TestSuite(
  [
    suite_betweenness_centrality_clustering,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

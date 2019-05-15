from __future__ import absolute_import, division, print_function

from boost_adaptbx import graph
from boost_adaptbx.graph import metric

import unittest
from six.moves import range


class TestBrandesBetweennessCentrality(unittest.TestCase):

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


  def undirected_results(self):

    return ( [ 0, 0, 6, 6, 0, 0 ], [ 1, 4, 4, 9, 1, 4, 4 ] )


  def directed_results(self):

    return ( [ 0, 0, 4, 3, 0, 0 ], [ 1, 3, 3, 6, 1, 1, 4 ] )


  def build_and_test(self, g, results):

    ( vds, eds ) = self.graph_build( g )
    ( exp_vcs, exp_ecs ) = results
    ( vcmap, ecmap ) = metric.brandes_betweenness_centrality( graph = g )

    self.assertEqual( len( vcmap ), len( exp_vcs ) )
    self.assertTrue( all( vd in vcmap for vd in vds ) )
    self.assertEqual( exp_vcs, [ int( vcmap[ vd ] ) for vd in vds ] )

    self.assertEqual( len( ecmap ), len( exp_ecs ) )
    self.assertTrue( all( ed in ecmap for ed in eds ) )
    self.assertEqual( exp_ecs, [ int( ecmap[ ed ] ) for ed in eds ] )


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
      self.build_and_test( g, self.undirected_results() )


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
      self.build_and_test( g, self.undirected_results() )


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
      self.build_and_test( g, self.undirected_results() )


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
      self.build_and_test( g, self.undirected_results() )


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
      self.build_and_test( g, self.directed_results() )


suite_brandes_betweenness_centrality = unittest.TestLoader().loadTestsFromTestCase(
  TestBrandesBetweennessCentrality
  )


alltests = unittest.TestSuite(
  [
    suite_brandes_betweenness_centrality,
    ]
  )


def load_tests(loader, tests, pattern):

  return alltests


if __name__ == "__main__":
  unittest.TextTestRunner( verbosity = 2 ).run( alltests )

from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext( "boost_adaptbx_graph_min_cut_max_flow_ext" )

def group_vertices(graph, parities):

  assert len( parities ) == graph.num_vertices()
  source = set()
  target = set()

  for ( v, p ) in zip( graph.vertices(), parities ):
    if p:
      source.add( v )

    else:
      target.add( v )

  return ( source, target )


class stoer_wagner_min_cut(object):
  """
  Object method version of stoer_wagner_min_cut
  """

  def __init__(self, graph):

    ( self.weight, self.parities ) = ext.stoer_wagner_min_cut( graph = graph )
    self.graph = graph


  @property
  def cutsets(self):

    return group_vertices( graph = self.graph, parities = self.parities )


class boykov_kolmogorov_max_flow(object):
  """
  Object method version of boykov_kolmogorov_max_flow
  """

  def __init__(self, graph, reverse_edge_map, source, sink):

    ( self.weight, self.parities ) = ext.boykov_kolmogorov_max_flow(
      graph = graph,
      reverse_edge_map = reverse_edge_map,
      source = source,
      sink = sink,
      )
    self.graph = graph
    self.source = source
    self.sink = sink


  @property
  def cutsets(self):

    ( source, sink ) = group_vertices( graph = self.graph, parities = self.parities )

    assert self.source in source
    assert self.sink in sink

    return ( source, sink )


class boykov_kolmogorov_min_st_cut(object):
  """
  Use Boykov-Kolmogorov algorithm to calculate min st-cut for undirected graphs
  """

  def __init__(self, graph, source, sink):

    from boost_adaptbx import graph as graphmod
    dirgraph = graphmod.adjacency_list(
      graph_type = "directed",
      vertex_type = "vector",
      edge_type = "vector",
      )

    vertex_for = {}

    for vertex in graph.vertices():
      vertex_for[ vertex ] = dirgraph.add_vertex( label = vertex )

    reverse_edge_map = {}

    for edge in graph.edges():
      s = vertex_for[ graph.source( edge = edge ) ]
      t = vertex_for[ graph.target( edge = edge ) ]
      w = graph.edge_weight( edge = edge )
      ( ed, success ) = dirgraph.add_edge( vertex1 = s, vertex2 = t, weight = w )
      assert success
      ( red, success )  = dirgraph.add_edge( vertex1 = t, vertex2 = s, weight = w )
      assert success
      reverse_edge_map[ ed ] = red
      reverse_edge_map[ red ] = ed

    self.result = boykov_kolmogorov_max_flow(
      graph = dirgraph,
      reverse_edge_map = reverse_edge_map,
      source = vertex_for[ source ],
      sink = vertex_for[ sink ],
      )


  @property
  def weight(self):

    return self.result.weight


  @property
  def cutsets(self):

    ( source, sink ) = self.result.cutsets

    return (
      set( self.result.graph.vertex_label( vertex = v ) for v in source ),
      set( self.result.graph.vertex_label( vertex = v ) for v in sink ),
      )

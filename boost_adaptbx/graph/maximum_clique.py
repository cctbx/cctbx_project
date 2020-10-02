from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
from functools import reduce
ext = bp.import_ext( "boost_adaptbx_graph_maximum_clique_ext" )
import operator

rascal = ext.rascal
greedy = ext.greedy
bron_kerbosch_all_cliques = ext.bron_kerbosch_all_cliques

def compatibility_graph(
  first,
  second,
  vertex_equality = operator.eq,
  edge_equality = operator.eq,
  output_graph_vertex = "vector",
  output_graph_edge = "vector",
  ):

  from boost_adaptbx import graph
  result = graph.adjacency_list(
    graph_type = "undirected",
    vertex_type = output_graph_vertex,
    edge_type = output_graph_edge,
    )
  vertex_from = {}

  import itertools

  for ( left, right ) in itertools.product( first.vertices(), second.vertices() ):
    if vertex_equality( first.vertex_label( left ), second.vertex_label( right ) ):
      label = ( left, right )
      vertex_from[ label ] = result.add_vertex( label = label )

  edge_weight_first = dict(
    ( frozenset( [ first.source( edge = e ), first.target( edge = e ) ] ), first.edge_weight( edge = e ) )
    for e in first.edges()
    )
  edge_weight_second = dict(
    ( frozenset( [ second.source( edge = e ), second.target( edge = e ) ] ), second.edge_weight( edge = e ) )
    for e in second.edges()
    )

  for ( ( f_l, s_l ), ( f_r, s_r ) ) in itertools.combinations( vertex_from, 2 ):
    if f_l == f_r or s_l == s_r:
      continue

    id_f = frozenset( [ f_l, f_r ] )
    id_s = frozenset( [ s_l, s_r ] )
    present_f = id_f in edge_weight_first
    present_s = id_s in edge_weight_second

    if (
      present_f
      and present_s
      and edge_equality( edge_weight_first[ id_f ], edge_weight_second[ id_s ] )
      ):
      result.add_edge(
        vertex1 = vertex_from[ ( f_l, s_l ) ],
        vertex2 = vertex_from[ ( f_r, s_r ) ],
        )

    elif not present_f and not present_s:
      result.add_edge(
        vertex1 = vertex_from[ ( f_l, s_l ) ],
        vertex2 = vertex_from[ ( f_r, s_r ) ],
        )

  return result


def selected_subgraph(graph, vertices):

  subgraph = graph.__class__()
  ext.selected_subgraph( graph = graph, subgraph = subgraph, iterable = vertices )
  return subgraph


def bessonov_projection_bound(graph, partition):
  """
  Simple upper bound by projecting back to factor graphs

  Vertices of the compatibility graph has to have properties 'first' and 'second'
  the vertex descriptors of the original graphs
  """

  descriptors = reduce( operator.add, partition, [] )
  properties = [ graph.vertex_label( vertex = v ) for v in descriptors ]

  return min(
    len( set( v[0]for v in properties ) ),
    len( set( v[1] for v in properties ) ),
    )


class raymond_gradiner_willet_projection_bound(object):
  """
  Simple upper bound by projecting back to factor graphs, and matching vertices
  according to their properties

  Vertices of the compatibility graph has to have properties 'first' and 'second'
  the vertex descriptors of the original graphs. The original graphs should
  have a vertex_label property.

  This technique is possibly very efficient for chemical line graphs
  """

  def __init__(self, first, second):

    self.first = first
    self.second = second


  def vertex_classification(self, graph, vertices):

    vertices_tagged = {}

    for v in vertices:
      label = graph.vertex_label( vertex = v )
      vertices_tagged.setdefault( label, [] ).append( v )

    return vertices_tagged


  def __call__(self, graph, partition):

    descriptors = reduce( operator.add, partition, [] )
    properties = [ graph.vertex_label( vertex = v ) for v in descriptors ]

    first_tagged = self.vertex_classification(
      graph = self.first,
      vertices = set( v[0] for v in properties ),
      )
    second_tagged = self.vertex_classification(
      graph = self.second,
      vertices = set( v[1] for v in properties ),
      )

    bound = 0

    for v in set( first_tagged ).intersection( second_tagged ):
      bound += min( len( first_tagged[ v ] ), len( second_tagged[ v ] ) )

    return bound

from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "boost_adaptbx_graph_utility_ext" )

def copy_graph(graph):

  copy = graph.__class__()
  ext.copy_graph( source = graph, target = copy )
  return copy


def copy_graph_and_map_vertices(graph):

  copy = graph.__class__()
  mapping = ext.copy_graph_and_map_vertices( source = graph, target = copy )
  return ( copy, mapping )


class undirected_graph_vertex_contraction(object):

  def __init__(self, joint_vertex_label, joint_edge_weight):

    self.joint_vertex_label = joint_vertex_label
    self.joint_edge_weight = joint_edge_weight


  def __call__(self, graph, v1, v2):

    v1_edge_to = {}

    for ed in graph.out_edges( vertex = v1 ):
      target = graph.target( edge = ed )
      assert target not in v1_edge_to
      v1_edge_to[ target ] = ed

    create_edge_between = []

    for ed in graph.out_edges( vertex = v2 ):
      target = graph.target( edge = ed )

      if target == v1 or target == v2:
        continue

      this_edge_weight = graph.edge_weight( edge = ed )

      if target in v1_edge_to:
        parallel_edge = v1_edge_to[ target ]
        graph.set_edge_weight(
          edge = parallel_edge,
          weight = self.joint_edge_weight(
            graph.edge_weight( edge = parallel_edge ),
            this_edge_weight,
            )
          )

      else:
        create_edge_between.append( ( target, this_edge_weight ) )

    graph.set_vertex_label(
      vertex = v1,
      label = self.joint_vertex_label(
        graph.vertex_label( vertex = v1 ),
        graph.vertex_label( vertex = v2 ),
        )
      )

    for ( other, weight ) in create_edge_between:
      ( ed, success ) = graph.add_edge( vertex1 = v1, vertex2 = other, weight = weight )
      assert success

    graph.remove_vertex( vertex = v2 )

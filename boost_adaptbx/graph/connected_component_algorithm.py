from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "boost_adaptbx_graph_connected_component_algorithm_ext" )


def connected_components(graph):

  result = {}

  for ( desc, component ) in ext.connected_components( graph = graph ):
    result.setdefault( component, [] ).append( desc )

  return list(result.values())

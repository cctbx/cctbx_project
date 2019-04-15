from __future__ import division

import boost.python
ext = boost.python.import_ext( "boost_adaptbx_graph_connected_component_algorithm_ext" )


def connected_components(graph):

  result = {}

  for ( desc, component ) in ext.connected_components( graph = graph ):
    result.setdefault( component, [] ).append( desc )

  return list(result.values())


from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
ext = bp.import_ext( "boost_adaptbx_graph_ext" )
from boost_adaptbx_graph_ext import *

_abbreviation_for_component = {
  "set": "set",
  "list": "list",
  "vector": "vect",
  }

_abbreviation_for_type = {
  "undirected": "undir",
  "directed": "dir",
  }


def adjacency_list(
  graph_type = "undirected",
  vertex_type = "vector",
  edge_type = "set"
  ):

  if graph_type not in _abbreviation_for_type:
    raise ValueError("Unknown graph_type: '%s'" % graph_type)

  if vertex_type not in _abbreviation_for_component:
    raise ValueError("Unknown vertex_type: '%s'" % vertex_type)

  if edge_type not in _abbreviation_for_component:
    raise ValueError("Unknown edge_type: '%s'" % edge_type)

  typename = "graph_adjlist_%s_%s_%s" % (
    _abbreviation_for_type[ graph_type ],
    _abbreviation_for_component[ vertex_type ],
    _abbreviation_for_component[ edge_type ],
    )

  try:
    return getattr( ext, typename )()

  except AttributeError:
    raise NotImplementedError

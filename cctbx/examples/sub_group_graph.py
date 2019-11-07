"""
Construct all subgroup graphs and their relations between them from a single space group.
"""
from __future__ import absolute_import, division, print_function


from cctbx import sgtbx
from cctbx.sgtbx import pointgroup_tools
import sys

def run(sg1):
  sg_high = sgtbx.space_group_info( sg1  ).group()
  sg_low  = sgtbx.space_group_info( "p1" ).group()
  graph_object =  pointgroup_tools.point_group_graph( sg_low, sg_high, False,True)
  graph_object.graph.show()
  out = """ """
  out += "digraph f { "
  out += "rankdir=LR"
  for pg in graph_object.graph.node_objects:
    for next_pg in graph_object.graph.edge_objects[ pg ]:
      pg = pg.replace( "\"","''" )
      next_pg = next_pg.replace( "\"","''" )
      out += "\""+pg+"\" -> \""+next_pg+"\" ;"
  out += "}"
  cmnd = """dot -Tpng > sg_graph.png << EOF
%s
EOF
  """%(out)
  print()
  print("Command for dot (graphviz package) to show relations between groups: ")
  print(cmnd)

if __name__=="__main__":
  run( sys.argv[1] )

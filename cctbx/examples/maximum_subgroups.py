"""
Construct all subgroup graphs and their relations between them from a single space group.
"""


from cctbx import sgtbx
from cctbx.sgtbx import show_cosets
from cctbx.sgtbx import pointgroup_tools
import sys,os


def reverse_dict( dict ):
  new_dict = {}
  for item in dict:
    for value in dict[item]:
      if value is not None:
        if new_dict.has_key( value ):
          tmp = new_dict[ value ]
          tmp.append( item )
          new_dict.update( {value:tmp} )
        else:
          new_dict.update( {value:[item]} )
  return new_dict

def get_maximum_subgroup( sg_name, reverse_graph ):
  subgroups = reverse_graph[ sg_name ]
  maximal = {}
  for sg in subgroups:
    maximal.update( {sg:True} )
  result = []
  for trial_sg in subgroups:
    tmp = reverse_graph[ trial_sg ]
    is_trial_sg_a_subgroup_of_items_in_subgroups=False
    for item in tmp:
      if item in subgroups:
        maximal.update( {item:False} )
        is_trial_sg_a_subgroup_of_subgroups=True
  for item in maximal:
    if maximal[item]:
      result.append( item )
  return result





def create_all_subgroups( sg1 ):
  sg_high = sgtbx.space_group_info( sg1  ).group()
  sg_low  = sgtbx.space_group_info( "p1" ).group()
  graph_object =  pointgroup_tools.point_group_graph( sg_low, sg_high, False,True)
  highest_sg = str( sgtbx.space_group_info( sg1  ) )
  rev_dict = reverse_dict( graph_object.graph.o )
  maximum_subgroups = get_maximum_subgroup( highest_sg, rev_dict )

  print "Subgroups of input space groups which can be constructed by introducing one single operator (and group completion) in the subgroup:"
  for sg in rev_dict[ highest_sg ]:
    line = "       "
    line += sg+(30-len(sg))*" "+str(graph_object.graph.edge_objects[ sg ][highest_sg])+(90-len( str(graph_object.graph.edge_objects[ sg ][highest_sg]) ))*" "
    print line

  print
  print "Maximum subgroup detected in the full sub-group-graph: "
  for sg in maximum_subgroups:
    line = "       "
    line += sg
    print line

  print
  print
  print
  print " Cosets for each maximum sub-group and the input space group are listed:"
  for sg in maximum_subgroups:
    print "-----------------------------------------------------------------"
    show_cosets.run( sg,highest_sg )
    print "-----------------------------------------------------------------"
    print
    print
    print
    print





def run(sg1):
  create_all_subgroups( sg1 )


if __name__=="__main__":
  run( sys.argv[1] )

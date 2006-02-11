from __future__ import generators

from cctbx.array_family import flex
from cctbx import uctbx
from cctbx import adptbx
from cctbx import sgtbx
from cctbx.sgtbx import cosets
from cctbx import maptbx
from cctbx import crystal
from cctbx import miller
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from libtbx.utils import Sorry
from libtbx import table_utils
from scitbx.python_utils import graph_tools

import sys,os

class sub_super_point_group_relations(object):
  def __init__(self,
               sg_low,
               sg_high):
    assert ( sg_low == sg_low.build_derived_point_group() )
    assert ( sg_high == sg_high.build_derived_point_group() )
    self.sg_low = sg_low
    self.sg_high = sg_high

    self.symops=[]
    self.grouped_symops = []
    self.grouped_index = []
    self.sg_groups = []
    self.left_over_symops= []

    self.get_symops_from_supergroup_that_are_not_in_subgroup(
      sg_high)

    self.assemble_symops()
    self.find_left_over_symops()

  def get_symops_from_supergroup_that_are_not_in_subgroup(self,sg_high):
    t_den = self.sg_low.t_den()
    r_den = self.sg_low.r_den()

    coset = cosets.left_decomposition(g=sg_high,
                                      h=self.sg_low)
    for set in coset.partitions[1:]:
      if (set[0].r().determinant()>0):
        self.symops.append( set[0] )

  def assemble_symops(self):
    t_den = self.sg_low.t_den()
    r_den = self.sg_low.r_den()

    # loop over all symops
    for item in range(len(self.symops)):
      symop = self.symops[item]

      tmp_symop = []
      tmp_index = []

      # multiply in the symop in the space group
      check_sg = sgtbx.space_group(self.sg_low)

      check_sg.expand_smx( symop.new_denominators(r_den, t_den) )

      # make sure that the new spacegroup is different from the old one
      # otherwise one might cycle in eternal loops in later applications
      assert check_sg.order_z() > self.sg_low.order_z()

      # Check if this SG is allready in the list
      if not ( check_sg in self.sg_groups):
        # add sg to list
        self.sg_groups.append( check_sg )
        tmp_symop.append( symop )
        tmp_index.append( item )

        # check if the other symops generate the same sg please
        for check_ops,check_item in zip(self.symops,
                                        range(len(self.symops))):
          check_again_sg = sgtbx.space_group( check_sg )
          check_again_sg.expand_smx(check_ops.new_denominators(r_den,t_den))
          if (check_again_sg == check_sg):
            # add symop to list if it is not in there yet
            if not (check_ops in tmp_symop):
              tmp_symop.append( check_ops )
            # add index to list
            if not (check_item in tmp_index):
              tmp_index.append( check_item )


        self.grouped_symops.append( tmp_symop )
        self.grouped_index.append( tmp_index )

  def find_left_over_symops(self):
    # this function gives the left over symops after
    # assuming a certain supergroup
    for set, group in zip( self.grouped_symops,
                           self.sg_groups ):
      if len(set)>0:
        select = []
        for item in range(len(self.symops)):
          symop = self.symops[item]
          if not(symop in set):
            select.append( symop )
        self.left_over_symops.append( select)
      else:
        self.left_over_symops.append( [] )

  def return_next_set(self):
    for set in self.grouped_symops:
      yield set

  def return_next_index(self):
    for iset in self.grouped_index:
      yield iset

  def return_next_sg(self):
    for sg in self.sg_groups:
      yield sg

  def return_next_left_over_set(self):
    for missing_set in self.left_over_symops:
      yield missing_set


  def show(self, out=None):
    if out is None:
      out = sys.stdout

    print >> out, "Input subgroup      : %s"%(
        sgtbx.space_group_info( group=self.sg_low ) )
    print >> out, "Input lattice group : %s"%(
      sgtbx.space_group_info( group=self.sg_high ) )
    print >> out
    print >> out
    for set,group,leftover in zip(self.grouped_symops,
                                  self.sg_groups,
                                  self.left_over_symops):
      assert( len(set)+len(leftover)==len(self.symops) )
      print >> out,   "Supergroup : %s"%( sgtbx.space_group_info( group=group ) )
      print >> out,   "             Used symops:"
      for symop in set:
        print >> out, "             (%s)  "%( symop.r().as_hkl() )
      print
      print >> out,   "             Left over symops:"
      for symop in leftover:
        if symop is not None:
          print >> out, "             (%s)  "%( symop.r().as_hkl() )
        else:
          print >> out, "             None"
      print
      print




class edge_object(object):
  def __init__(self, used, unused):
    # This object characterises a spacegroup transformation
    # by listing: used symops, unused symops
    self.symops_used = used
    self.symops_unused = unused

  def return_used(self):
    for symop in self.symops_used:
      yield symop

  def return_unused(self):
    for symop in self.symops_unused:
      yield symop

  def __repr__(self):
    repr = str()
    repr += "  using: "
    for symop in self.symops_used:
      repr+="("+symop.r().as_hkl()+")  "
    repr +="  symops left: " +str( len(self.symops_unused) )
    return repr

  def __str__(self):
    repr = str()
    repr += "  using: "
    for symop in self.symops_used:
      repr+="("+symop.r().as_hkl()+")  "
    repr +="  symops left: " +str( len(self.symops_unused) )
    return repr


class point_group_graph(object):
  def __init__(self,
               sg_low,
               sg_high):

    # It is rather import (i think) to make sure
    # that point groups are supplied. This might prevent later surprises.
    # hopfully.

    low_point_group_check = (
      sg_low ==
      sg_low.build_derived_point_group() )
    if not low_point_group_check:
      raise Sorry("Input spacegroup not a point group")

    high_point_group_check = (
      sg_high ==
      sg_high.build_derived_point_group() )
    if not high_point_group_check:
      raise Sorry("Input spacegroup not a point group")


    self.graph = graph_tools.graph()
    self.sg_low = sg_low
    self.sg_high = sg_high

    self.queue = [] # the queue used in building the space_group
    self.build_it()
    del self.queue # you have no business looking at this object, so I delete it .. ;-)

    self.graph.assert_is_clean()

  def build_it(self):
    # we start by putting the spacegroup on the queue
    self.queue.append( self.sg_low )

    while len(self.queue) > 0 :
      this_sg = self.queue.pop( 0 )
      self.make_and_place_nodes_and_connections( this_sg )

  def make_and_place_nodes_and_connections(self, input_sg):
    # make the object and name please
    object = sgtbx.space_group_info( group=input_sg )
    name = str(object)
    sg_relations = sub_super_point_group_relations(
      input_sg,
      self.sg_high)

    # loop over the possible outgoing edges
    edge_list = {}
    for possible_super_sg, used_symops, unused_symops \
        in zip( sg_relations.return_next_sg(),
                sg_relations.return_next_set(),
                sg_relations.return_next_left_over_set() ):

      # This is enough info to make connections from the given node
      edge = edge_object( used =used_symops,
                          unused = unused_symops)
      edge_list.update(
        { str(sgtbx.space_group_info(group=possible_super_sg)) : edge }
      )

      # place the sg's generated on the queue
      if not possible_super_sg in self.queue:
        self.queue.append( possible_super_sg )

    # place/insert the node with the proper connections please
    self.graph.insert_node( name = name,
                            edge_object = edge_list,
                            node_object = object )

  def remove_point_group_and_its_super_groups_from_graph(self, group ):
    if not( group == group.build_derived_point_group() ):
      raise Sorry( "Provided group is not a point group" )

    # please find all super groups of the given group
    # and remove them from the graph
    #
    # I think the easiest way is just to find all nodes with
    # no outgoing nodes, and determine all possible paths between them
    # Then group the found pg's and remove them from the graph
    end_nodes = []

    for trial_node in self.graph.edge_objects:
      n_out = len( self.graph.edge_objects[ trial_node ] )
      if n_out == 0:
        end_nodes.append( trial_node )

    # now see if there is a path between the given group and the end points
    to_be_removed = []
    group_name = str(sgtbx.space_group_info(group=group))

    for trial_end_node in end_nodes:
      tmp_paths = self.graph.find_all_paths( group_name, trial_end_node )
      for path in tmp_paths:
        for sg in path:
          if not (sg in to_be_removed):
            to_be_removed.append( sg )
    for sg in to_be_removed:
      self.graph.remove_node( sg )
    self.graph.assert_is_clean()




class find_compatible_space_groups(object):
  def __init__(self,
               likely_point_group,
               xtal_sg,
               unit_cell,
               chiral=True,
               sys_abs_rule=True):
    #
    #
    # It is assumed that the point group supplied is comming
    # from the graph outlined above
    #
    # It is also assumed that the cell provided is NOT niggli reduced
    #
    #
    #

    self.likely_point_group = likely_point_group

    self.x_sg = xtal_sg
    self.x_uc = unit_cell
    self.xs = crystal.symmetry( unit_cell=self.x_uc,
                                space_group=self.x_sg )
    self.cb_op_xs_to_nigli = self.xs.change_of_basis_op_to_minimum_cell()
    self.cb_op_l_pg_to_ref = sgtbx.space_group_info(
      group = self.likely_point_group).type().cb_op()

    self.xs_new = None

    self.point_group_compatible_sg = []
    self.is_chiral = []
    self.get_space_groups_compatible_with_likely_point_group( )

    self.get_new_xtal_symmetry()

    self.xs_new_abs_list = self.make_sys_abs_list( self.xs_new.space_group() )

    self.lpg_sg_candidates_sys_abs_list = []
    self.xs_abs_in_sg_from_lpg_abs = []
    self.sg_from_lpg_abs_in_xs_abs = []
    self.allowed = []
    self.likely_sg = []

    for sg_trial, is_chiral in zip(self.point_group_compatible_sg,
                                   self.is_chiral):

      trial_abs_list = self.make_sys_abs_list( sg_trial )
      self.lpg_sg_candidates_sys_abs_list.append( trial_abs_list )
      xs_in_sg_lpg = self.compare_abs_lists( self.xs_new_abs_list,
                                             trial_abs_list )
      sg_lpg_in_xs = self.compare_abs_lists( trial_abs_list,
                                             self.xs_new_abs_list)
      tmp_allowed = False

      if sys_abs_rule:
        if xs_in_sg_lpg:
          tmp_allowed=True

      self.allowed.append( tmp_allowed )
      if tmp_allowed:
        self.likely_sg.append( sg_trial )


  def get_space_groups_compatible_with_likely_point_group(self):
    # loop over all sg's
    for space_group_number in xrange(1,231):
      trial_group = sgtbx.space_group_info( space_group_number ).group()
      if trial_group.build_derived_point_group() == self.likely_point_group.\
            change_basis( self.cb_op_l_pg_to_ref)  :
        self.is_chiral.append( trial_group.is_chiral() )
        self.point_group_compatible_sg.append( trial_group )

  def get_new_xtal_symmetry(self):
    # apply first the nigli operater etc
    self.xs_new = self.xs.change_basis( self.cb_op_xs_to_nigli )
    self.xs_new = self.xs_new.change_basis( self.cb_op_l_pg_to_ref )

  def make_sys_abs_list(self, space_group):
    max_index=(6,6,6)
    sg_type = sgtbx.space_group_info( 1 ).type()
    tmp_miller_list = miller.index_generator( sg_type,
                                              False,
                                              max_index ).to_array()
    abs_list = flex.miller_index()
    for miller_index in tmp_miller_list:
      if space_group.is_sys_absent( miller_index ):
        abs_list.append( miller_index )
    return abs_list

  def compare_abs_lists(self,
                        abs_list_sub,
                        abs_list_super):
    # check sysb abs
    contains=True
    for sub_item in abs_list_sub:
      if not (sub_item in abs_list_super):
        contains=False
    return contains

  def show(self):
    print "Input space group  : ", sgtbx.space_group_info(group=self.x_sg)
    print "Input unit cell    : ", self.x_uc.parameters()
    print "Likely point group : ", sgtbx.space_group_info(group=self.likely_point_group)
    print
    print "New unit cell      : ", self.xs_new.unit_cell().parameters()
    print "Possible space groups "
    for sg in self.likely_sg:
      print sgtbx.space_group_info( group=sg )

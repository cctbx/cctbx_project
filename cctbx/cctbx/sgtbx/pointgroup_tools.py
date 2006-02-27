from __future__ import generators

from cctbx.array_family import flex
from cctbx import uctbx
from cctbx import adptbx
from cctbx import sgtbx
from cctbx.crystal.find_best_cell import alternative_find_best_cell
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


def reference_setting_choices(space_group):
  cyclic_permutations =  ['x,y,z',
                          'y,z,x',
                          'z,x,y' ]

  adams_group = sgtbx.space_group_info(
    group=space_group.build_derived_group(False,False) )

  space_group = sgtbx.space_group_info(group=space_group)

  # please check that we have something in reference setting
  # just to make sure that thne thing is used for it's original purpose
  assert space_group.is_reference_setting()

  info = []

  identity_op = sgtbx.change_of_basis_op( 'x,y,z' ).c().r()

  for cyclic_permutation in cyclic_permutations:

    cob_op = sgtbx.change_of_basis_op( cyclic_permutation )
    transformed_adams_group = adams_group.change_basis( cob_op )
    transformed_space_group = space_group.change_basis( cob_op )

    cob_to_ref_sg = transformed_space_group.\
                    change_of_basis_op_to_reference_setting()
    cob_to_ref_pg = transformed_adams_group.\
                    change_of_basis_op_to_reference_setting()


    adams_norm = False
    space_norm = False

    # check if the rotation part of the cb_op to ref is
    # the identity operator

    # if hall symbols are equal, sg's are equal
    if ( identity_op == cob_to_ref_pg.c().r() ):
      adams_norm=True

    if ( identity_op == cob_to_ref_sg.c().r() ):
      space_norm=True

    info_tuple = (cob_op, cob_to_ref_sg, adams_norm, space_norm)
    info.append( info_tuple )

  possible_additional_transforms = []
  # we have to of course take into account the identity operator
  possible_additional_transforms.append( info[0][0]*info[0][1] )
  for ii in info:
    if ii[2]: # should fall in the adams normalizer
      if not ii[3]: # should NOT fall in the space normalizer
        # cob should ONLY be applied on unit cell, not to the sg.
        possible_additional_transforms.append( ii[0] )

  return possible_additional_transforms


def coset_lookup(pg_low,
                 pg_high):
  coset = cosets.left_decomposition(g=pg_high,
                                    h=pg_low)
  full_cosets = []

  for set in coset.partitions:
    if (set[0].r().determinant()>0):
      tmp = []
      for symop in set:
        tmp.append( symop.r().as_hkl() )
      full_cosets.append( tmp )
  return full_cosets


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
               pg_low,
               pg_high):

    # It is rather import (i think) to make sure
    # that point groups are supplied. This might prevent later surprises.
    # hopefully.

    low_point_group_check = (
      pg_low ==
      pg_low.build_derived_point_group() )
    if not low_point_group_check:
      raise Sorry("Input spacegroup not a point group")

    high_point_group_check = (
      pg_high ==
      pg_high.build_derived_point_group() )
    if not high_point_group_check:
      raise Sorry("Input spacegroup not a point group")


    self.graph = graph_tools.graph()
    self.pg_low = pg_low
    self.pg_high = pg_high

    self.queue = [] # the queue used in building the space_group
    self.build_it()
    del self.queue # you have no business looking at this object, so I delete it .. ;-)

    self.graph.assert_is_clean()

  def build_it(self):
    # we start by putting the spacegroup on the queue
    self.queue.append( self.pg_low )

    while len(self.queue) > 0 :
      this_sg = self.queue[ 0 ]
      self.queue.remove( this_sg )

      self.make_and_place_nodes_and_connections( this_sg )

  def make_and_place_nodes_and_connections(self, input_sg):
    # make the object and name please
    object = sgtbx.space_group_info( group=input_sg )
    name = str(object)
    sg_relations = sub_super_point_group_relations(
      input_sg,
      self.pg_high)

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
    # print object, type( object)
    self.graph.insert_node( name = name,
                            edge_object = edge_list,
                            node_object = object )

  def remove_point_group_and_its_super_groups_from_graph(self, group_name ):
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
               likely_pointgroup,
               xtal_sg,
               unit_cell,
               sys_abs_flag=True):

    self.x_sg = xtal_sg
    self.x_uc = unit_cell

    self.xs = crystal.symmetry( self.x_uc,
                                space_group=self.x_sg)

    self.cb_op_xs_to_niggli = self.xs.change_of_basis_op_to_niggli_cell()

    self.x_lpg = likely_pointgroup
    self.cb_op_lpg_to_ref_set = sgtbx.space_group_info(
      group=self.x_lpg).change_of_basis_op_to_reference_setting()

    self.point_group_compatible_sg = []
    self.is_chiral = []
    self.get_space_groups_compatible_with_likely_point_group()

    self.allowed_under_pg_and_sys_abs = []


    for sg in self.point_group_compatible_sg:
      additional_cb_ops = reference_setting_choices( sg )
      for add_cb_op in additional_cb_ops: # not optimal, but EASY
         # make a new xs object please
        new_xs = self.xs.change_basis(add_cb_op *
                                      self.cb_op_lpg_to_ref_set *
                                      self.cb_op_xs_to_niggli )
        if sys_abs_flag:

          # get sys abs for this setting please
          current_xtal_sys_abs = self.make_sys_abs_list( new_xs.space_group() )
          # get sys abs for current sg please
          trial_sg_abs = self.make_sys_abs_list( sg )
          # check if current abs are 'in' trial abs
          inside = self.compare_abs_lists( current_xtal_sys_abs,
                                           trial_sg_abs )
          if inside:
            final_xs = crystal.symmetry( unit_cell = new_xs.unit_cell(),
                                         space_group = sg,
                                         assert_is_compatible_unit_cell=False )
            best_cell = alternative_find_best_cell(
              final_xs.unit_cell(),
              final_xs.space_group() )
            final_xs = best_cell.return_best_xs()
            # please get the total cb_op from input cell to best cell
            tmp_cb_op = best_cell.return_change_of_basis_op_to_best_cell()
            final_cb_op = ( tmp_cb_op*
                            add_cb_op*
                            self.cb_op_lpg_to_ref_set*
                            self.cb_op_xs_to_niggli)
            self.allowed_under_pg_and_sys_abs.append( (final_xs,
                                                        final_cb_op) )


  def get_space_groups_compatible_with_likely_point_group(self):
    # loop over all standard sg's
    for space_group_number in xrange(1,231):
      trial_group = sgtbx.space_group_info( space_group_number ).group()

      if trial_group.build_derived_group(False,False) \
             == self.x_lpg.change_basis( self.cb_op_lpg_to_ref_set ):
        self.is_chiral.append( trial_group.is_chiral() )
        self.point_group_compatible_sg.append( trial_group )

  def make_sys_abs_list(self, space_group):
    max_index=(5,5,5)
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

  def show(self, out=None):
    if out == None:
      out=sys.stdout

    print >> out, "Input space group  : ", sgtbx.space_group_info(group=self.x_sg)
    print >> out, "Input unit cell    : ", self.x_uc.parameters()
    print >> out, "Likely point group : ", sgtbx.space_group_info(group=self.lpg)
    print >> out
    print >> out, "Possible crystal symmetries: "
    for sx_cb in self.allowed_under_pg_and_sys_abs:
      sx_cb.show(out)
      print >> out


class node_object(object):
  def __init__(self,
               point_group_info,
               allowed_xtal_syms):
    self.point_group_info = point_group_info
    self.allowed_xtal_syms = allowed_xtal_syms


class space_group_graph_from_cell_and_sg(object):
  def __init__(self,
               unit_cell,
               sg_low,
               max_delta=5.0):
    # No primitive settings assumed
    self.unit_cell = unit_cell
    self.sg_low = sg_low
    self.pg_low = sg_low.build_derived_point_group()
    self.xs_low = crystal.symmetry( unit_cell = self.unit_cell,
                                    space_group=self.sg_low)

    self.pg_low_prim_set = self.pg_low.change_basis(
      self.xs_low.change_of_basis_op_to_niggli_cell() )

    self.xs_prim_set = self.xs_low.change_basis(
      self.xs_low.change_of_basis_op_to_niggli_cell() )

    self.pg_high = sgtbx.lattice_symmetry.group(
      self.xs_prim_set.unit_cell(),
      max_delta=max_delta)


    self.pg_graph = point_group_graph( self.pg_low_prim_set,
                                       self.pg_high )

    self.coset_table = coset_lookup( self.pg_low_prim_set,
                                     self.pg_high )

    new_dict = {}
    for pg_name in self.pg_graph.graph.node_objects:
      # for each possible point group, find a set of allowed space_groups
      pg_object = self.pg_graph.graph.node_objects[ pg_name ].group()

      # now find out the possible space groups please
      tmp_allowed_sgs = find_compatible_space_groups(
        pg_object, # point group
        self.sg_low, # sg of given system
        self.unit_cell) # uc of given system

      new_node = node_object(
        self.pg_graph.graph.node_objects[ pg_name ],
        tmp_allowed_sgs.allowed_under_pg_and_sys_abs)
      # now it is good to replace the node with the one we just made
      new_dict.update( {pg_name:
                        new_node} )
    #
    self.pg_graph.graph.node_objects = new_dict


  def return_likely_sg_and_cell(self):
    for pg in self.pg_graph.graph.node_objects:
      object = self.pg_graph.graph.node_objects[ pg ]
      for sg in object.allowed_xtal_syms:
        yield ( sg[0].space_group(), sg[0].unit_cell() )

  def show(self, out=None):

    if out==None:
      out=sys.stdout
    print >> out, "------------------------"
    print >> out, "Vertices and their edges"
    print >> out, "------------------------"
    print >> out
    for pg in self.pg_graph.graph.node_objects:
      print >> out, "Point group  ", pg, "  is a maximal subgroup of :"
      if (len( self.pg_graph.graph.o[ str( pg ) ] )==0):
         print >> out, "  * None"
      else:
        for edge in self.pg_graph.graph.o[ str( pg ) ]:
          print >> out, "  *", edge
      print >> out

    print >> out
    print >> out
    print >> out, "-------------------------"
    print >> out, "Transforming point groups"
    print >> out, "-------------------------"
    print >> out
    for pg in self.pg_graph.graph.node_objects:
      for next_pg in self.pg_graph.graph.edge_objects[ pg ]:
        print >> out, "From", pg, "  to ", next_pg, " using :"
        for symop in self.pg_graph.graph.edge_objects[ pg ][ next_pg].\
                return_used():
          print >> out, "  * ",symop.r().as_hkl()
        print >> out

    print >> out
    print >> out
    print >> out, "----------------------"
    print >> out, "Compatible spacegroups"
    print >> out, "----------------------"
    print >> out
    print >> out, "Spacegroups compatible with a specified point group "
    print >> out, "**and** with the systematic absenses specified by the "
    print >> out, "input space group, are listed below."
    print >> out

    for pg in self.pg_graph.graph.node_objects:
      print >> out, "Spacegroup candidates in point group %s:"%(pg)
      for trial_sym in self.pg_graph.graph.node_objects[ pg ].allowed_xtal_syms:
        trial_sg = trial_sym[0].space_group_info()
        print >> out, "  *", trial_sg,
        print >> out, " %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f"\
            %(trial_sym[0].unit_cell().parameters()[0],
              trial_sym[0].unit_cell().parameters()[1],
              trial_sym[0].unit_cell().parameters()[2],
              trial_sym[0].unit_cell().parameters()[3],
              trial_sym[0].unit_cell().parameters()[4],
              trial_sym[0].unit_cell().parameters()[5])
      print >> out


  def graphviz_pg_graph(self, out=None ):
    if out==None:
      out=sys.stdout
    print >> out, "digraph f { "
    print >> out, "rankdir=LR"
    for pg in self.pg_graph.graph.node_objects:
      for next_pg in self.pg_graph.graph.edge_objects[ pg ]:
        pg = pg.replace( "\"","''" )
        next_pg = next_pg.replace( "\"","''" )

        print >> out, "\"",pg,"\" -> \"", next_pg, "\" ;"
    print >> out, "}"

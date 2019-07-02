from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx.crystal.find_best_cell import alternative_find_best_cell
from cctbx.sgtbx import cosets
from cctbx import crystal
from cctbx import miller
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from scitbx import matrix
from scitbx.python_utils import graph_tools
from libtbx.utils import Sorry
import sys
from six.moves import range
from six.moves import zip


def reference_setting_choices(space_group):
  # we used to have
  cyclic_permutations =  ['x,y,z',
                          'y,z,x',
                          'z,x,y' ]




  adams_group = sgtbx.space_group_info(
    group=space_group.build_derived_group(False, False))

  space_group = sgtbx.space_group_info(group=space_group)

  # please check that we have something in reference setting
  # just to make sure that thne thing is used for it's original purpose
  assert space_group.is_reference_setting()

  info = []

  identity_op = sgtbx.change_of_basis_op('x,y,z').c().r()

  for cyclic_permutation in cyclic_permutations:

    cob_op = sgtbx.change_of_basis_op(cyclic_permutation)
    transformed_adams_group = adams_group.change_basis(cob_op)
    transformed_space_group = space_group.change_basis(cob_op)

    cob_to_ref_sg = transformed_space_group.\
                    change_of_basis_op_to_reference_setting()
    cob_to_ref_pg = transformed_adams_group.\
                    change_of_basis_op_to_reference_setting()


    adams_norm = False
    space_norm = False

    # check if the rotation part of the cb_op to ref is
    # the identity operator

    # if hall symbols are equal, sg's are equal
    if (identity_op == cob_to_ref_pg.c().r()):
      adams_norm=True

    if (identity_op == cob_to_ref_sg.c().r()):
      space_norm=True

    info_tuple = (cob_op, cob_to_ref_sg, adams_norm, space_norm)
    info.append(info_tuple)

  possible_additional_transforms = []
  # we have to of course take into account the identity operator
  possible_additional_transforms.append(info[0][0]*info[0][1])
  for ii in info:
    if ii[2]: # should fall in the adams normalizer
      if not ii[3]: # should NOT fall in the space normalizer
        # cob should ONLY be applied on unit cell, not to the sg.
        possible_additional_transforms.append(ii[0])

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
        tmp.append(symop.r().as_hkl())
      full_cosets.append(tmp)
  return full_cosets


class sub_super_point_group_relations(object):
  def __init__(self,
               sg_low,
               sg_high,
               enforce_point_groups=True):
    self.enforce_point_groups=enforce_point_groups
    if enforce_point_groups:
      assert (sg_low == sg_low.build_derived_point_group())
      assert (sg_high == sg_high.build_derived_point_group())
    self.sg_low = sg_low
    self.sg_high = sg_high

    self.symops=[]
    self.grouped_symops = []
    self.grouped_index = []
    self.sg_groups = []
    self.left_over_symops= []

    self.get_symops_from_supergroup_that_are_not_in_subgroup(sg_high)

    self.assemble_symops()
    self.find_left_over_symops()


  def get_symops_from_supergroup_that_are_not_in_subgroup(self, sg_high):
    coset = cosets.left_decomposition(g=sg_high,
                                      h=self.sg_low)
    for set in coset.partitions[1:]:
      if self.enforce_point_groups:
        if set[0].r().determinant()>0:
          self.symops.append(set[0])
      else:
        self.symops.append(set[0])

  def assemble_symops(self):
    t_den = self.sg_low.t_den()
    r_den = self.sg_low.r_den()

    # loop over all symops
    for item, symop in enumerate(self.symops):
      tmp_symops = []
      tmp_indices = []

      # multiply in the symop in the space group
      check_sg = sgtbx.space_group(self.sg_low)

      check_sg.expand_smx(symop.new_denominators(r_den, t_den))
      # Check if this SG is already in the list
      assert check_sg != self.sg_low
      if check_sg not in self.sg_groups:
        # add sg to list
        self.sg_groups.append(check_sg)
        tmp_symops.append(symop)
        tmp_indices.append(item)

        # check if the other symops generate the same sg please
        for check_item, check_op in enumerate(self.symops):
          if check_sg.contains(check_op.new_denominators(r_den,t_den)):
            # add symop to list if it is not in there yet
            if check_op not in tmp_symops:
              tmp_symops.append(check_op)
            # add index to list
            if check_item not in tmp_indices:
              tmp_indices.append(check_item)

        self.grouped_symops.append(tmp_symops)
        self.grouped_index.append(tmp_indices)

  def find_left_over_symops(self):
    # this function gives the left over symops after
    # assuming a certain supergroup
    for set, group in zip(self.grouped_symops,
                          self.sg_groups):
      if len(set)>0:
        select = []
        for item, symop in enumerate(self.symops):
          if symop not in set:
            select.append(symop)
        self.left_over_symops.append(select)
      else:
        self.left_over_symops.append([])

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

    print((
      "Input subgroup      : %s" % sgtbx.space_group_info(group=self.sg_low)), file=out)
    print((
      "Input lattice group : %s" % sgtbx.space_group_info(group=self.sg_high)), file=out)
    print(file=out)
    print(file=out)
    for set,group,leftover in zip(self.grouped_symops,
                                  self.sg_groups,
                                  self.left_over_symops):
      assert(len(set)+len(leftover)==len(self.symops))
      print((
        "Supergroup : %s" % sgtbx.space_group_info(group=group)), file=out)
      print("             Used symops:", file=out)
      for symop in set:
        print("             (%s)    "%(symop.r().as_hkl()), file=out)
      print(file=out)
      print("             Left over symops:", file=out)
      for symop in leftover:
        if symop is not None:
          print("             (%s)  "%(symop.r().as_hkl()), file=out)
        else:
          print("             None", file=out)
      print(file=out)
      print(file=out)




class edge_object(object):
  def __init__(self, used, unused,as_xyz=False):
    # This object characterises a spacegroup transformation
    # by listing: used symops, unused symops
    self.symops_used = used
    self.symops_unused = unused
    self.as_xyz = as_xyz

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
      if self.as_xyz:
        repr+="("+symop.as_xyz()+")  "
      else:
        repr+="("+symop.r().as_hkl()+")  "
    repr +="  symops left: " +str(len(self.symops_unused))
    return repr

  def __str__(self):
    repr = str()
    repr += "  using: "
    for symop in self.symops_used:
      if self.as_xyz:
        repr+="("+symop.as_xyz()+")  "
      else:
        repr+="("+symop.r().as_hkl()+")  "
    repr +="  symops left: " +str(len(self.symops_unused))
    return repr



class point_group_graph(object):
  def __init__(self,
               pg_low,
               pg_high,
               enforce_point_group=True,
               as_xyz = False):

    # It is rather import (i think) to make sure
    # that point groups are supplied. This might prevent later surprises.
    # hopefully.
    self.as_xyz = as_xyz

    low_point_group_check = (pg_low == pg_low.build_derived_point_group())
    if enforce_point_group:
      if not low_point_group_check:
        raise Sorry("Input spacegroup not a point group")

    high_point_group_check = (pg_high == pg_high.build_derived_point_group())

    if enforce_point_group:
      if not high_point_group_check:
        raise Sorry("Input spacegroup not a point group")

    self.assert_pg = enforce_point_group
    self.graph = graph_tools.graph()
    self.pg_low = pg_low
    self.pg_high = pg_high

    self.queue = [] # the queue used in building the space_group
    self.build_it()
    del self.queue # you have no business looking at this object,
                   # so I delete it .. ;-)

    self.graph.assert_is_clean()

  def build_it(self):
    # we start by putting the spacegroup on the queue
    self.queue.append(self.pg_low)

    while len(self.queue) > 0 :
      this_sg = self.queue.pop(0)
      self.make_and_place_nodes_and_connections(this_sg)

  def make_and_place_nodes_and_connections(self, input_sg):
    # make the object and name please
    object = sgtbx.space_group_info(group=input_sg)
    name = str(object)
    sg_relations = sub_super_point_group_relations(
      input_sg,
      self.pg_high,
      self.assert_pg)

    # loop over the possible outgoing edges
    edge_list = {}
    for possible_super_sg, used_symops, unused_symops \
        in zip(sg_relations.return_next_sg(),
               sg_relations.return_next_set(),
               sg_relations.return_next_left_over_set()):

      # This is enough info to make connections from the given node
      edge = edge_object(used =used_symops,
                         unused = unused_symops,
                         as_xyz = self.as_xyz)
      edge_list[str(sgtbx.space_group_info(group=possible_super_sg)) ] = edge

      # place the sg's generated on the queue
      if not possible_super_sg in self.queue:
        self.queue.append(possible_super_sg)

    # place/insert the node with the proper connections please
    # print object, type(object)
    self.graph.insert_node(name = name,
                          edge_object = edge_list,
                          node_object = object)

  def remove_point_group_and_its_super_groups_from_graph(self, group_name):
    # please find all super groups of the given group
    # and remove them from the graph
    #
    # I think the easiest way is just to find all nodes with
    # no outgoing nodes, and determine all possible paths between them
    # Then group the found pg's and remove them from the graph
    end_nodes = []

    for trial_node in self.graph.edge_objects:
      n_out = len(self.graph.edge_objects[trial_node])
      if n_out == 0:
        end_nodes.append(trial_node)

    # now see if there is a path between the given group and the end points
    to_be_removed = []

    for trial_end_node in end_nodes:
      tmp_paths = self.graph.find_all_paths(group_name, trial_end_node)
      for path in tmp_paths:
        for sg in path:
          if not (sg in to_be_removed):
            to_be_removed.append(sg)
    for sg in to_be_removed:
      self.graph.remove_node(sg)
    self.graph.assert_is_clean()


  def reverse_dict(self):
      new_dict = {}
      for item in self.graph.o:
        for value in self.graph.o[item]:
          if value is not None:
            if value in new_dict:
              tmp = new_dict[value]
              tmp.append(item)
              new_dict[value] = tmp
            else:
              new_dict[value] = [item]
      return new_dict

  def get_maximal_subgroup(self, sg_name):
      subgroups = []
      reverse_graph = self.reverse_dict()
      if sg_name in reverse_graph:
        subgroups = reverse_graph[sg_name]
      maximal = {}
      for sg in subgroups:
        maximal[sg] = True
      result = []
      for trial_sg in subgroups:
        tmp = {}
        if trial_sg in reverse_graph:
          tmp = reverse_graph[trial_sg]
        is_trial_sg_a_subgroup_of_items_in_subgroups=False
        for item in tmp:
          if item in subgroups:
            maximal[item] = False
            is_trial_sg_a_subgroup_of_subgroups=True
      for item in maximal:
        if maximal[item]:
          result.append(item)
      return result




class find_compatible_space_groups(object):
  def __init__(self,
               likely_pointgroup=None,
               xtal_sg=None,
               unit_cell=None,
               sys_abs_flag=True,
               miller_array=None):
    # we have the choice of supplynig either a dataset
    # or just cell, symmetry and likely point group
    # when supplynig data, an attempt will bhe made to detemine
    # the most likely spacegroup

    self.miller_array = None
    self.x_sg = None
    self.x_uc = None
    self.x_lpg = None

    if (xtal_sg is None) or (unit_cell is None):
      assert miller_array is not None
      self.miller_array = miller_array
      assert likely_pointgroup is None
      self.miller_array = miller_array
      self.x_sg = miller_array.space_group()
      self.x_uc = miller_array.unit_cell()
      self.x_lpg = miller_array.space_group().build_derived_group(True, True)

    if miller_array is None:
      assert xtal_sg is not None
      assert unit_cell is not None
      assert likely_pointgroup is not None
      self.x_lpg = likely_pointgroup
      self.x_sg = xtal_sg
      self.x_uc = unit_cell

    self.xs = crystal.symmetry(self.x_uc,
                              space_group=self.x_sg)

    self.cb_op_xs_to_niggli = self.xs.change_of_basis_op_to_niggli_cell()

    self.cb_op_lpg_to_ref_set = sgtbx.space_group_info(
      group=self.x_lpg).change_of_basis_op_to_reference_setting()

    self.point_group_compatible_sg = []
    self.is_chiral = []
    self.get_space_groups_compatible_with_likely_point_group()

    self.allowed_under_pg_and_sys_abs = []


    for sg in self.point_group_compatible_sg:
      additional_cb_ops = reference_setting_choices(sg)
      for add_cb_op in additional_cb_ops: # not optimal, but EASY
         # make a new xs object please
        new_xs = self.xs.change_basis(add_cb_op *
                                      self.cb_op_lpg_to_ref_set *
                                      self.cb_op_xs_to_niggli)
        if sys_abs_flag:

          # get sys abs for this setting please
          current_xtal_sys_abs = self.make_sys_abs_list(new_xs.space_group())
          # get sys abs for current sg please
          trial_sg_abs = self.make_sys_abs_list(sg)
          # check if current abs are 'in' trial abs
          inside = self.compare_abs_lists(current_xtal_sys_abs,
                                          trial_sg_abs)
          if inside:
            final_xs = crystal.symmetry(unit_cell = new_xs.unit_cell(),
                                        space_group = sg,
                                        assert_is_compatible_unit_cell=False)
            best_cell = alternative_find_best_cell(
              final_xs.unit_cell(),
              final_xs.space_group())
            final_xs = best_cell.return_best_xs()
            # please get the total cb_op from input cell to best cell
            tmp_cb_op = best_cell.return_change_of_basis_op_to_best_cell()
            final_cb_op = ( tmp_cb_op*
                            add_cb_op*
                            self.cb_op_lpg_to_ref_set*
                            self.cb_op_xs_to_niggli)
            self.allowed_under_pg_and_sys_abs.append( (final_xs,
                                                        final_cb_op) )
    #-------------------------------------------
    self.full_mask=None
    if self.miller_array is not None:
      self.full_mask = flex.bool(self.miller_array.data().size(), True)
      # True when present, false when absent
      for xs_and_cb_op in self.allowed_under_pg_and_sys_abs:
        xs_and_cb_op[0].show_summary()
        xs_and_cb_op[1].as_xyz()

  def find_absent_indices():
    print()

  def get_space_groups_compatible_with_likely_point_group(self):
    # loop over all standard sg's
    for space_group_number in range(1,231):
      trial_group = sgtbx.space_group_info(space_group_number).group()

      if trial_group.build_derived_group(False,False) \
             == self.x_lpg.change_basis(self.cb_op_lpg_to_ref_set):
        self.is_chiral.append(trial_group.is_chiral())
        self.point_group_compatible_sg.append(trial_group)

  def make_sys_abs_list(self, space_group):
    max_index=(5,5,5)
    sg_type = sgtbx.space_group_info(1).type()
    tmp_miller_list = miller.index_generator(sg_type,
                                             False,
                                             max_index).to_array()
    abs_list = flex.miller_index()
    for miller_index in tmp_miller_list:
      if space_group.is_sys_absent(miller_index):
        abs_list.append(miller_index)
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

  def full_sys_abs_mask(self):
    # here we make a
    print()
    print()



  def show(self, out=None):
    if out == None:
      out=sys.stdout

    print("Input space group  : ", sgtbx.space_group_info(
      group=self.x_sg), file=out)
    print("Input unit cell    : ", self.x_uc.parameters(), file=out)
    print("Likely point group : ", sgtbx.space_group_info(
      group=self.lpg), file=out)
    print(file=out)
    print("Possible crystal symmetries: ", file=out)
    for sx_cb in self.allowed_under_pg_and_sys_abs:
      sx_cb.show(out)
      print(file=out)


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
               max_delta=5.0,
               ):
    # No primitive settings assumed
    self.unit_cell = unit_cell
    self.sg_low = sg_low
    self.pg_low = sg_low.build_derived_point_group()
    self.xs_low = crystal.symmetry(unit_cell = self.unit_cell,
                                   space_group=self.sg_low)

    self.pg_low_prim_set = self.pg_low.change_basis(
      self.xs_low.change_of_basis_op_to_niggli_cell())

    self.xs_prim_set = self.xs_low.change_basis(
      self.xs_low.change_of_basis_op_to_niggli_cell())

    self.pg_high = sgtbx.lattice_symmetry.group(
      self.xs_prim_set.unit_cell(),
      max_delta=max_delta)


    self.pg_graph = point_group_graph(self.pg_low_prim_set,
                                      self.pg_high)

    self.coset_table = coset_lookup(self.pg_low_prim_set,
                                    self.pg_high)

    new_dict = {}
    for pg_name in self.pg_graph.graph.node_objects:
      # for each possible point group, find a set of allowed space_groups
      pg_object = self.pg_graph.graph.node_objects[pg_name].group()

      # now find out the possible space groups please
      tmp_allowed_sgs = find_compatible_space_groups(
        likely_pointgroup=pg_object, # point group
        xtal_sg=self.sg_low, # sg of given system
        unit_cell=self.unit_cell) # uc of given system

      new_node = node_object(
        self.pg_graph.graph.node_objects[pg_name],
        tmp_allowed_sgs.allowed_under_pg_and_sys_abs)
      # now it is good to replace the node with the one we just made
      new_dict[pg_name] = new_node
    #
    self.pg_graph.graph.node_objects = new_dict


  def return_likely_sg_and_cell(self):
    for pg in self.pg_graph.graph.node_objects:
      object = self.pg_graph.graph.node_objects[ pg ]
      for sg in object.allowed_xtal_syms:
        yield (sg[0].space_group(), sg[0].unit_cell())

  def show(self, out=None):

    if out==None:
      out=sys.stdout
    print("------------------------", file=out)
    print("Vertices and their edges", file=out)
    print("------------------------", file=out)
    print(file=out)
    for pg in self.pg_graph.graph.node_objects:
      print("Point group  ", pg, "  is a maximal subgroup of :", file=out)
      if (len(self.pg_graph.graph.o[ str(pg) ])==0):
        print("  * None", file=out)
      else:
        for edge in self.pg_graph.graph.o[ str(pg) ]:
          print("  *", edge, file=out)
      print(file=out)

    print(file=out)
    print(file=out)
    print("-------------------------", file=out)
    print("Transforming point groups", file=out)
    print("-------------------------", file=out)
    print(file=out)
    for pg in self.pg_graph.graph.node_objects:
      for next_pg in self.pg_graph.graph.edge_objects[ pg ]:
        print("From", pg, "  to ", next_pg, " using :", file=out)
        for symop in self.pg_graph.graph.edge_objects[ pg ][ next_pg]\
                                        .return_used():
          print("  * ",symop.r().as_hkl(), file=out)
        print(file=out)

    print(file=out)
    print(file=out)
    print("----------------------", file=out)
    print("Compatible spacegroups", file=out)
    print("----------------------", file=out)
    print(file=out)
    print("Spacegroups compatible with a specified point group ", file=out)
    print("**and** with the systematic absenses specified by the ", file=out)
    print("input space group, are listed below.", file=out)
    print(file=out)

    for pg in self.pg_graph.graph.node_objects:
      print("Spacegroup candidates in point group %s:"%(pg), file=out)
      for trial_sym in self.pg_graph.graph.node_objects[ pg ]\
                                    .allowed_xtal_syms:
        trial_sg = trial_sym[0].space_group_info()
        print("  *", trial_sg, end=' ', file=out)
        print(" %5.2f %5.2f %5.2f %5.2f %5.2f %5.2f"\
            %(trial_sym[0].unit_cell().parameters()[0],
              trial_sym[0].unit_cell().parameters()[1],
              trial_sym[0].unit_cell().parameters()[2],
              trial_sym[0].unit_cell().parameters()[3],
              trial_sym[0].unit_cell().parameters()[4],
              trial_sym[0].unit_cell().parameters()[5]), file=out)
      print(file=out)


  def graphviz_pg_graph(self, out=None):
    if out==None:
      out=sys.stdout
    print("digraph f { ", file=out)
    print("rankdir=LR", file=out)
    for pg in self.pg_graph.graph.node_objects:
      for next_pg in self.pg_graph.graph.edge_objects[ pg ]:
        pg = pg.replace("\"","''")
        next_pg = next_pg.replace("\"","''")

        print("\"",pg,"\" -> \"", next_pg, "\" ;", file=out)
    print("}", file=out)


def compatible_symmetries(point_group):
  """ Primitive setting assumed """
  for op in point_group:
    r = op.r()
    order = r.order()
    if r.info().type() == 1: continue
    yield op
    invariants = [ matrix.col(u) for u in r.info().basis_of_invariant() ]
    if len(invariants) == 2:
      t1, t2 = invariants
      invariants.extend((t1 + t2, t1 - t2))
    translations = []
    for t in invariants:
      t = sgtbx.tr_vec(t, order).mod_short()
      if not t.is_zero() and t not in translations:
        translations.append(t)
    for t in translations:
      yield sgtbx.rt_mx(r, t.new_denominator(sgtbx.sg_t_den))

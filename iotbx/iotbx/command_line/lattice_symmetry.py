#! /usr/bin/env python

# Comments by Phil Evans, MRC-LMB, Cambridge, U.K.

from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.sgtbx import lattice_symmetry
from cctbx.sgtbx import bravais_types
from iotbx.option_parser import iotbx_option_parser
from scitbx.array_family import flex
import math
import sys

def run():
  command_line = (iotbx_option_parser(
    usage="iotbx.lattice_symmetry [options] [centring_type_symbol]",
    description="Example: iotbx.lattice_symmetry"
               +" --unit_cell=12,12,12.1,89,90,92 F")
    .enable_symmetry_comprehensive()
    .option(None, "--delta",
      action="store",
      type="float",
      default=3.,
      dest="delta",
      help="angular tolerance in degrees")
  ).process(max_nargs=1)
  # Pick up symmetry object
  input_symmetry = command_line.symmetry
  # Check that we have what we need
  if (input_symmetry.unit_cell() is None):
    print
    print "***********************************"
    print "Please specify unit cell parameters"
    print "***********************************"
    print
    command_line.parser.show_help()
    return
  if (len(command_line.args) > 0):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="Hall: %s 1" % command_line.args[0])
  elif (input_symmetry.space_group_info() is None):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="P 1")
  # Do it
  show(input_symmetry, command_line.options.delta)

def show(input_symmetry, max_delta):
  print
  print "Input"
  print "====="
  print
  input_symmetry.show_summary()
  print
  print "Angular tolerance: %.3f degrees" % max_delta
  print
  print "Similar symmetries"
  print "=================="
  print
  # Get cell reduction operator
  cb_op_inp_minimum = input_symmetry.change_of_basis_op_to_minimum_cell()

  # New symmetry object with changed basis
  minimum_symmetry = input_symmetry.change_basis(cb_op_inp_minimum)

  # Get highest symmetry compatible with lattice
  lattice_group = lattice_symmetry.group(
    minimum_symmetry.unit_cell(), max_delta=max_delta)
  lattice_group_info = sgtbx.space_group_info(group=lattice_group)

  # Get list of sub-spacegroups
  subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()

  # Order sub-groups
  sort_values = flex.double()
  for group in subgrs:
    order_z = group.order_z()
    space_group_number = sgtbx.space_group_type(group, 00000).number()
    assert 1 <= space_group_number <= 230
    sort_values.append(order_z*1000+space_group_number)
  perm = flex.sort_permutation(sort_values, 0001)

  # Loop sub-groups in sorted order
  for i_subgr in perm:
    acentric_subgroup = subgrs[i_subgr]
    # Add centre of inversion to acentric lattice symmetry
    centric_group = sgtbx.space_group(acentric_subgroup)
    centric_group.expand_inv(sgtbx.tr_vec((0,0,0)))
    # Make symmetry object: unit-cell + space-group
    # The unit cell is potentially modified to be exactly compatible
    # with the space group symmetry.
    subsym = crystal.symmetry(
      unit_cell=minimum_symmetry.unit_cell(),
      space_group=centric_group,
      assert_is_compatible_unit_cell=00000)
    # Convert subgroup to reference setting
    cb_op_minimum_ref = subsym.space_group_info().type().cb_op()
    ref_subsym = subsym.change_basis(cb_op_minimum_ref)
    # Ignore unwanted groups
    if (not str(ref_subsym.space_group_info()) in bravais_types.centric):
      continue
    # Choose best setting for monoclinic and orthorhombic systems
    cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell()
    best_subsym = ref_subsym.change_basis(cb_op_best_cell)
    # Total basis transformation
    cb_op_inp_best = cb_op_best_cell * cb_op_minimum_ref * cb_op_inp_minimum

    subsym.space_group_info().show_summary(
      prefix="Symmetry in minimum-lengths cell: ")
    print "      Input minimum-lengths cell:", minimum_symmetry.unit_cell()
    print "           Symmetry-adapted cell:", subsym.unit_cell()
    best_subsym.space_group_info().show_summary(
      prefix="            Conventional setting: ")
    print "                       Unit cell:", best_subsym.unit_cell()
    print "                 Change of basis:", cb_op_inp_best.c()
    print "                         Inverse:", cb_op_inp_best.c_inv()
    print "      Maximal angular difference: %.3f degrees" % (
      lattice_symmetry.find_max_delta(
        minimum_cell=minimum_symmetry.unit_cell(),
        group=acentric_subgroup))
    print
  # End loop over subgroups

if (__name__ == "__main__"):
  run()

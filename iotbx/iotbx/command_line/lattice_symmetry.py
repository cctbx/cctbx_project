#! /usr/bin/env python

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
  input_symmetry = command_line.symmetry
  if (input_symmetry.unit_cell() is None):
    print
    print "***********************************"
    print "Please specify unit cell parameters"
    print "***********************************"
    print
    sys.stdout.write(command_line.parser.format_help())
    return
  if (len(command_line.args) > 0):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="Hall: %s 1" % command_line.args[0])
  elif (input_symmetry.space_group_info() is None):
    input_symmetry = crystal.symmetry(
      unit_cell=input_symmetry.unit_cell(),
      space_group_symbol="P 1")
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
  cb_op_inp_minimum = input_symmetry.change_of_basis_op_to_minimum_cell()
  minimum_symmetry = input_symmetry.change_basis(cb_op_inp_minimum)
  lattice_group = lattice_symmetry.group(
    minimum_symmetry.unit_cell(), max_delta=max_delta)
  lattice_group_info = sgtbx.space_group_info(group=lattice_group)
  subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()
  sort_values = flex.double()
  for group in subgrs:
    order_z = group.order_z()
    space_group_number = sgtbx.space_group_type(group, 00000).number()
    assert 1 <= space_group_number <= 230
    sort_values.append(order_z*1000+space_group_number)
  perm = flex.sort_permutation(sort_values, 0001)
  for i_subgrs in perm:
    subsym = crystal.symmetry(
      unit_cell=minimum_symmetry.unit_cell(),
      space_group=subgrs[i_subgrs],
      assert_is_compatible_unit_cell=00000)
    cb_op_minimum_ref = subsym.space_group_info().type().cb_op()
    ref_subsym = subsym.change_basis(cb_op_minimum_ref)
    if (not str(ref_subsym.space_group_info()) in bravais_types.acentric):
      continue
    cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell()
    best_subsym = ref_subsym.change_basis(cb_op_best_cell)
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
        group=subsym.space_group()))
    print

if (__name__ == "__main__"):
  run()

from cctbx import crystal
from cctbx import sgtbx
from cctbx.sgtbx import subgroups
from cctbx.sgtbx import lattice_symmetry
from iotbx.option_parser import iotbx_option_parser
from scitbx.array_family import flex
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
  print "Input"
  print "====="
  print
  input_symmetry.show_summary()
  print
  print "Angular tolerance:", command_line.options.delta, "degrees"
  print
  print "Similar symmetries"
  print "=================="
  print
  cb_op_inp_niggli = input_symmetry.change_of_basis_op_to_niggli_cell()
  niggli_symmetry = input_symmetry.change_basis(cb_op_inp_niggli)
  lattice_group = lattice_symmetry.group(
    niggli_symmetry.unit_cell(), max_delta=command_line.options.delta)
  lattice_group_info = sgtbx.space_group_info(group=lattice_group)
  subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()
  order_z = flex.double()
  for group in subgrs:
    order_z.append(group.order_z())
  perm = flex.sort_permutation(order_z, 0001)
  for i_subgrs in perm:
    subsym = crystal.symmetry(
      unit_cell=niggli_symmetry.unit_cell(),
      space_group=subgrs[i_subgrs],
      assert_is_compatible_unit_cell=00000)
    subsym.space_group_info().show_summary(
      prefix="Symmetry in Niggli cell: ")
    print "      Input Niggli cell:", niggli_symmetry.unit_cell()
    print "  Symmetry-adapted cell:", subsym.unit_cell()
    cb_op_niggli_ref = subsym.space_group_info().type().cb_op()
    ref_subsym = subsym.change_basis(cb_op_niggli_ref)
    cb_op_best_cell = ref_subsym.change_of_basis_op_to_best_cell()
    best_subsym = ref_subsym.change_basis(cb_op_best_cell)
    best_subsym.space_group_info().show_summary(
      prefix="   Conventional setting: ")
    print "              Unit cell:", best_subsym.unit_cell()
    cb_op_inp_best = cb_op_best_cell * cb_op_niggli_ref * cb_op_inp_niggli
    print "        Change of basis:", cb_op_inp_best.c()
    print "                Inverse:", cb_op_inp_best.c_inv()
    print

if (__name__ == "__main__"):
  run()

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
  crystal_symmetry = command_line.symmetry
  if (crystal_symmetry.unit_cell() is None):
    print "***********************************"
    print "Please specify unit cell parameters"
    print "***********************************"
    print
    sys.stdout.write(command_line.parser.format_help())
    return
  if (len(command_line.args) > 0):
    crystal_symmetry = crystal.symmetry(
      unit_cell=crystal_symmetry.unit_cell(),
      space_group_symbol="Hall: %s 1" % command_line.args[0])
  elif (crystal_symmetry.space_group_info() is None):
    crystal_symmetry = crystal.symmetry(
      unit_cell=crystal_symmetry.unit_cell(),
      space_group_symbol="P 1")
  print "Input:"
  print crystal_symmetry.unit_cell(),
  crystal_symmetry.space_group_info().show_summary()
  print "angular tolerance:", command_line.options.delta, "degrees"
  print
  print "Similar symmetries:"
  niggli_cb_op = crystal_symmetry.change_of_basis_op_to_niggli_cell()
  niggli_symmetry = crystal_symmetry.change_basis(niggli_cb_op)
  lattice_group = lattice_symmetry.group(
    niggli_symmetry.unit_cell(), max_delta=command_line.options.delta)
  lattice_group_info = sgtbx.space_group_info(group=lattice_group)
  subgrs = subgroups.subgroups(lattice_group_info).groups_parent_setting()
  order_z = flex.double()
  for group in subgrs:
    order_z.append(group.order_z())
  perm = flex.sort_permutation(order_z, 0001)
  for i_subgrs in perm:
    group = subgrs[i_subgrs]
    subsym = crystal.symmetry(
      unit_cell=niggli_symmetry.unit_cell(),
      space_group=group,
      assert_is_compatible_unit_cell=00000)
    try: subsym = subsym.change_basis(niggli_cb_op.inverse())
    except: reference_cell = niggli_symmetry.unit_cell()
    else: reference_cell = crystal_symmetry.unit_cell()
    if (subsym.unit_cell().is_similar_to(reference_cell)):
      print "*",
    else:
      print " ",
    print subsym.unit_cell(),
    subsym.space_group_info().show_summary()

if (__name__ == "__main__"):
  run()

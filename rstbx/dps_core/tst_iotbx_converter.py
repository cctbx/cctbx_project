from cctbx import crystal,sgtbx,uctbx
from iotbx.command_line import lattice_symmetry
from rstbx.dps_core.lepage import iotbx_converter
from libtbx.test_utils import show_diff

def run_one_example():
  # symmetry from 3ged
  crystal_symmetry = crystal.symmetry(
    unit_cell=uctbx.unit_cell((124.287,124.287,162.608,90.0,90.0,90.0)),
    space_group=sgtbx.space_group_info('I 41 2 2').group())

  cb_op_input_to_primitive = crystal_symmetry.change_of_basis_op_to_primitive_setting()
  primitive_symmetry = crystal_symmetry.primitive_setting()
  cb_op_input_to_minimum = primitive_symmetry.change_of_basis_op_to_minimum_cell() * cb_op_input_to_primitive

  subgroup_list = iotbx_converter(
    crystal_symmetry.change_basis(cb_op_input_to_minimum).unit_cell(),max_delta=3.0,
    bravais_types_only=False,
    space_group_symbol= str(crystal_symmetry.change_basis(cb_op_input_to_minimum).space_group_info()),
    force_minimum=False,interest_focus="input_symmetry")

  for subgroup in subgroup_list:
    subgroup['best_subsym'].show_summary()
    print

  # for comparision
  subgroup_list = lattice_symmetry.metric_subgroups(
    crystal_symmetry,3.0,bravais_types_only=False,
    best_monoclinic_beta=False).result_groups

  for subgroup in subgroup_list:
    subgroup['best_subsym'].show_summary()
    print

expected_output = """Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 41 2 2 (No. 98)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 41 (No. 80)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 21 21 21 (No. 24)

Unit cell: (162.608, 175.768, 175.768, 90, 90, 90)
Space group: F 2 2 2 (No. 22)

Unit cell: (124.287, 162.608, 124.287, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 175.768, 119.725, 90, 94.4545, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 175.768, 119.725, 90, 94.4545, 90)
Space group: I 1 2 1 (No. 5)

Unit cell: (119.725, 119.725, 119.725, 94.4545, 117.462, 117.462)
Space group: P 1 (No. 1)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 4/m m m (No. 139)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I 4/m (No. 87)

Unit cell: (124.287, 124.287, 162.608, 90, 90, 90)
Space group: I m m m (No. 71)

Unit cell: (162.608, 175.768, 175.768, 90, 90, 90)
Space group: F m m m (No. 69)

Unit cell: (175.768, 162.608, 124.287, 90, 135, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (162.608, 175.768, 119.725, 90, 132.773, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (204.667, 124.287, 124.287, 90, 127.392, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (204.667, 124.287, 124.287, 90, 127.392, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (162.608, 175.768, 119.725, 90, 132.773, 90)
Space group: C 1 2/m 1 (No. 12)

Unit cell: (119.725, 119.725, 119.725, 117.462, 117.462, 94.4545)
Space group: P -1 (No. 2)

"""
if __name__=="__main__":
  import StringIO,sys
  S = StringIO.StringIO()
  sys.stdout = S
  run_one_example()
  sys.stdout=sys.__stdout__

  assert not show_diff(S.getvalue(),expected_output)
  print "OK"

from __future__ import absolute_import, division, print_function
from iotbx import crystal_symmetry_from_any

def exercise():
  space_group = "P1"
  unit_cell = "100,100,100,90,90,90"
  symmetry = unit_cell+","+space_group
  sp_space_group = " P 1 "
  sp_symmetry = unit_cell+","+sp_space_group
  killer="100,P1"
  killer2=space_group+","+unit_cell
  killer3="100 100 100 90 90 90"
  killer4=killer3+" "+sp_space_group

  cs = crystal_symmetry_from_any.from_string(space_group)
  assert(cs is not None)
  assert(cs.unit_cell() is None)
  assert(str(cs.space_group_info())=='P 1')

  cs = crystal_symmetry_from_any.from_string(unit_cell)
  assert(cs is not None)
  assert(cs.unit_cell().parameters()==(100.0,100.0,100.0,90.0,90.0,90.0))
  assert(cs.space_group_info() is None)

  cs = crystal_symmetry_from_any.from_string(symmetry)
  assert(cs is not None)
  assert(cs.unit_cell().parameters()==(100.0,100.0,100.0,90.0,90.0,90.0))
  assert(str(cs.space_group_info())=='P 1')

  cs = crystal_symmetry_from_any.from_string(sp_space_group)
  assert(cs is not None)
  assert(cs.unit_cell() is None)
  assert(str(cs.space_group_info())=='P 1')

  cs = crystal_symmetry_from_any.from_string(sp_symmetry)
  assert(cs is not None)
  assert(cs.unit_cell().parameters()==(100.0,100.0,100.0,90.0,90.0,90.0))
  assert(str(cs.space_group_info())=='P 1')

  cs = crystal_symmetry_from_any.from_string(killer)
  assert(cs is None)

  cs = crystal_symmetry_from_any.from_string(killer2)
  assert(cs is None)

  cs = crystal_symmetry_from_any.from_string(killer3)
  assert(cs is None)

  cs = crystal_symmetry_from_any.from_string(killer4)
  assert(cs is None)

  print("OK")

if __name__=="__main__":
  exercise()

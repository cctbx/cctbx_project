from __future__ import absolute_import, division, print_function
import mmtbx.model
import iotbx.pdb
from libtbx.utils import format_cpu_times

pdb_str_1 = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
TER
"""
pdb_str_2 = """\
CRYST1   20.000   30.000   40.000  90.00  90.00  90.00 P 1 21 1
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
TER
"""
pdb_str_3= """\
CRYST1   40.000   40.000   60.000  90.00  90.00  90.00 I 4
ATOM      1   N  ILE A  40       1.000   1.000   1.000  1.00162.33           C
ATOM      2  CA  LEU A  40      94.618  -5.253  91.582  1.00 87.10           C
TER
"""


def exercise_symmetry():
  inp = iotbx.pdb.input(lines=pdb_str_1, source_info=None)
  a = mmtbx.model.manager(model_input=inp)
  inp = iotbx.pdb.input(lines=pdb_str_2, source_info=None)
  b = mmtbx.model.manager(model_input=inp)
  inp = iotbx.pdb.input(lines=pdb_str_3, source_info=None)
  c = mmtbx.model.manager(model_input=inp)
  c.set_unit_cell_crystal_symmetry(b.crystal_symmetry())
  c.set_shift_cart((23,1,7))

  print("A:",a.crystal_symmetry(),"A unit cell",a.unit_cell_crystal_symmetry())
  print("\nB:",b.crystal_symmetry(),"B unit cell:",
    b.unit_cell_crystal_symmetry())
  print("\nC:",c.crystal_symmetry(),"C unit cell:",
     c.unit_cell_crystal_symmetry())

  d = a.deep_copy()

  print('\nD:',d.crystal_symmetry(),"D unit cell:",
     d.unit_cell_crystal_symmetry())
  d.set_symmetry_and_shift_to_match_other(c)
  print('\nD matched to C:',d.crystal_symmetry(),
     "D matched to C unit cell:",d.unit_cell_crystal_symmetry())

  assert d.crystal_symmetry().is_similar_symmetry(c.crystal_symmetry())
  assert d.shift_cart() == c.shift_cart()
  assert d.unit_cell_crystal_symmetry().is_similar_symmetry(
     c.unit_cell_crystal_symmetry())

  assert not d.crystal_symmetry().is_similar_symmetry(a.crystal_symmetry())
  assert not a.unit_cell_crystal_symmetry()

def run():
  exercise_symmetry()
  print(format_cpu_times())

if (__name__ == "__main__"):
  run()

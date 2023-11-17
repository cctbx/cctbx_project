from __future__ import absolute_import, division, print_function

from libtbx.utils import format_cpu_times


def exercise_load_pdb_structure():
  from iotbx.pdb.fetch import load_pdb_structure
  h, xrs = load_pdb_structure('1yjp')
  assert h.atoms_size() == 66, h.atoms_size()
  assert xrs.scatterers().size() == 66, xrs.scatterers_size()

if (__name__ == "__main__"):
  exercise_load_pdb_structure()
  print(format_cpu_times())

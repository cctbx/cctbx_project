from __future__ import absolute_import, division, print_function
# This script generates a list of non-standard space group settings.
# The settings are used for testing.
#
# usage: python make_settings.py > settings.py

from cctbx import sgtbx
from six.moves import range

def run():
  settings = [0]
  for i in range(1, 231): settings.append({})

  list_cb_op = []
  for xyz in ("x,y,z", "z,x,y", "y,z,x"):
    list_cb_op.append(sgtbx.change_of_basis_op(sgtbx.rt_mx(xyz)))

  n_built = 0
  for i in sgtbx.space_group_symbol_iterator():
    hall_symbol = i.hall()
    for z in "PABCIRHF":
      hall_z = hall_symbol[0] + z + hall_symbol[2:]
      for cb_op in list_cb_op:
        group = sgtbx.space_group(hall_z).change_basis(cb_op)
        sg_type = group.type()
        settings[sg_type.number()][sg_type.lookup_symbol()] = 0
        n_built += 1
  print("# n_built =", n_built)

  n_non_redundant = 0
  print("settings = (")
  for i in range(1, 231):
    print("#", i)
    symbols = list(settings[i].keys())
    symbols.sort()
    for s in symbols:
      print("'" + s + "',")
      n_non_redundant += 1
  print(")")
  print("# n_non_redundant =", n_non_redundant)

if (__name__ == "__main__"):
  run()

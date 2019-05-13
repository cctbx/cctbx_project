from __future__ import division
from __future__ import print_function
from cctbx import crystal

def run():
  symmetry = crystal.symmetry(
    unit_cell=(11, 12, 13, 90, 100, 90),
    space_group_symbol="C 2")
  symmetry.show_summary()
  for s in symmetry.space_group(): print(s)

if (__name__ == "__main__"):
  run()

from __future__ import absolute_import, division, print_function
from scitbx import regular_grid_on_unit_sphere

def run():
  r = regular_grid_on_unit_sphere.rosca(m=9, hemisphere=False)
  assert r.size() == 649
  r = regular_grid_on_unit_sphere.rosca(m=9, hemisphere=True)
  assert r.size() == 361

if (__name__ == "__main__"):
  run()

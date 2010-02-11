from __future__ import division

from cctbx.array_family import flex
from cctbx import maptbx, masks, uctbx
from libtbx.test_utils import approx_equal

def exercise_flood_fill():
  uc = uctbx.unit_cell('10 10 10 90 90 90')
  gridding = maptbx.crystal_gridding(
    unit_cell=uc,
    pre_determined_n_real=(5,5,5))
  corner_cube = (0,4,20,24,100,104,120,124) # cube across all 8 corners
  channel = (12,37,38,39,42,43,62,63,67,68,87,112)
  data = flex.int(flex.grid(gridding.n_real()))
  for i in (corner_cube + channel): data[i] = 1
  flood_fill = masks.flood_fill(data)
  assert data.count(0) == 105
  for i in corner_cube: assert data[i] == 2
  for i in channel: assert data[i] == 3
  assert approx_equal(
    flood_fill.averaged_frac_coords(), ((-0.1, -0.1, -0.1), (-0.5, 7/15, 0.5)))
  assert approx_equal(
    flood_fill.averaged_indices(), ((-0.5, -0.5, -0.5), (-2.5, 7/3, 2.5)))
  assert flood_fill.n_voids() == 2
  assert approx_equal(flood_fill.grid_points_per_void(), (8, 12))
  if 0:
    from crys3d import wx_map_viewer
    wx_map_viewer.display(raw_map=data.as_double(), unit_cell=uc, wires=False)

def run():
  exercise_flood_fill()
  print "OK"

if __name__ == '__main__':
  run()

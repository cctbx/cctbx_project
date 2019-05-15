from __future__ import absolute_import, division, print_function

from cctbx.array_family import flex
from cctbx import maptbx, masks, uctbx
from libtbx.test_utils import approx_equal
from scitbx import math, matrix
from six.moves import range

def exercise_flood_fill():
  uc = uctbx.unit_cell('10 10 10 90 90 90')
  for uc in (uctbx.unit_cell('10 10 10 90 90 90'),
             uctbx.unit_cell('9 10 11 87 91 95')):
    gridding = maptbx.crystal_gridding(
      unit_cell=uc,
      pre_determined_n_real=(5,5,5))
    corner_cube = (0,4,20,24,100,104,120,124) # cube across all 8 corners
    channel = (12,37,38,39,42,43,62,63,67,68,87,112)
    data = flex.int(flex.grid(gridding.n_real()))
    for i in (corner_cube + channel): data[i] = 1
    flood_fill = masks.flood_fill(data, uc)
    assert data.count(0) == 105
    for i in corner_cube: assert data[i] == 2
    for i in channel: assert data[i] == 3
    assert approx_equal(flood_fill.centres_of_mass(),
                        ((-0.5, -0.5, -0.5), (-2.5, 7/3, 2.5)))
    assert approx_equal(flood_fill.centres_of_mass_frac(),
                        ((-0.1, -0.1, -0.1), (-0.5, 7/15, 0.5)))
    assert approx_equal(flood_fill.centres_of_mass_cart(),
                        uc.orthogonalize(flood_fill.centres_of_mass_frac()))
    assert flood_fill.n_voids() == 2
    assert approx_equal(flood_fill.grid_points_per_void(), (8, 12))
    if 0:
      from crys3d import wx_map_viewer
      wx_map_viewer.display(raw_map=data.as_double(), unit_cell=uc, wires=False)
    #
    gridding = maptbx.crystal_gridding(
      unit_cell=uc,
      pre_determined_n_real=(10,10,10))
    data = flex.int(flex.grid(gridding.n_real()))
    # parallelogram
    points = [(2,4,5),(3,4,5),(4,4,5),(5,4,5),(6,4,5),
              (3,5,5),(4,5,5),(5,5,5),(6,5,5),(7,5,5),
              (4,6,5),(5,6,5),(6,6,5),(7,6,5),(8,6,5)]
    points_frac = flex.vec3_double()
    for p in points:
      data[p] = 1
      points_frac.append([p[i]/gridding.n_real()[i] for i in range(3)])
    points_cart = uc.orthogonalize(points_frac)
    flood_fill = masks.flood_fill(data, uc)
    assert data.count(2) == 15
    assert approx_equal(flood_fill.centres_of_mass_frac(), ((0.5,0.5,0.5),))
    pai_cart = math.principal_axes_of_inertia(
      points=points_cart, weights=flex.double(points_cart.size(),1.0))
    F = matrix.sqr(uc.fractionalization_matrix())
    O = matrix.sqr(uc.orthogonalization_matrix())
    assert approx_equal(
      pai_cart.center_of_mass(), flood_fill.centres_of_mass_cart()[0])
    assert approx_equal(
      flood_fill.covariance_matrices_cart()[0],
      (F.transpose() * matrix.sym(
        sym_mat3=flood_fill.covariance_matrices_frac()[0]) * F).as_sym_mat3())
    assert approx_equal(
      pai_cart.inertia_tensor(), flood_fill.inertia_tensors_cart()[0])
    assert approx_equal(pai_cart.eigensystem().vectors(),
                        flood_fill.eigensystems_cart()[0].vectors())
    assert approx_equal(pai_cart.eigensystem().values(),
                        flood_fill.eigensystems_cart()[0].values())
  return

def run():
  exercise_flood_fill()
  print("OK")

if __name__ == '__main__':
  run()

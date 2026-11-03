from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cctbx import crystal

# map_manager.density_at_sites_cart near the edges of a map.
# For a boxed map (not full size, so wrapping=False) the density is defined
# only between the first and last sampled grid planes on each axis (grid
# units 0 to n-1). Sites beyond them, including the strip between the last
# plane and the box edge, read 0; sites on the boundary planes keep their
# values, also when fractionalization round-off puts them a little outside.
# A full-size map with wrapping=False reads 0 only outside fractional
# coordinates 0 to 1, as before boxed maps were cut off at n-1: no round-off
# tolerance, and periodic interpolation inside, including the last grid
# interval.

def linear_value(g, lower_bounds):
  '''Value of the test map at grid coordinates g of a box with lower_bounds'''
  x, y, z = [g[i] + lower_bounds[i] for i in range(3)]
  return 1 + 0.1 * x + 0.2 * y + 0.3 * z

def make_map_manager(unit_cell_parameters, wrapping):
  '''Full map, 30 x 40 x 50 grid, with density linear in the grid indices'''
  from iotbx.map_manager import map_manager
  n_real = (30, 40, 50)
  map_data = flex.double(flex.grid(n_real), 0)
  for i in range(n_real[0]):
    for j in range(n_real[1]):
      for k in range(n_real[2]):
        map_data[(i, j, k)] = linear_value((i, j, k), (0, 0, 0))
  cs = crystal.symmetry(unit_cell_parameters, "P1")
  return map_manager(map_data = map_data,
    unit_cell_grid = n_real,
    unit_cell_crystal_symmetry = cs,
    wrapping = wrapping)

def density_at_grid(mm, g):
  '''Density at grid coordinates g (floats) of map_manager mm'''
  n = mm.map_data().all()
  uc = mm.crystal_symmetry().unit_cell()
  site_frac = [g[i] / n[i] for i in range(3)]
  site_cart = uc.orthogonalize(site_frac)
  return mm.density_at_sites_cart(flex.vec3_double([site_cart]))[0]

def check_non_wrapping_box(mm, lower_bounds):
  '''Boxed map: values inside grid units 0 to n-1 on faces, edges and
  corners match the linear map; sites beyond them read 0'''
  assert not mm.wrapping()
  assert mm.map_data().origin() == (0, 0, 0)
  n = mm.map_data().all()
  tiny = 1.e-9
  for axis in range(3):
    last = n[axis] - 1
    inside = [0, tiny, -tiny, 0.5, last - 0.5, last, last - tiny, last + tiny]
    outside = [-0.5, -1.e-3, last + 1.e-3, last + 0.5, n[axis] - tiny,
      n[axis], n[axis] + 0.5]
    # Faces: one coordinate at or near a boundary, others in the middle
    for x in inside + outside:
      g = [n[i] // 2 + 0.25 for i in range(3)]
      g[axis] = x
      value = density_at_grid(mm, g)
      if x in inside:
        assert approx_equal(value, linear_value(g, lower_bounds), eps=1.e-6), (
          lower_bounds, axis, x)
      else:
        assert value == 0, (lower_bounds, axis, x, value)
  # Edges and corners: each coordinate on the first plane, on the last
  # plane, in the middle, or in the strip beyond the last plane
  choices = []
  for i in range(3):
    last = n[i] - 1
    choices.append([(0, True), (last, True), (n[i] // 2 + 0.25, True),
      (-1.e-12, True), (last + 1.e-12, True), (last + 0.5, False)])
  for cx in choices[0]:
    for cy in choices[1]:
      for cz in choices[2]:
        g = [cx[0], cy[0], cz[0]]
        value = density_at_grid(mm, g)
        if cx[1] and cy[1] and cz[1]:
          assert approx_equal(value, linear_value(g, lower_bounds), eps=1.e-6)
        else:
          assert value == 0, (lower_bounds, g, value)

def check_round_off_on_boundary_planes(mm, lower_bounds):
  '''Sites built from exact fractional coordinates 0 and (n-1)/n, which can
  come back from fractionalize a little outside the map, keep their values'''
  n = mm.map_data().all()
  uc = mm.crystal_symmetry().unit_cell()
  for axis in range(3):
    for plane in (0, n[axis] - 1):
      site_frac = [0.5, 0.5, 0.5]
      site_frac[axis] = plane / n[axis]
      g = [site_frac[i] * n[i] for i in range(3)]
      value = mm.density_at_sites_cart(
        flex.vec3_double([uc.orthogonalize(site_frac)]))[0]
      assert approx_equal(value, linear_value(g, lower_bounds), eps=1.e-6), (
        lower_bounds, axis, plane, value)

def check_full_size_non_wrapping(mm):
  '''Full-size map with wrapping=False: periodic interpolation, set to 0
  only for fractional coordinates below 0 or above 1, with no tolerance'''
  from cctbx.maptbx import real_space_target_simple_per_site
  assert not mm.wrapping()
  assert mm.is_full_size()
  n = mm.map_data().all()
  uc = mm.crystal_symmetry().unit_cell()
  sites_cart = flex.vec3_double()
  for axis in range(3):
    last = n[axis] - 1
    for x in [-0.5, -1.e-3, -1.e-9, 0, 0.5, last - 0.5, last, last + 1.e-3,
        last + 0.5, n[axis] - 1.e-9, n[axis], n[axis] + 1.e-3, n[axis] + 0.5]:
      g = [n[i] // 2 + 0.25 for i in range(3)]
      g[axis] = x
      sites_cart.append(uc.orthogonalize([g[i] / n[i] for i in range(3)]))
    site_frac = [0.5, 0.5, 0.5]
    for plane in (0, last):
      site_frac[axis] = plane / n[axis]
      sites_cart.append(uc.orthogonalize(site_frac))
  values = mm.density_at_sites_cart(sites_cart)
  expected = real_space_target_simple_per_site(
    unit_cell = uc,
    density_map = mm.map_data(),
    sites_cart = sites_cart)
  for i, site_frac in enumerate(uc.fractionalize(sites_cart)):
    if min(site_frac) < 0 or max(site_frac) > 1:
      expected[i] = 0
  assert values.all_eq(expected)
  # The last grid interval blends the last plane with plane 0
  g = [n[i] // 2 for i in range(3)]
  g[0] = n[0] - 0.5
  value = mm.density_at_sites_cart(flex.vec3_double(
    [uc.orthogonalize([g[i] / n[i] for i in range(3)])]))[0]
  g_last = [n[0] - 1, g[1], g[2]]
  g_first = [0, g[1], g[2]]
  assert approx_equal(value, 0.5 * (linear_value(g_last, (0, 0, 0)) +
    linear_value(g_first, (0, 0, 0))), eps=1.e-6)

def exercise_non_wrapping():
  '''Edge behavior of full-size and boxed maps with wrapping=False'''
  from cctbx.maptbx.box import with_bounds
  for unit_cell_parameters in [(30, 40, 50, 90, 90, 90),
       (30, 40, 50, 80, 100, 110)]:
    mm = make_map_manager(unit_cell_parameters, wrapping = False)
    check_full_size_non_wrapping(mm)
    for lower_bounds, upper_bounds in [
        ((0, 0, 0), (9, 9, 9)),
        ((5, 5, 5), (14, 14, 14)),
        ((3, 7, 11), (16, 19, 30)),
        ((20, 30, 40), (29, 39, 49))]:
      box_mm = with_bounds(mm, lower_bounds = lower_bounds,
        upper_bounds = upper_bounds).map_manager()
      check_non_wrapping_box(box_mm, lower_bounds)
      check_round_off_on_boundary_planes(box_mm, lower_bounds)

def run():
  exercise_non_wrapping()

if (__name__ == "__main__"):
  run()
  print("OK")

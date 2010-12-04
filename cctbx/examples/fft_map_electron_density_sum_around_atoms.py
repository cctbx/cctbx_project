from __future__ import division
from cctbx import miller
from cctbx import uctbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from cctbx import maptbx
from scitbx import matrix
from libtbx.math_utils import ifloor
from libtbx.test_utils import approx_equal
import random

if (1):
  random.seed(0)
  flex.set_random_seed(0)

def exercise(xray_structure, d_min, resolution_factor, rho_sum_radius):
  n_real = []
  n_half_plus = []
  n_half_minus = []
  s2 = d_min * resolution_factor * 2
  for l in xray_structure.unit_cell().parameters()[:3]:
    nh = ifloor(l / s2)
    n_real.append(2*nh+1)
    n_half_plus.append(nh)
    n_half_minus.append(-nh)
  n_real = tuple(n_real)
  n_real_product = matrix.col(n_real).product()
  crystal_gridding = maptbx.crystal_gridding(
    unit_cell=xray_structure.unit_cell(),
    space_group_info=xray_structure.space_group_info(),
    pre_determined_n_real=n_real)
  miller_indices = flex.miller_index()
  miller_indices.reserve(n_real_product)
  for h in flex.nested_loop(n_half_minus, n_half_plus, open_range=False):
    miller_indices.append(h)
  assert miller_indices.size() == n_real_product
  #
  miller_set = miller.set(
    crystal_symmetry=xray_structure,
    anomalous_flag=True,
    indices=miller_indices).sort(by_value="resolution")
  assert miller_set.indices()[0] == (0,0,0)
  f_calc = miller_set.structure_factors_from_scatterers(
    xray_structure=xray_structure,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  #
  number_of_miller_indices = []
  rho_sums_around_atoms = []
  for f in [f_calc, f_calc.resolution_filter(d_min=d_min)]:
    assert f.indices()[0] == (0,0,0)
    number_of_miller_indices.append(f.indices().size())
    fft_map = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=f)
    assert fft_map.n_real() == n_real
    rho = fft_map.real_map_unpadded() / n_real_product
    assert approx_equal(flex.sum(rho), f_calc.data()[0])
    rsaa = []
    for site_cart in xray_structure.sites_cart():
      gias = maptbx.grid_indices_around_sites(
        unit_cell=xray_structure.unit_cell(),
        fft_n_real=n_real,
        fft_m_real=n_real,
        sites_cart=flex.vec3_double([site_cart]),
        site_radii=flex.double([rho_sum_radius]))
      rsaa.append(flex.sum(rho.as_1d().select(gias)))
    rho_sums_around_atoms.append(rsaa)
  #
  def get_rho_sums(i):
    return " ".join(["%3.1f" % v for v in rho_sums_around_atoms[i]])
  print "%3.1f %4.2f %-12s %5d %5d" % (
    d_min,
    resolution_factor,
    n_real,
    number_of_miller_indices[0],
    number_of_miller_indices[1]) \
      + " | " + get_rho_sums(0) \
      + " | " + get_rho_sums(1)

def run():
  unit_cell = uctbx.unit_cell((10,10,10,90,90,90))
  rho_sum_radius = 2
  xrs0 = random_structure.xray_structure(
    unit_cell=unit_cell,
    elements=["C"]*1,
    sites_frac=[(0,0,0)])
  xrs1 = random_structure.xray_structure(
    unit_cell=unit_cell,
    elements=["C"]*3,
    min_distance=3)
  for xray_structure in [xrs0, xrs1]:
    xray_structure.show_scatterers()
    print "res fac  grid       hkl all Ewald | sum rho all | Ewald"
    for d_min in [4, 3, 2, 1]:
      for resolution_factor in [1/2, 1/3, 1/4]:
        exercise(
          xray_structure=xray_structure,
          d_min=d_min,
          resolution_factor=resolution_factor,
          rho_sum_radius=rho_sum_radius)
    print
  print "OK"

if (__name__ == "__main__"):
  run()

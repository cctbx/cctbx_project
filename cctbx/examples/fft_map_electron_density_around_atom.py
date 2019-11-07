from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import maptbx
from cctbx.adptbx import u_as_b
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.math_utils import ifloor
from libtbx.test_utils import approx_equal
from six.moves import range

def rho_stats(
      xray_structure,
      d_min,
      resolution_factor,
      electron_sum_radius,
      zero_out_f000):
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
  if (zero_out_f000):
    f_calc.data()[0] = 0j
  #
  unit_cell_volume = xray_structure.unit_cell().volume()
  voxel_volume = unit_cell_volume / n_real_product
  number_of_miller_indices = []
  rho_max = []
  electron_sums_around_atoms = []
  densities_along_x = []
  for f in [f_calc, f_calc.resolution_filter(d_min=d_min)]:
    assert f.indices()[0] == (0,0,0)
    number_of_miller_indices.append(f.indices().size())
    fft_map = miller.fft_map(
      crystal_gridding=crystal_gridding,
      fourier_coefficients=f)
    assert fft_map.n_real() == n_real
    rho = fft_map.real_map_unpadded() / unit_cell_volume
    assert approx_equal(voxel_volume*flex.sum(rho), f_calc.data()[0])
    if (xray_structure.scatterers().size() == 1):
      assert flex.max_index(rho) == 0
      rho_max.append(rho[0])
    else:
      rho_max.append(flex.max(rho))
    site_cart = xray_structure.sites_cart()[0]
    gias = maptbx.grid_indices_around_sites(
      unit_cell=xray_structure.unit_cell(),
      fft_n_real=n_real,
      fft_m_real=n_real,
      sites_cart=flex.vec3_double([site_cart]),
      site_radii=flex.double([electron_sum_radius]))
    electron_sums_around_atoms.append(
      flex.sum(rho.as_1d().select(gias))*voxel_volume)
    #
    a = xray_structure.unit_cell().parameters()[0]
    nx = n_real[0]
    nxh = nx//2
    x = []
    y = []
    for ix in range(-nxh,nxh+1):
      x.append(a*ix/nx)
      y.append(rho[(ix%nx,0,0)])
    densities_along_x.append((x,y))
  #
  print("%3.1f %4.2f %-12s %5d %5d | %6.3f %6.3f | %6.3f %6.3f | %4.2f %5.1f" % (
      d_min,
      resolution_factor,
      n_real,
      number_of_miller_indices[0],
      number_of_miller_indices[1],
      electron_sums_around_atoms[0],
      electron_sums_around_atoms[1],
      rho_max[0],
      rho_max[1],
      f_calc.data()[0].real,
      u_as_b(xray_structure.scatterers()[0].u_iso)))
  #
  return densities_along_x

table_header = """\
                          hkl     |    electrons  |     rho max
res fac  grid           all Ewald |    all  Ewald |    all  Ewald | F000  Biso\
"""

def build_xray_structure_with_carbon_along_x(a, b_iso, x=[0]):
  result = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(a,a,a,90,90,90),
      space_group_symbol="P1"))
  for i,v in enumerate(x):
    result.add_scatterer(xray.scatterer(
      label="C%d" % (i+1), scattering_type="C", site=(v,0,0), b=b_iso))
  reg = result.scattering_type_registry(table="n_gaussian", d_min=1/12)
  g = reg.as_type_gaussian_dict()["C"]
  assert g.n_terms() == 5
  assert not g.use_c()
  return result

def loop_res_fac(b, electron_sum_radius=2, zero_out_f000=False):
  xray_structure = build_xray_structure_with_carbon_along_x(a=10, b_iso=b)
  xray_structure.show_scatterers()
  print(table_header)
  for d_min in [4, 3, 2, 1]:
    for resolution_factor in [1/2, 1/3, 1/4]:
      rho_stats(
        xray_structure=xray_structure,
        d_min=d_min,
        resolution_factor=resolution_factor,
        electron_sum_radius=electron_sum_radius,
        zero_out_f000=False)

def run(args):
  assert len(args) == 0
  for b in [0, 5, 20]:
    loop_res_fac(b=b)
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

from __future__ import division
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import maptbx
from cctbx.adptbx import u_as_b
from cctbx.array_family import flex
from scitbx import matrix
from libtbx.math_utils import ifloor
from libtbx.test_utils import approx_equal

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
    assert flex.max_index(rho) == 0
    rho_max.append(rho[0])
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
    dax = []
    for ix in xrange(-nxh,nxh+1):
      dax.append((a*ix/nx, rho[(ix%nx,0,0)]))
    densities_along_x.append(dax)
  #
  print \
    "%3.1f %4.2f %-12s %5d %5d | %6.3f %6.3f | %6.3f %6.3f | %4.2f %5.1f" % (
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
      u_as_b(xray_structure.scatterers()[0].u_iso))
  #
  return densities_along_x

table_header = """\
                          hkl     |    electrons  |     rho max
res fac  grid           all Ewald |    all  Ewald |    all  Ewald | F000  Biso\
"""

def build_xray_structure_with_carbon_at_origin(a, b):
  result = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=(a,a,a,90,90,90),
      space_group_symbol="P1"))
  result.add_scatterer(xray.scatterer(
    scattering_type="C", site=(0,0,0), b=b))
  return result

def loop_res_fac(b, electron_sum_radius=2, zero_out_f000=False):
  xray_structure = build_xray_structure_with_carbon_at_origin(a=10, b=b)
  xray_structure.show_scatterers()
  print table_header
  for d_min in [4, 3, 2, 1]:
    for resolution_factor in [1/2, 1/3, 1/4]:
      densities_along_x = rho_stats(
        xray_structure=xray_structure,
        d_min=d_min,
        resolution_factor=resolution_factor,
        electron_sum_radius=electron_sum_radius,
        zero_out_f000=False)

def text_only():
  for b in [0, 5, 20]:
    loop_res_fac(b=b)

def matplotlib_plots(dry_run):
  if (not dry_run):
    import matplotlib.pyplot as plt
    fig = plt.figure()
  print table_header
  a = 5
  gr = 3
  gc = 4
  def p(i, r, f, b, z, ym):
    xs = build_xray_structure_with_carbon_at_origin(a=a, b=b)
    densities_along_x = rho_stats(
      xray_structure=xs,
      d_min=r,
      resolution_factor=f,
      electron_sum_radius=2,
      zero_out_f000=bool(z))
    if (dry_run):
      return
    ax = fig.add_subplot(gr, gc, i)
    for dax,code in reversed(zip(densities_along_x, ["b-", "r-"])):
      x_values = [x for x,y in dax]
      y_values = [y for x,y in dax]
      ax.plot(x_values, y_values, code)
    ax.axis([-a/2, a/2, -1, ym])
    if (i > 8):
      for label in ax.get_xticklabels():
        label.set_fontsize(8)
    else:
      for label in ax.get_xticklabels():
        label.set_visible(False)
    if (i % 4 == 1):
      for label in ax.get_yticklabels():
        label.set_fontsize(8)
    else:
      for label in ax.get_yticklabels():
        label.set_visible(False)
    j = len(y_values)//2
    y0 = y_values[j]
    if (y0 > ym):
      ax.text(0.5, ym*0.8, "h=%.2f" % y0, fontsize=8, color="b")
  i = 1
  for r in [1, 2, 3]:
    for f in [1/3, 1/4]:
      p(i, r, f, 0, 0, 10)
      i += 1
    p(i, r, 1/4, 20, 0, 10)
    i += 1
    p(i, r, 1/4, 20, 1, 10)
    i += 1
  if (not dry_run):
    plt.show()

def run(args):
  assert args in [[], ["plots"]]
  dry_run = False
  if (len(args) == 0):
    text_only()
    dry_run = True
  matplotlib_plots(dry_run=dry_run)
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

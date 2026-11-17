from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import sgtbx
from cctbx.development import random_structure
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
import time

random.seed(0)
flex.set_random_seed(0)

def grads_and_curvs(f_obs, xray_structure):
  """Least-squares target gradients and curvatures of xray_structure."""
  f_calc = f_obs.structure_factors_from_scatterers(
    xray_structure=xray_structure,
    algorithm="direct",
    cos_sin_table=False).f_calc()
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=False,
    obs_type="F",
    obs=f_obs.data(),
    weights=flex.double(f_obs.size(), 1),
    r_free_flags=None,
    f_calc=f_calc.data(),
    derivatives_depth=2,
    scale_factor=1)
  return xray_structure.grads_and_curvs_target_simple(
    miller_indices=f_obs.indices(),
    da_db=ls.gradients_work(),
    daa_dbb_dab=ls.hessians_work())

def exercise_site_curvs(space_group_symbol):
  """
  site_curvs holds, per scatterer, the packed upper triangle of the
  curvature block of the site parameters: its diagonal is the site part of
  curvs, and every entry matches finite differences of grads.
  """
  xs = random_structure.xray_structure(
    sgtbx.space_group_info(space_group_symbol),
    elements=("N", "C", "O", "S"),
    volume_per_atom=100,
    min_distance=2.0,
    general_positions_only=True,
    random_u_iso=True)
  f_obs = abs(xs.structure_factors(
    anomalous_flag=False, d_min=1.5, algorithm="direct").f_calc())
  # off the truth, so that the gradients and the model term are not zero
  xs.shake_sites_in_place(rms_difference=0.3)
  gact = grads_and_curvs(f_obs, xs)
  n = xs.scatterers().size()
  n_params = 7  # site (3), u_iso, occupancy, fp, fdp on a general position
  assert gact.curvs.size() == n_params * n
  assert gact.site_curvs.size() == 6 * n
  h = 1.e-5
  for i in range(n):
    block = gact.site_curvs[6 * i:6 * i + 6]
    assert approx_equal(
      [block[0], block[3], block[5]],
      gact.curvs[n_params * i:n_params * i + 3])
    scale = max(1., flex.max(flex.abs(block)))
    k = 0
    for ip in range(3):
      for jp in range(ip, 3):
        g = []
        for delta in [h, -h]:
          xp = xs.deep_copy_scatterers()
          sc = xp.scatterers()[i]
          site = list(sc.site)
          site[jp] += delta
          sc.site = tuple(site)
          g.append(grads_and_curvs(f_obs, xp).grads[n_params * i + ip])
        fd = (g[0] - g[1]) / (2 * h)
        assert approx_equal(block[k], fd, eps=1.e-5 * scale), (
          space_group_symbol, i, ip, jp, block[k], fd)
        k += 1

def run():
  for symbol in ["P 1", "P 1 21 1", "P 21 21 21"]:
    exercise_site_curvs(symbol)
  print("OK")

if __name__ == "__main__":
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time() - t0))

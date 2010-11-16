from __future__ import division
from cctbx.maptbx import real_space_target_and_gradients
from cctbx import sgtbx
from cctbx.array_family import flex
from cctbx.development import random_structure
import random, time
from cctbx.maptbx import real_space_refinement_simple
from libtbx.test_utils import approx_equal

if (1):
  random.seed(0)
  flex.set_random_seed(0)


def run():
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info  = sgtbx.space_group_info("P212121"),
    elements          = ["N","C","O","S","P"]*10,
    volume_per_atom   = 50)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.8)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  print "Start:", start_error
  map_current, miller_array, crystal_gridding = compute_map(
    xray_structure = xrs_sh)
  for step in [miller_array.d_min()/4]*5:
    if(1):
      minimized = real_space_target_and_gradients.minimization(
        xray_structure              = xrs_sh,
        miller_array                = miller_array,
        crystal_gridding            = crystal_gridding,
        map_target                  = map_target,
        max_iterations              = 500,
        min_iterations              = 25,
        step                        = step,
        geometry_restraints_manager = None,
        target_type                 = "diff_map")
      xrs_sh = minimized.xray_structure
      map_current = minimized.map_current
      final_error = flex.mean(xrs.distances(other = minimized.xray_structure))
    if(0):
      minimized = real_space_refinement_simple.lbfgs(
        sites_cart=xrs_sh.sites_cart(),
        density_map=map_target,
        unit_cell=xrs_sh.unit_cell(),
        geometry_restraints_manager=None,
        real_space_gradients_delta=step)
      xrs_sh = xrs_sh.replace_sites_cart(minimized.sites_cart)
      final_error = flex.mean(xrs.distances(other = xrs_sh))

    print "Final:", final_error
  assert approx_equal(start_error, 0.8, 1.e-3)
  assert final_error < 1.e-4
  print "OK"

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %8.3f"%(time.time()-t0)

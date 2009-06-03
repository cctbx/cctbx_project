from __future__ import division
from cctbx.maptbx import real_space_refinement_simple_2
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from cctbx.development import random_structure
import random, time
import sys

def run(args):
  assert len(args) == 0
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map
  xrs = random_structure.xray_structure(
    space_group_info = sgtbx.space_group_info("P1"),
    unit_cell        = (10, 20, 30, 70, 80, 120),
    elements         =(("O","N","C")*10),
    volume_per_atom  = 50,
    min_distance     = 2,
    u_iso            = adptbx.b_as_u(10.),
    use_u_iso        = True)
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.5)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  print "Start:", start_error
  map_current, miller_array, crystal_gridding = compute_map(
    xray_structure = xrs_sh)
  steps1 = [1.0,0.5,0.25,0.1,0.05,0.01,0.001]
  steps1.reverse()
  steps2 = [1.0,0.5,0.25,0.1,0.05,0.01,0.001]
  for step in steps1+steps2:
    minimized = real_space_refinement_simple_2.minimization(
      xray_structure   = xrs_sh,
      miller_array     = miller_array,
      crystal_gridding = crystal_gridding,
      map_target       = map_target,
      max_iterations   = 500,
      min_iterations   = 25,
      step             = step)
    xrs_sh = minimized.xray_structure
    map_current = minimized.map_current
    final_error = flex.mean(xrs.distances(other = minimized.xray_structure))
    print "Final:", final_error
  assert start_error >= 0.5
  assert final_error < 0.021
  print "OK"

if (__name__ == "__main__"):
  t0 = time.time()
  run(args=sys.argv[1:])
  print "Time: %8.3f"%(time.time()-t0)

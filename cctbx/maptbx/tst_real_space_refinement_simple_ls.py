from __future__ import division
from cctbx.maptbx import real_space_refinement_simple_ls
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.array_family import flex
from cctbx.development import random_structure
import random, time
import sys, os
import iotbx.pdb
import libtbx.load_env
import mmtbx.utils
from cctbx.maptbx import real_space_refinement_simple

def run():
  if (1):
    random.seed(0)
    flex.set_random_seed(0)
  def compute_map(xray_structure, d_min=1.5, resolution_factor=1./4):
    fc = xray_structure.structure_factors(d_min = d_min).f_calc()
    fft_map = fc.fft_map(resolution_factor=resolution_factor)
    fft_map.apply_sigma_scaling()
    result = fft_map.real_map_unpadded()
    return result, fc, fft_map  
  pdbfn = libtbx.env.find_in_repositories(
        relative_path="phenix_regression/pdb/1akg.pdb", test=os.path.isfile)
  xrs = iotbx.pdb.input(file_name = pdbfn).xray_structure_simple()
  model = mmtbx.utils.model_simple(pdb_file_names = [pdbfn],
    scattering_table="wk1995", use_elbow = False)     
  map_target,tmp,tmp = compute_map(xray_structure = xrs)
  xrs_sh = xrs.deep_copy_scatterers()
  xrs_sh.shake_sites_in_place(mean_distance=0.5)
  start_error = flex.mean(xrs.distances(other = xrs_sh))
  print "Start:", start_error
  map_current, miller_array, crystal_gridding = compute_map(
    xray_structure = xrs_sh)
  #for step in steps1+steps2:
  geometry_restraints_manager = None
  if 1:
    geometry_restraints_manager = model.restraints_manager
  for step in [miller_array.d_min()/4]*5:
    if(1):
      minimized = real_space_refinement_simple_ls.minimization(
        xray_structure              = xrs_sh,
        miller_array                = miller_array,
        crystal_gridding            = crystal_gridding,
        map_target                  = map_target,
        max_iterations              = 500,
        min_iterations              = 25,
        step                        = step,
        geometry_restraints_manager = geometry_restraints_manager)
      xrs_sh = minimized.xray_structure
      map_current = minimized.map_current
      final_error = flex.mean(xrs.distances(other = minimized.xray_structure))
    if(0):
      minimized = real_space_refinement_simple.lbfgs(
        sites_cart=xrs_sh.sites_cart(),
        density_map=map_target,
        unit_cell=xrs_sh.unit_cell(),
        geometry_restraints_manager=geometry_restraints_manager,
        real_space_gradients_delta=step)
      xrs_sh = xrs_sh.replace_sites_cart(minimized.sites_cart)
      final_error = flex.mean(xrs.distances(other = xrs_sh))
    
    print "Final:", final_error
  assert start_error >= 0.5
  assert final_error < 0.00035
  print "OK"

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print "Time: %8.3f"%(time.time()-t0)
  print "fft time: ", real_space_refinement_simple_ls.FFT_TIME

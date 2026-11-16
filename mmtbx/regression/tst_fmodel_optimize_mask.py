from __future__ import absolute_import, division, print_function
import random
import iotbx.pdb
import mmtbx.f_model
import mmtbx.masks
import mmtbx.model
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from mmtbx.secondary_structure.build import ss_idealization as ssb

def get_xray_structure():
  random.seed(0)
  flex.set_random_seed(0)
  h = ssb.secondary_structure_from_sequence(
    pdb_str=ssb.alpha_helix_str, sequence="A"*12)
  pdb_inp = iotbx.pdb.input(
    source_info=None, lines=h.as_pdb_string(crystal_symmetry=None))
  m = mmtbx.model.manager(model_input=pdb_inp, log=null_out())
  m.add_crystal_symmetry_if_necessary(box_cushion=6)
  return m.get_xray_structure()

def make_f_obs(xrs, r_solv, r_shrink, k_sol, b_sol, d_min=2.5):
  """Error-free data from the model plus a flat bulk solvent computed with
  the given mask parameters."""
  f_calc = xrs.structure_factors(d_min=d_min).f_calc()
  mp = mmtbx.masks.mask_master_params.extract()
  mp.solvent_radius = r_solv
  mp.shrink_truncation_radius = r_shrink
  f_mask = mmtbx.masks.manager(miller_array=f_calc,
    mask_params=mp).shell_f_masks(xray_structure=xrs, force_update=True)[0]
  ss = 1./flex.pow2(f_calc.d_spacings().data())/4.
  data = f_calc.data() + k_sol*flex.exp(-b_sol*ss)*f_mask.data()
  f_obs = f_calc.customized_copy(data=flex.abs(data))
  f_obs.set_observation_type_xray_amplitude()
  flags = f_obs.generate_r_free_flags(fraction=0.1, max_free=None)
  return f_obs, flags

def r_work_with_mask(f_obs, flags, xrs, r_solv, r_shrink):
  mp = mmtbx.masks.mask_master_params.extract()
  mp.solvent_radius = r_solv
  mp.shrink_truncation_radius = r_shrink
  fm = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=xrs, mask_params=mp)
  fm.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
  return fm.r_work()

def optimize_from_default(f_obs, flags, xrs, nproc):
  fm = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
    xray_structure=xrs)
  assert approx_equal(fm.mask_params.solvent_radius, 1.1)
  assert approx_equal(fm.mask_params.shrink_truncation_radius, 0.9)
  fm.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
  r_start = fm.r_work()
  fm.optimize_mask(out=null_out(), nproc=nproc)
  mp = fm.mask_params
  r_end = fm.r_work()
  assert r_end < r_start, (r_end, r_start)
  # the manager must be left with scales fitted for the chosen mask
  fm.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
  assert approx_equal(r_end, fm.r_work(), 1.e-4)
  assert mp.shrink_truncation_radius < mp.solvent_radius
  assert mp.shrink_truncation_radius > 0
  return mp.solvent_radius, mp.shrink_truncation_radius, r_end

def exercise_recovers_known_mask():
  # constant k_mask: the fast scaler represents the truth exactly, so the
  # true mask parameters are the unique minimum and must be found
  xrs = get_xray_structure()
  f_obs, flags = make_f_obs(xrs, r_solv=1.3, r_shrink=1.0, k_sol=0.35,
    b_sol=0)
  results = [optimize_from_default(f_obs, flags, xrs, nproc=n) for n in (1, 2)]
  for r_solv, r_shrink, r_end in results:
    assert approx_equal(r_solv, 1.3, 0.051), r_solv
    assert approx_equal(r_shrink, 1.0, 0.051), r_shrink
  assert approx_equal(results[0], results[1], 1.e-6)

def exercise_realistic_b_sol():
  # k_mask decaying with resolution: scales must be re-fitted for every
  # trial mask, otherwise a nonsense mask wins (observed: 0.4/0.0)
  xrs = get_xray_structure()
  f_obs, flags = make_f_obs(xrs, r_solv=0.8, r_shrink=0.5, k_sol=0.35,
    b_sol=50)
  r_true = r_work_with_mask(f_obs, flags, xrs, 0.8, 0.5)
  r_solv, r_shrink, r_end = optimize_from_default(f_obs, flags, xrs, nproc=1)
  assert r_end < r_true + 0.001, (r_end, r_true)
  assert approx_equal(r_solv, 0.8, 0.11), r_solv
  assert approx_equal(r_shrink, 0.5, 0.11), r_shrink

if (__name__ == "__main__"):
  exercise_recovers_known_mask()
  exercise_realistic_b_sol()
  print("OK")

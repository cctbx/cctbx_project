from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.f_model
import mmtbx.bulk_solvent
from libtbx import group_args
from libtbx.test_utils import approx_equal
import random

def build_fmodel(n_atoms=50, d_min=1.8, seed=0, space_group="P 21 21 21"):
  random.seed(seed)
  flex.set_random_seed(seed)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(space_group),
    elements                = (("O", "N", "C") * n_atoms),
    volume_per_atom         = 200,
    min_distance            = 1.5,
    general_positions_only  = True,
    random_u_iso            = True)
  fc = x.structure_factors(d_min=d_min, algorithm="direct").f_calc()
  f_obs = abs(fc)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags)
  fmodel.update_all_scales()
  return fmodel

def synthetic_llgi_data(fmodel, seed=1, feff_scale=1.0):
  """ feff deliberately rescaled/perturbed relative to f_obs (not merely
  a copy of it) so a test can tell whether a result actually depends on
  feff, rather than being coincidentally identical to the f_obs-based
  answer. """
  f_obs = fmodel.f_obs()
  n = f_obs.size()
  rnd = random.Random(seed)
  dobs = f_obs.array(data=flex.double([0.5 + 0.4 * rnd.random()
    for i in range(n)]))
  feff_data = f_obs.data() * feff_scale * flex.double(
    [0.85 + 0.3 * rnd.random() for i in range(n)])
  feff = f_obs.array(data=feff_data)
  teps = f_obs.array(data=flex.double(n, 1.0))
  resn = f_obs.array(data=flex.double(n, 1.0))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def exercise_llgi_target_active_gating():
  # llgi_target_active() requires BOTH target_name=="llgi" AND
  # llgi_data attached -- neither alone is sufficient.
  fmodel = build_fmodel(seed=9)
  assert not fmodel.llgi_target_active()  # neither

  llgi_data = synthetic_llgi_data(fmodel, seed=9)
  fmodel.set_llgi_data(llgi_data)
  assert not fmodel.llgi_target_active()  # llgi_data but target=ml

  fmodel2 = build_fmodel(seed=9)
  fmodel2._target_name = "llgi"  # bypass set_target_name's own Sorry
  assert fmodel2.llgi_data() is None
  assert not fmodel2.llgi_target_active()  # target=llgi but no data

  fmodel.set_target_name("llgi")
  assert fmodel.llgi_target_active()  # both

def exercise_r_work_r_free_r_all_bins_stay_f_obs_based_always():
  # r_work()/r_free()/r_all()/bins() must NEVER change meaning based on
  # target_name/llgi_data -- they stay f_obs-based always, even with
  # target=llgi and llgi_data attached AND systematically different from
  # f_obs (feff_scale=1.4). This is the key regression this test guards:
  # internal self-consistency checks throughout mmtbx.bulk_solvent.
  # f_model_all_scales (bss's own "did update_core() actually take
  # effect" assert) and elsewhere call fmodel.r_work()/r_all() directly
  # and require the f_obs-based answer regardless of target_name -- see
  # (for target=llgi, phenix.refine makes FEFF the f_obs() array itself).
  fmodel = build_fmodel(seed=20)
  llgi_data = synthetic_llgi_data(fmodel, seed=21, feff_scale=1.4)
  fmodel.set_llgi_data(llgi_data)

  r_work_before = fmodel.r_work()
  r_free_before = fmodel.r_free()
  r_all_before = fmodel.r_all()
  bins_before = fmodel.bins()

  fmodel.set_target_name("llgi")
  assert fmodel.llgi_target_active()

  assert approx_equal(fmodel.r_work(), r_work_before, eps=1.e-12)
  assert approx_equal(fmodel.r_free(), r_free_before, eps=1.e-12)
  assert approx_equal(fmodel.r_all(), r_all_before, eps=1.e-12)
  bins_after = fmodel.bins()
  assert len(bins_after) == len(bins_before)
  for b_before, b_after in zip(bins_before, bins_after):
    assert approx_equal(b_before.r, b_after.r, eps=1.e-12)
    assert approx_equal(b_before.fo_mean, b_after.fo_mean, eps=1.e-12)

  # And they must independently match a direct f_obs-based computation
  # (not just "unchanged from before" -- confirms they never touched
  # feff at all, not merely that some other bug happened to cancel out).
  r_work_direct = abs(mmtbx.bulk_solvent.r_factor(
    fmodel.f_obs_work().data(),
    fmodel.f_model_scaled_with_k1_w().data(), 1.0))
  assert approx_equal(fmodel.r_work(), r_work_direct, eps=1.e-10)

def exercise_info_in_llgi_mode():
  # In llgi mode info() reports the ordinary r_work()/r_free()/r_all()
  # and takes its likelihood statistics (FOM, phase error, D and V in
  # place of alpha and beta) from the LLGI fit. info() builds a real
  # target_functor() internally, which requires sigmaa/scatfrac to be
  # attached, so run the real per-macrocycle estimator first.
  fmodel = build_fmodel(n_atoms=70, d_min=1.7, seed=18)
  llgi_data = synthetic_llgi_data(fmodel, seed=19, feff_scale=1.3)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()

  info = fmodel.info(n_bins=5)
  assert info._llgi
  assert len(info.bins) > 0
  assert approx_equal(info.r_work, fmodel.r_work(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all(), eps=1.e-10)
  mch = fmodel.map_calculation_helper_llgi()
  assert approx_equal(info.ml_phase_error,
    flex.mean(fmodel.phase_errors_llgi(mch)), eps=1.e-10)
  assert 0 < info.fom_work_mean <= 1

def exercise_info_stays_f_obs_based_without_llgi():
  # Ordinary ml-target info() (no llgi_data at all) must be completely
  # unaffected -- same numbers as always.
  fmodel = build_fmodel(seed=23)
  info = fmodel.info(n_bins=5)
  assert approx_equal(info.r_work, fmodel.r_work(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all(), eps=1.e-10)

def exercise():
  exercise_llgi_target_active_gating()
  exercise_r_work_r_free_r_all_bins_stay_f_obs_based_always()
  exercise_info_in_llgi_mode()
  exercise_info_stays_f_obs_based_without_llgi()
  print("OK")

if (__name__ == "__main__"):
  exercise()

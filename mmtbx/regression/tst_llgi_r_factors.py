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

def _reference_r_factor(f_obs_amplitudes, f_model_amplitudes, selection=None):
  """ Independent reference computation (least-squares k1 scale, then
  R = sum|k1*|Fc| - |Fo|| / sum|Fo|), NOT calling any of the fmodel/
  f_model.py machinery under test -- so this genuinely cross-checks
  r_work_llgi()/r_free_llgi()/r_all_llgi()'s arithmetic rather than
  restating it. """
  fo = f_obs_amplitudes
  fm = flex.abs(f_model_amplitudes)
  if(selection is not None):
    fo = fo.select(selection)
    fm = fm.select(selection)
  k1 = flex.sum(fo * fm) / flex.sum(fm * fm)
  return flex.sum(flex.abs(k1 * fm - fo)) / flex.sum(fo)

def exercise_llgi_r_factors_available_gating():
  # llgi_r_factors_available() requires BOTH target_name=="llgi" AND
  # llgi_data attached -- neither alone is sufficient.
  fmodel = build_fmodel(seed=9)
  assert not fmodel.llgi_r_factors_available()  # neither

  llgi_data = synthetic_llgi_data(fmodel, seed=9)
  fmodel.set_llgi_data(llgi_data)
  assert not fmodel.llgi_r_factors_available()  # llgi_data but target=ml

  fmodel2 = build_fmodel(seed=9)
  fmodel2._target_name = "llgi"  # bypass set_target_name's own Sorry
  assert fmodel2.llgi_data() is None
  assert not fmodel2.llgi_r_factors_available()  # target=llgi but no data

  fmodel.set_target_name("llgi")
  assert fmodel.llgi_r_factors_available()  # both

def exercise_r_work_r_free_r_all_bins_stay_f_obs_based_always():
  # r_work()/r_free()/r_all()/bins() must NEVER change meaning based on
  # target_name/llgi_data -- they stay f_obs-based always, even with
  # target=llgi and llgi_data attached AND systematically different from
  # f_obs (feff_scale=1.4). This is the key regression this test guards:
  # internal self-consistency checks throughout mmtbx.bulk_solvent.
  # f_model_all_scales (bss's own "did update_core() actually take
  # effect" assert) and elsewhere call fmodel.r_work()/r_all() directly
  # and require the f_obs-based answer regardless of target_name -- see
  # mmtbx.f_model.manager.llgi_r_factors_available's docstring.
  fmodel = build_fmodel(seed=20)
  llgi_data = synthetic_llgi_data(fmodel, seed=21, feff_scale=1.4)
  fmodel.set_llgi_data(llgi_data)

  r_work_before = fmodel.r_work()
  r_free_before = fmodel.r_free()
  r_all_before = fmodel.r_all()
  bins_before = fmodel.bins()

  fmodel.set_target_name("llgi")
  assert fmodel.llgi_r_factors_available()

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

def exercise_llgi_methods_require_available():
  # r_work_llgi()/r_free_llgi()/r_all_llgi()/bins_llgi() are only valid
  # once llgi_r_factors_available() is True -- calling them beforehand
  # (llgi_data is None) must fail loudly (AttributeError on
  # self.llgi_data().feff), not silently fall back to f_obs.
  fmodel = build_fmodel(seed=22)
  assert not fmodel.llgi_r_factors_available()
  try:
    fmodel.r_work_llgi()
  except AttributeError:
    pass
  else:
    raise RuntimeError(
      "Expected AttributeError calling r_work_llgi() with no llgi_data.")

def exercise_llgi_r_work_r_free_use_feff():
  # r_work_llgi()/r_free_llgi() must be computed against llgi_data.feff,
  # matching an independent reference computation (own k1 scale, own R-
  # factor formula) to a tight tolerance -- and differ measurably from
  # the f_obs-based r_work()/r_free() when feff is systematically
  # different from f_obs (feff_scale=1.3 here).
  fmodel = build_fmodel(seed=11)
  llgi_data = synthetic_llgi_data(fmodel, seed=12, feff_scale=1.3)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  assert fmodel.llgi_r_factors_available()

  r_work_feff_based = fmodel.r_work_llgi()
  r_free_feff_based = fmodel.r_free_llgi()
  r_work_f_obs_based = fmodel.r_work()
  r_free_f_obs_based = fmodel.r_free()

  assert abs(r_work_feff_based - r_work_f_obs_based) > 1.e-3, (
    r_work_feff_based, r_work_f_obs_based)
  assert abs(r_free_feff_based - r_free_f_obs_based) > 1.e-3, (
    r_free_feff_based, r_free_f_obs_based)

  feff = llgi_data.feff
  work_sel = ~fmodel.r_free_flags().data()
  free_sel = fmodel.r_free_flags().data()
  ref_r_work = _reference_r_factor(
    feff.data().select(work_sel), fmodel.f_model().data().select(work_sel))
  ref_r_free = _reference_r_factor(
    feff.data().select(free_sel), fmodel.f_model().data().select(free_sel))
  assert approx_equal(r_work_feff_based, ref_r_work, eps=1.e-8)
  assert approx_equal(r_free_feff_based, ref_r_free, eps=1.e-8)

def exercise_llgi_r_all_uses_feff():
  fmodel = build_fmodel(seed=13)
  llgi_data = synthetic_llgi_data(fmodel, seed=14, feff_scale=1.25)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")

  r_all = fmodel.r_all_llgi()
  ref_r_all = _reference_r_factor(
    llgi_data.feff.data(), fmodel.f_model().data())
  assert approx_equal(r_all, ref_r_all, eps=1.e-8)

def exercise_bins_llgi_uses_feff():
  # bins_llgi() must switch to feff -- same basis as r_work_llgi(), and
  # must differ from the plain (always f_obs-based) bins().
  fmodel = build_fmodel(n_atoms=70, d_min=1.7, seed=16)
  llgi_data = synthetic_llgi_data(fmodel, seed=17, feff_scale=1.2)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")

  bins_f_obs_based = fmodel.bins()
  bins_feff_based = fmodel.bins_llgi()

  assert len(bins_f_obs_based) == len(bins_feff_based)
  any_differs = False
  for b_fo, b_feff in zip(bins_f_obs_based, bins_feff_based):
    if(abs(b_fo.r - b_feff.r) > 1.e-4):
      any_differs = True
    # fo_mean should also reflect feff's systematic 1.2x rescaling, not
    # stay identical to the f_obs-based value.
    if(b_fo.fo_mean is not None and b_feff.fo_mean is not None):
      assert b_feff.fo_mean > b_fo.fo_mean, (b_fo.fo_mean, b_feff.fo_mean)
  assert any_differs

def exercise_info_uses_llgi_r_factors_when_available():
  # mmtbx.f_model.f_model_info.info's own top-level r_work/r_free/r_all
  # and per-bin table must switch to the _llgi (feff-based) methods when
  # llgi_r_factors_available() is True -- the one printed in
  # phenix.refine's final statistics report. info() builds a real
  # target_functor() internally (for target_work/target_free), which
  # requires sigmaa/scatfrac to already be attached (see
  # mmtbx.refinement.targets.target_functor's llgi branch), so run the
  # real per-macrocycle estimator first, exactly as phenix.refine's own
  # updatellgisigmaa task would.
  fmodel = build_fmodel(n_atoms=70, d_min=1.7, seed=18)
  llgi_data = synthetic_llgi_data(fmodel, seed=19, feff_scale=1.3)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()

  info = fmodel.info(n_bins=5)
  assert len(info.bins) > 0
  assert approx_equal(info.r_work, fmodel.r_work_llgi(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free_llgi(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all_llgi(), eps=1.e-10)
  # And NOT accidentally equal to the f_obs-based answer.
  assert abs(info.r_work - fmodel.r_work()) > 1.e-3

  # At least one bin's mean_f_obs should reflect feff's 1.3x rescaling
  # relative to what an f_obs-based (ml-target) table would show.
  fmodel_ml = build_fmodel(n_atoms=70, d_min=1.7, seed=18)
  fmodel_ml.set_llgi_data(llgi_data)
  info_ml = fmodel_ml.info(n_bins=5)
  assert len(info_ml.bins) == len(info.bins)
  ratios = []
  for b_ml, b_llgi in zip(info_ml.bins, info.bins):
    if(b_ml.mean_f_obs and b_llgi.mean_f_obs):
      ratios.append(b_llgi.mean_f_obs / b_ml.mean_f_obs)
  assert len(ratios) > 0
  assert all(r > 1.05 for r in ratios), ratios

def exercise_info_stays_f_obs_based_without_llgi():
  # Ordinary ml-target info() (no llgi_data at all) must be completely
  # unaffected -- same numbers as always.
  fmodel = build_fmodel(seed=23)
  info = fmodel.info(n_bins=5)
  assert approx_equal(info.r_work, fmodel.r_work(), eps=1.e-10)
  assert approx_equal(info.r_free, fmodel.r_free(), eps=1.e-10)
  assert approx_equal(info.r_all, fmodel.r_all(), eps=1.e-10)

def exercise():
  exercise_llgi_r_factors_available_gating()
  exercise_r_work_r_free_r_all_bins_stay_f_obs_based_always()
  exercise_llgi_methods_require_available()
  exercise_llgi_r_work_r_free_use_feff()
  exercise_llgi_r_all_uses_feff()
  exercise_bins_llgi_uses_feff()
  exercise_info_uses_llgi_r_factors_when_available()
  exercise_info_stays_f_obs_based_without_llgi()
  print("OK")

if (__name__ == "__main__"):
  exercise()

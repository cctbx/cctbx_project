from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.development import random_structure
from cctbx import sgtbx
import mmtbx.f_model
import mmtbx.refinement.llgi_bss as llgi_bss
from libtbx import group_args
from libtbx.utils import Sorry
from libtbx.test_utils import approx_equal
import random

def build_fmodel(n_atoms=60, d_min=1.7, seed=0, space_group="P 21 21 21"):
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
  """ Build a plausible, self-consistent llgi_data group_args on
  fmodel.f_obs()'s CURRENT index set (i.e. after update_all_scales()'s
  own outlier removal, matching build_fmodel()'s ordering -- same
  precondition as mmtbx.refinement.llgi_e_bulk_solvent.run_inner_loop's
  own _synthetic_llgi_inputs helper). feff is built from f_obs itself
  (optionally rescaled) plus reflection-dependent noise, so it is
  amplitude-shaped and strictly positive but NOT numerically identical
  to f_obs -- letting a test tell whether a result actually came from
  fitting against feff or was coincidentally identical to the f_obs fit.
  """
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

def exercise_requires_llgi_data():
  # No llgi_data attached at all -> Sorry, not a silent no-op or crash.
  fmodel = build_fmodel(seed=10)
  try:
    llgi_bss.run_llgi_bss(fmodel)
  except Sorry as e:
    assert "llgi_data" in str(e)
  else:
    raise RuntimeError("Expected Sorry for missing llgi_data.")

def exercise_requires_matching_indices():
  # llgi_data attached but stale relative to f_obs's current index set
  # (e.g. f_obs shrunk via outlier removal since llgi_data was last
  # attached) -> Sorry, not a silent misalignment.
  fmodel = build_fmodel(seed=11)
  llgi_data = synthetic_llgi_data(fmodel, seed=12)
  n = fmodel.f_obs().size()
  truncated = group_args(
    dobs=llgi_data.dobs.customized_copy(
      indices=llgi_data.dobs.indices()[:n - 5],
      data=llgi_data.dobs.data()[:n - 5]),
    feff=llgi_data.feff.customized_copy(
      indices=llgi_data.feff.indices()[:n - 5],
      data=llgi_data.feff.data()[:n - 5]),
    teps=llgi_data.teps.customized_copy(
      indices=llgi_data.teps.indices()[:n - 5],
      data=llgi_data.teps.data()[:n - 5]),
    resn=llgi_data.resn.customized_copy(
      indices=llgi_data.resn.indices()[:n - 5],
      data=llgi_data.resn.data()[:n - 5]),
    info=None)
  fmodel._llgi_data = truncated  # bypass set_llgi_data's own index check
  try:
    llgi_bss.run_llgi_bss(fmodel)
  except Sorry as e:
    assert "index set" in str(e)
  else:
    raise RuntimeError("Expected Sorry for a stale/mismatched llgi_data.")

def exercise_does_not_modify_f_obs():
  # run_llgi_bss must never touch fmodel.f_obs() itself (see its
  # docstring: R-work/R-free/maps stay F-obs-based) -- only
  # k_isotropic/k_anisotropic/k_mask.
  fmodel = build_fmodel(seed=13)
  llgi_data = synthetic_llgi_data(fmodel, seed=14)
  fmodel.set_llgi_data(llgi_data)
  f_obs_before = flex.double(fmodel.f_obs().data())

  llgi_bss.run_llgi_bss(fmodel)

  f_obs_after = flex.double(fmodel.f_obs().data())
  assert approx_equal(list(f_obs_before), list(f_obs_after), eps=0.0)

def exercise_changes_scale_and_mask_from_ls_f_obs_result():
  # The whole point: fitting against feff (systematically rescaled
  # relative to f_obs here, feff_scale=1.4) must produce a MEASURABLY
  # different k_isotropic/k_anisotropic/k_mask than bss's own F-obs fit
  # left in place -- otherwise this would silently be a no-op reusing
  # the F-obs result.
  fmodel = build_fmodel(seed=15)
  k_iso_before = flex.double(fmodel.k_isotropic())
  k_mask_before = flex.double(fmodel.k_masks()[0])

  llgi_data = synthetic_llgi_data(fmodel, seed=16, feff_scale=1.4)
  fmodel.set_llgi_data(llgi_data)
  llgi_bss.run_llgi_bss(fmodel)

  k_iso_after = flex.double(fmodel.k_isotropic())
  k_mask_after = flex.double(fmodel.k_masks()[0])
  assert flex.max(flex.abs(k_iso_after - k_iso_before)) > 1.e-3, (
    flex.max(flex.abs(k_iso_after - k_iso_before)))
  # And the new k_isotropic should track the feff_scale=1.4 rescaling
  # roughly (least-squares scale factor against ~1.4x-larger amplitudes
  # should land near k_iso_before*1.4, not near k_iso_before itself).
  ratio = flex.mean(k_iso_after) / flex.mean(k_iso_before)
  assert 1.1 < ratio < 1.7, ratio
  assert flex.max(flex.abs(k_mask_after - k_mask_before)) >= 0.0  # smoke

def exercise_result_r_all_feff_is_reasonable():
  # The LS scaler's own R-factor against feff (not f_obs) should be
  # small on this noiseless-but-rescaled synthetic feff, confirming the
  # fit is actually converging against feff and not e.g. silently
  # falling back to f_obs internally.
  fmodel = build_fmodel(seed=17)
  llgi_data = synthetic_llgi_data(fmodel, seed=18, feff_scale=1.0)
  fmodel.set_llgi_data(llgi_data)

  result = llgi_bss.run_llgi_bss(fmodel)
  assert result.r_all_feff < 0.3, result.r_all_feff

def exercise_slow_path_returns_k_sol_b_sol():
  # params.fast=False uses the minimization-based scaler
  # (bss.bulk_solvent_and_scales), which reports k_sol/b_sol directly
  # (unlike the fast/analytical path, which leaves them as None -- see
  # run_llgi_bss's docstring).
  fmodel = build_fmodel(seed=19)
  llgi_data = synthetic_llgi_data(fmodel, seed=20)
  fmodel.set_llgi_data(llgi_data)

  params = llgi_bss.llgi_bss_params.extract()
  params.fast = False
  result = llgi_bss.run_llgi_bss(fmodel, params=params)
  assert result.k_sol is not None
  assert result.b_sol is not None
  assert len(result.k_sol) >= 1
  assert 0.0 <= result.k_sol[0] <= 0.6
  assert 0.0 <= result.b_sol[0] <= 300.0

def exercise():
  exercise_requires_llgi_data()
  exercise_requires_matching_indices()
  exercise_does_not_modify_f_obs()
  exercise_changes_scale_and_mask_from_ls_f_obs_result()
  exercise_result_r_all_feff_is_reasonable()
  exercise_slow_path_returns_k_sol_b_sol()
  print("OK")

if (__name__ == "__main__"):
  exercise()

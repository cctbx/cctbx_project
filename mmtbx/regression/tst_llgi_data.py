from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx import sgtbx
from cctbx.array_family import flex
import mmtbx.f_model
from libtbx.utils import Sorry
import random

def build_fmodel():
  random.seed(0)
  flex.set_random_seed(0)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 4"),
    elements                = (("O","N","C")*5),
    volume_per_atom         = 200,
    min_distance            = 1.5,
    general_positions_only  = True,
    random_u_iso            = True)
  fc = x.structure_factors(d_min=2.5, algorithm="direct").f_calc()
  f_obs = abs(fc)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags)
  return fmodel

def make_llgi_arrays(f_obs, resize=None):
  """ Build synthetic dobs/feff/teps/resn arrays on (a possibly resized
  subset of) f_obs's index set, matching the shapes get_llgi_data would
  produce, without going through an actual reflection file. """
  from libtbx import group_args
  indices = f_obs.indices()
  if(resize is not None):
    indices = indices[:resize]
    f_obs = f_obs.customized_copy(indices=indices,
      data=f_obs.data()[:resize])
  n = indices.size()
  dobs = f_obs.array(data=flex.double(n, 0.7))
  feff = f_obs.array(data=f_obs.data())
  teps = f_obs.array(data=flex.double(n, 1.0))
  resn = f_obs.array(data=flex.double(n, 1.0))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def exercise_set_llgi_data_ok():
  fmodel = build_fmodel()
  assert fmodel.llgi_data() is None
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  fmodel.set_llgi_data(llgi_data)
  assert fmodel.llgi_data() is not None
  assert fmodel.llgi_data().dobs.indices().all_eq(fmodel.f_obs().indices())

def exercise_set_llgi_data_mismatched_indices():
  # An array on a different (smaller) index set must be rejected rather
  # than silently accepted -- set_llgi_data() requires callers to have
  # already completed/common_set every array against f_obs (see
  # phenix.refinement.llgi_data.get_llgi_data, which does this).
  fmodel = build_fmodel()
  n_full = fmodel.f_obs().indices().size()
  assert n_full > 5
  llgi_data = make_llgi_arrays(fmodel.f_obs(), resize=n_full - 2)
  try:
    fmodel.set_llgi_data(llgi_data)
  except Sorry:
    pass
  else:
    assert False, "expected Sorry to be raised"

def exercise_set_llgi_data_missing_component():
  from libtbx import group_args
  fmodel = build_fmodel()
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  incomplete = group_args(
    dobs=llgi_data.dobs, feff=llgi_data.feff, teps=None,
    resn=llgi_data.resn, info=None)
  try:
    fmodel.set_llgi_data(incomplete)
  except Sorry:
    pass
  else:
    assert False, "expected Sorry to be raised"

def exercise_set_target_name_llgi_requires_data():
  fmodel = build_fmodel()
  assert fmodel.llgi_data() is None
  try:
    fmodel.set_target_name("llgi")
  except Sorry:
    pass
  else:
    assert False, "expected Sorry to be raised"
  # once attached, set_target_name should succeed
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  assert fmodel.target_name == "llgi"
  # target_attributes() must resolve the newly-registered "llgi" family
  attr = fmodel.target_attributes()
  assert attr.family == "llgi"
  assert attr.specialization is None

def exercise_target_functor_reports_missing_sigmaa_scatfrac():
  # The sigmaA/ScatFrac estimator is not yet implemented (see
  # doc/llgi_target_design.md sec. 5); target_functor() must fail with a
  # clear Sorry, not an opaque AttributeError/TypeError, when llgi_data
  # lacks sigmaa/scatfrac (which it always will until that estimator
  # exists and attaches them).
  fmodel = build_fmodel()
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  try:
    fmodel.target_functor()
  except Sorry:
    pass
  else:
    assert False, "expected Sorry to be raised"

def exercise_update_llgi_sigmaa_scatfrac_enables_target_functor():
  # Once mmtbx.refinement.llgi_sigmaa's estimator exists (it does now --
  # see tst_llgi_sigmaa.py for its own correctness tests), calling
  # fmodel.update_llgi_sigmaa_scatfrac() should attach sigmaa/scatfrac to
  # llgi_data and let target_functor() succeed, completing the gap
  # exercise_target_functor_reports_missing_sigmaa_scatfrac documents.
  fmodel = build_fmodel()
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  assert getattr(fmodel.llgi_data(), "sigmaa", None) is None
  fmodel.update_llgi_sigmaa_scatfrac()
  assert getattr(fmodel.llgi_data(), "sigmaa", None) is not None
  assert getattr(fmodel.llgi_data(), "scatfrac", None) is not None
  # dobs/feff/teps/resn must survive the update unchanged.
  assert fmodel.llgi_data().dobs.indices().all_eq(fmodel.f_obs().indices())
  result = fmodel.target_functor()
  core_result = result(compute_gradients=True)
  # A finite target_work is enough here (correctness of the underlying
  # estimator and target math is covered by tst_llgi_sigmaa.py and
  # cctbx/xray/targets/tst_llgi.py respectively) -- this test is only
  # about the plumbing between update_llgi_sigmaa_scatfrac() and
  # target_functor() actually working end to end.
  import math
  assert math.isfinite(core_result.core_result.target_work())

def exercise_update_llgi_sigmaa_scatfrac_requires_llgi_data():
  fmodel = build_fmodel()
  assert fmodel.llgi_data() is None
  try:
    fmodel.update_llgi_sigmaa_scatfrac()
  except Sorry:
    pass
  else:
    assert False, "expected Sorry to be raised"

def exercise_llgi_data_survives_select_and_update_all_scales():
  # Real bug, found running target=llgi against a real (non-synthetic)
  # dataset for the first time: manager.select() (used internally by
  # manager.remove_outliers(), itself called from every ordinary
  # f_model_all_scales.compute() / bulk-solvent-and-scaling pass -- not
  # a rare code path) constructs a *fresh* manager(...) object directly,
  # separately from the self.__init__(...) re-init path
  # set_target_name/_validate_and_set_llgi_data's docstrings already
  # discuss. Before llgi_data was threaded through select()'s manager(...)
  # call the same way abcd already was, this silently dropped llgi_data
  # on the very first bulk-solvent-and-scaling step of a real refinement,
  # producing a "no LLGI data provided" Sorry despite it having been
  # attached and even used successfully one step earlier -- a much more
  # commonly hit path than the low-resolution-outlier self.__init__(...)
  # branch this test's sibling coverage (implicitly, via
  # exercise_update_llgi_sigmaa_scatfrac_enables_target_functor) does not
  # exercise at all.
  fmodel = build_fmodel()
  llgi_data = make_llgi_arrays(fmodel.f_obs())
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  assert fmodel.llgi_data() is not None
  # select() with an all-True selection (a deep_copy, in effect) must
  # still carry llgi_data (all components, including sigmaa/scatfrac)
  # through to the new manager.
  selection = flex.bool(fmodel.f_obs().data().size(), True)
  new_fmodel = fmodel.select(selection=selection, deep_copy_xray_structure=False)
  assert new_fmodel.llgi_data() is not None
  assert getattr(new_fmodel.llgi_data(), "sigmaa", None) is not None
  assert getattr(new_fmodel.llgi_data(), "scatfrac", None) is not None
  assert new_fmodel.llgi_data().dobs.indices().all_eq(
    new_fmodel.f_obs().indices())
  # target_functor() must still succeed on the selected copy without
  # needing update_llgi_sigmaa_scatfrac() to be called again.
  new_fmodel.target_functor()

def build_larger_fmodel(seed=20):
  # A bigger, higher-resolution structure than build_fmodel() above --
  # needed for update_llgi_e_bulk_solvent tests, since the E-scale bulk-
  # solvent inner loop's Chebyshev/kernel-based SigmaP fit and its
  # sigmaA spline fit both want more than the ~15-atom/2.5A structure
  # build_fmodel() uses for the simpler set_llgi_data-only tests.
  random.seed(seed)
  flex.set_random_seed(seed)
  x = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 21 21 21"),
    elements                = (("O", "N", "C") * 60),
    volume_per_atom         = 200,
    min_distance            = 1.5,
    general_positions_only  = True,
    random_u_iso            = True)
  fc = x.structure_factors(d_min=1.7, algorithm="direct").f_calc()
  f_obs = abs(fc)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  fmodel = mmtbx.f_model.manager(
    xray_structure = x,
    f_obs          = f_obs,
    r_free_flags   = r_free_flags)
  fmodel.update_all_scales()
  return fmodel

def make_e_scale_llgi_arrays(f_obs, seed=21):
  # Synthetic-but-plausible nacelle-like dobs/feff/resn (TEPS == 1
  # throughout -- tNCS is not supported, see phenix.refinement.
  # llgi_data.check_teps_no_tncs), sized against f_obs's CURRENT index
  # set (i.e. AFTER any outlier removal update_all_scales() may already
  # have performed -- see mmtbx.refinement.llgi_e_bulk_solvent.
  # run_inner_loop's docstring for why this ordering matters).
  from libtbx import group_args
  n = f_obs.size()
  epsilons = f_obs.epsilons().data().as_double()
  rnd = random.Random(seed)
  dobs = f_obs.array(
    data=flex.double([0.5 + 0.4 * rnd.random() for i in range(n)]))
  feff = f_obs.array(data=f_obs.data() * flex.double(
    [0.9 + 0.2 * rnd.random() for i in range(n)]))
  resn = f_obs.array(data=flex.sqrt(epsilons) * flex.double(
    [rnd.uniform(2.0, 6.0) for i in range(n)]))
  teps = f_obs.array(data=flex.double(n, 1.0))
  return group_args(dobs=dobs, feff=feff, teps=teps, resn=resn, info=None)

def exercise_update_llgi_e_bulk_solvent_updates_k_mask():
  fmodel = build_larger_fmodel(seed=22)
  llgi_data = make_e_scale_llgi_arrays(fmodel.f_obs(), seed=23)
  fmodel.set_llgi_data(llgi_data)

  import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bulk_solvent
  params = llgi_e_bulk_solvent.llgi_e_bulk_solvent_params.extract()
  # Explicitly the two-stage (LLGI-driven bulk solvent) path: this test
  # asserts k_mask CHANGES, which only happens when Stage 2 runs.
  # fix_bulk_solvent_from_ls now defaults to True, which deliberately
  # leaves k_mask untouched (see mmtbx.regression.
  # tst_llgi_e_bulk_solvent.exercise_fixed_bulk_solvent_leaves_k_mask_
  # untouched for the default path's own coverage).
  params.fix_bulk_solvent_from_ls = False
  params.max_inner_iterations = 3
  params.sigmaa_max_iterations = 20
  params.bulk_solvent_max_iterations = 15

  k_mask_before = flex.double(fmodel.k_masks()[0])
  result = fmodel.update_llgi_e_bulk_solvent(params=params)
  k_mask_after = flex.double(fmodel.k_masks()[0])

  assert result.n_iterations >= 1
  assert flex.max(flex.abs(k_mask_before - k_mask_after)) > 0
  # update_llgi_e_bulk_solvent must NOT touch llgi_data.sigmaa/scatfrac
  # (the F-scale curve consumed directly by the llgi coordinate-
  # refinement target, mmtbx/refinement/targets.py) -- the E-scale
  # sigmaA fit it computes internally is a genuinely different curve
  # and must stay separate. See that method's docstring.
  assert not hasattr(fmodel.llgi_data(), "sigmaa") or \
    fmodel.llgi_data().sigmaa is None
  assert not hasattr(fmodel.llgi_data(), "scatfrac") or \
    fmodel.llgi_data().scatfrac is None

def exercise_update_llgi_e_bulk_solvent_requires_llgi_data():
  fmodel = build_larger_fmodel(seed=24)
  try:
    fmodel.update_llgi_e_bulk_solvent()
  except Sorry as e:
    assert "LLGI data" in str(e)
  else:
    raise RuntimeError("Expected Sorry for missing LLGI data.")

def exercise_e_scale_bulk_solvent_phil_scope_parses():
  # The refinement.llgi_data.e_scale_bulk_solvent phil scope (phenix/
  # phenix/refinement/__init__.params) must parse and its "include scope
  # mmtbx.refinement.llgi_e_bulk_solvent.llgi_e_bulk_solvent_params"
  # directive must resolve to the same fields run_inner_loop expects
  # (n_sigmap_nodes, k_sol_min/max, etc.) alongside its own .enabled
  # switch -- exercised directly against mmtbx.refinement.
  # llgi_e_bulk_solvent's own phil scope (not the full phenix
  # __init__.params file, which is outside this module's reach) to catch
  # any field-name drift between the two without a phenix-side test.
  import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bulk_solvent
  extract = llgi_e_bulk_solvent.llgi_e_bulk_solvent_params.extract()
  for name in ("n_sigmap_nodes", "auto_kernel_number", "n_sigmaa_coeffs",
               "spline_degree", "sigmaa_max_iterations", "k_sol_min",
               "k_sol_max", "b_sol_min", "b_sol_max",
               "bulk_solvent_max_iterations", "max_inner_iterations",
               "convergence_tolerance", "sigmaa_curvature_weight",
               "fix_bulk_solvent_from_ls"):
    assert hasattr(extract, name), name

def exercise():
  exercise_set_llgi_data_ok()
  exercise_set_llgi_data_mismatched_indices()
  exercise_set_llgi_data_missing_component()
  exercise_set_target_name_llgi_requires_data()
  exercise_target_functor_reports_missing_sigmaa_scatfrac()
  exercise_update_llgi_sigmaa_scatfrac_enables_target_functor()
  exercise_update_llgi_sigmaa_scatfrac_requires_llgi_data()
  exercise_llgi_data_survives_select_and_update_all_scales()
  exercise_update_llgi_e_bulk_solvent_updates_k_mask()
  exercise_update_llgi_e_bulk_solvent_requires_llgi_data()
  exercise_e_scale_bulk_solvent_phil_scope_parses()
  print("OK")

if (__name__ == "__main__"):
  exercise()

from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.utils import Sorry
from mmtbx.regression.llgi_test_utils import build_fmodel

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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  # target_functor() must fail with a clear Sorry, not an opaque
  # AttributeError/TypeError, when llgi_data lacks sigmaa/scatfrac (i.e.
  # before update_llgi_sigmaa_scatfrac() has run).
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  # fmodel.update_llgi_sigmaa_scatfrac() attaches sigmaa/scatfrac to
  # llgi_data, after which target_functor() succeeds (the counterpart of
  # exercise_target_functor_reports_missing_sigmaa_scatfrac).
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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
  fmodel = build_fmodel(5, 2.5, space_group="P 4", update_scales=False)
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

def exercise_e_scale_sigmaa_phil_scope_parses():
  # The refinement.llgi_data.e_scale_sigmaa phil scope (phenix/phenix/
  # refinement/__init__.params) includes mmtbx.refinement.
  # llgi_e_bulk_solvent.llgi_e_bulk_solvent_params; check the fields the
  # E-scale sigmaA fit reads are all present, including the d_model
  # sub-scope, to catch field-name drift without a phenix-side test.
  import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bulk_solvent
  extract = llgi_e_bulk_solvent.llgi_e_bulk_solvent_params.extract()
  for name in ("sigmaa_model", "n_sigmap_nodes", "auto_kernel_number",
               "n_sigmaa_coeffs", "spline_degree", "sigmaa_max_iterations",
               "sigmaa_curvature_weight"):
    assert hasattr(extract, name), name
  for name in ("n_gaussian_terms", "max_iterations", "include_constant_term"):
    assert hasattr(extract.d_model_params, name), name

def exercise():
  exercise_set_llgi_data_ok()
  exercise_set_llgi_data_mismatched_indices()
  exercise_set_llgi_data_missing_component()
  exercise_set_target_name_llgi_requires_data()
  exercise_target_functor_reports_missing_sigmaa_scatfrac()
  exercise_update_llgi_sigmaa_scatfrac_enables_target_functor()
  exercise_update_llgi_sigmaa_scatfrac_requires_llgi_data()
  exercise_llgi_data_survives_select_and_update_all_scales()
  exercise_e_scale_sigmaa_phil_scope_parses()
  print("OK")

if (__name__ == "__main__"):
  exercise()

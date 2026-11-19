import random
import time
from io import StringIO
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.development import random_structure
from libtbx.test_utils import approx_equal
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.f_model
import mmtbx.masks
import mmtbx.bulk_solvent.bulk_solvent_and_scaling as bss

K_SOL  = 0.35
B_SOL  = 50.0
B_CART = [3., -5., 2., 0, 0, 0]

def get_params(edits=()):
  user = iotbx.phil.parse("\n".join(edits))
  working, unused = bss.master_params.fetch(
    source=user, track_unused_definitions=True)
  assert len(unused) == 0
  return working.extract()

def get_fmodel(symbol="P212121", twin_law=None, n_radial_shells=1,
               b_cart=B_CART):
  """
  fmodel with all scales at their initial values and f_obs calculated with
  known k_sol, b_sol, b_cart.
  """
  random.seed(0)
  flex.set_random_seed(0)
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               = (("O","N","C")*50),
    volume_per_atom        = 100,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)
  f_dummy = abs(xrs.structure_factors(d_min=3.0).f_calc())
  flags = f_dummy.generate_r_free_flags(fraction=0.1, max_free=99999999)
  mask_params = mmtbx.masks.mask_master_params.extract()
  mask_params.n_radial_shells = n_radial_shells
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_dummy,
    r_free_flags   = flags,
    xray_structure = xrs,
    mask_params    = mask_params)
  kbu = fmodel.fmodel_kbu()
  kbu.update(k_sols=[K_SOL]*n_radial_shells, b_sols=[B_SOL]*n_radial_shells,
    b_cart=b_cart)
  f_obs = abs(kbu.f_model)
  return mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = flags,
    xray_structure = xrs,
    mask_params    = mask_params,
    twin_law       = twin_law)

def exercise_mode_from_params():
  # mode=slow from params: exact answer, reported in the log
  fmodel = get_fmodel()
  log = StringIO()
  r = fmodel.update_all_scales(params=get_params(["mode=slow"]), log=log,
    remove_outliers=False)
  assert approx_equal(r.k_sol[0], K_SOL, 1.e-3)
  assert approx_equal(r.b_sol[0], B_SOL, 0.1)
  assert approx_equal(r.b_cart, B_CART, 0.01)
  assert fmodel.r_work() < 1.e-4
  text = log.getvalue()
  for line in [
    "Bulk-solvent and scaling: mode=slow",
    "  bulk_solvent=True anisotropic_scaling=True apply_back_trace=False",
    "  grid search: k_sol 0.00-0.60 step 0.10, b_sol 10.0-80.0 step 5.0",
    "  minimization: k_sol,b_sol=True b_cart=True (symmetry constraints=True)",
    "  bounds: k_sol 0.00-0.60, b_sol 0.0-300.0",
    "  fixed: none",
    "  k_sol:       0.350",
    "  b_sol:       50.0"]:
    assert text.count(line) == 1, (line, text)
  # default params: fast
  fmodel = get_fmodel()
  log = StringIO()
  r = fmodel.update_all_scales(params=get_params(), log=log,
    remove_outliers=False)
  text = log.getvalue()
  assert text.count("Bulk-solvent and scaling: mode=fast") == 1
  assert text.count("grid search") == 0
  assert text.count("  k_sol:") == 0
  assert r.k_sol is None
  r_fast = fmodel.r_work()
  # no params, no log: fast, same result
  fmodel = get_fmodel()
  fmodel.update_all_scales(remove_outliers=False)
  assert approx_equal(fmodel.r_work(), r_fast, 1.e-9)
  # legacy: explicit fast=False with default params is the slow protocol
  fmodel = get_fmodel()
  r = fmodel.update_all_scales(params=get_params(), fast=False,
    remove_outliers=False)
  assert approx_equal(r.k_sol[0], K_SOL, 1.e-3)
  # explicit fast=True against mode=slow: programmer error
  fmodel = get_fmodel()
  error = None
  try:
    fmodel.update_all_scales(params=get_params(["mode=slow"]), fast=True)
  except AssertionError as e:
    error = e
  assert error is not None

def exercise_log_fixed_and_off():
  fmodel = get_fmodel()
  log = StringIO()
  p = get_params(["mode=slow", "slow.k_sol_b_sol_grid_search=False",
    "slow.minimization_k_sol_b_sol=False", "slow.fix_k_sol=0.2",
    "slow.fix_b_sol=30"])
  r = fmodel.update_all_scales(params=p, log=log, remove_outliers=False)
  text = log.getvalue()
  for line in ["  grid search: off",
    "  minimization: k_sol,b_sol=False b_cart=True (symmetry constraints=True)",
    "  fixed: k_sol=0.2000 b_sol=30.0000",
    "  k_sol:       0.2000",
    "  b_sol:       30.0000"]:
    assert text.count(line) == 1, (line, text)

def exercise_sorry_through_fmodel():
  for edits, n_shells, text in [
      (["slow.k_sol_max=0.5"], 1, "mode=fast"),
      ([], 2, "mode=slow"),
      (["mode=slow", "slow.fix_k_sol=0.3"], 1, "slow.fix_k_sol")]:
    fmodel = get_fmodel(n_radial_shells=n_shells)
    sorry = None
    try:
      fmodel.update_all_scales(params=get_params(edits),
        remove_outliers=False)
    except Sorry as e:
      sorry = str(e)
    assert sorry is not None and sorry.count(text) > 0, (edits, sorry)

def exercise_several_masks_back_trace_skipped():
  fmodel = get_fmodel(n_radial_shells=2)
  log = StringIO()
  b_before = fmodel.xray_structure.extract_u_iso_or_u_equiv().deep_copy()
  r = fmodel.update_all_scales(params=get_params(["mode=slow"]), log=log,
    remove_outliers=False, apply_back_trace=True)
  text = log.getvalue()
  assert text.count("  apply_back_trace skipped: more than one mask") == 1
  assert len(r.k_sol) == 2
  assert approx_equal(b_before,
    fmodel.xray_structure.extract_u_iso_or_u_equiv())

def exercise_apply_back_trace_param():
  result = []
  for flag in [True, False]:
    fmodel = get_fmodel()
    u_before = flex.mean(fmodel.xray_structure.extract_u_iso_or_u_equiv())
    p = get_params(["mode=slow", "apply_back_trace=%s"%str(flag)])
    log = StringIO()
    r = fmodel.update_all_scales(params=p, remove_outliers=False,
      apply_back_trace=True, log=log)
    assert log.getvalue().count("apply_back_trace=%s"%str(flag)) == 1
    u_after = flex.mean(fmodel.xray_structure.extract_u_iso_or_u_equiv())
    result.append([adptbx.u_as_b(u_after-u_before), r.b_cart,
      fmodel.r_work()])
  trace = sum(B_CART[:3])/3.
  # True: atoms take the trace, b_cart becomes traceless
  assert approx_equal(result[0][0], trace, 0.01)
  assert approx_equal(sum(result[0][1][:3]), 0, 0.01)
  # False: params wins over the argument
  assert approx_equal(result[1][0], 0, 1.e-6)
  assert approx_equal(result[1][1], B_CART, 0.01)
  assert result[0][2] < 1.e-3 and result[1][2] < 1.e-4

def exercise_twin():
  # default parameters must run with twinned data
  fmodel = get_fmodel(symbol="P31", twin_law="-h,-k,l",
    b_cart=[3., 3., -6., 0, 0, 0])
  assert fmodel.twinned()
  log = StringIO()
  fmodel.update_all_scales(params=get_params(), log=log,
    remove_outliers=False)
  text = log.getvalue()
  assert text.count(
    "Bulk-solvent and scaling: twin protocol (mode and slow.* not used)") > 0
  # twin_f_model passes fast=False with default params: must run too
  fmodel.update_all_scales(params=get_params(), fast=False,
    remove_outliers=False)
  sorry = None
  try:
    fmodel.update_all_scales(params=get_params(["mode=slow"]),
      remove_outliers=False)
  except Sorry as e:
    sorry = str(e)
  assert sorry is not None and sorry.count("twin") > 0

def run():
  exercise_mode_from_params()
  exercise_log_fixed_and_off()
  exercise_sorry_through_fmodel()
  exercise_several_masks_back_trace_skipped()
  exercise_apply_back_trace_param()
  exercise_twin()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")

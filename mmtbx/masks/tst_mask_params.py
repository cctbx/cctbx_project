import random
import time
from io import StringIO
from cctbx.array_family import flex
from cctbx import sgtbx
from cctbx.development import random_structure
from libtbx.utils import Sorry
import iotbx.phil
import mmtbx.f_model
import mmtbx.masks

def get_params(edits=()):
  user = iotbx.phil.parse("\n".join(edits))
  working, unused = mmtbx.masks.mask_master_params.fetch(
    source=user, track_unused_definitions=True)
  assert len(unused) == 0, [str(u) for u in unused]
  return working.extract()

def get_xrs(symbol="P212121", elements=("O","N","C")*30, seed=0,
            volume_per_atom=100):
  random.seed(seed)
  flex.set_random_seed(seed)
  return random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info(symbol=symbol),
    elements               = elements,
    volume_per_atom        = volume_per_atom,
    min_distance           = 1.5,
    general_positions_only = True,
    random_u_iso           = True,
    random_occupancy       = False)

def f_mask(xrs, edits, d_min=2.0):
  ma = xrs.structure_factors(d_min=d_min).f_calc()
  mm = mmtbx.masks.manager(
    miller_array   = ma,
    xray_structure = xrs,
    mask_params    = get_params(edits))
  return mm.shell_f_masks()[0]

def mean_diff(a, b):
  return flex.mean(flex.abs(a.data()-b.data()))

def exercise_scope():
  p = mmtbx.masks.mask_master_params.extract()
  assert not hasattr(p, "verbose")
  for o in mmtbx.masks.mask_master_params.objects:
    assert o.help is not None and len(o.help) > 0, o.name
    if(o.name == "ignore_zero_occupancy_atoms"):
      assert o.help.lower().startswith("do not include"), o.help

def exercise_validate():
  cases = [
    (["n_real=40,48,54"], "mask.n_real"),
    (["step=None"], "mask.step"),
    (["use_resolution_based_gridding=True", "grid_step_factor=None"],
      "mask.grid_step_factor"),
    (["grid_step_factor=0"], "mask.grid_step_factor"),
    (["n_radial_shells=2", "use_asu_masks=False"], "mask.n_radial_shells"),
    (["solvent_radius=-1"], "mask.solvent_radius"),
    (["shrink_truncation_radius=-1"], "mask.shrink_truncation_radius"),
    (["step=0"], "mask.step"),
    (["step=None", "n_real=40,0,54"], "mask.n_real"),
    (["radial_shell_width=0", "n_radial_shells=2"],
      "mask.radial_shell_width"),
    (["n_radial_shells=0"], "mask.n_radial_shells"),
    (["mean_shift_for_mask_update=-1"], "mask.mean_shift_for_mask_update"),
    (["Fmask_res_high=-1"], "mask.Fmask_res_high"),
  ]
  for edits, text in cases:
    sorry = None
    try:
      mmtbx.masks.validate_mask_params(get_params(edits))
    except Sorry as e:
      sorry = str(e)
    assert sorry is not None, edits
    assert sorry.count(text) > 0, (edits, sorry)
  # legal
  mmtbx.masks.validate_mask_params(get_params())
  mmtbx.masks.validate_mask_params(get_params(["step=None",
    "n_real=40,48,54"]))
  mmtbx.masks.validate_mask_params(get_params(["step=None",
    "use_resolution_based_gridding=True"]))
  mmtbx.masks.validate_mask_params(get_params(["n_radial_shells=3"]))
  mmtbx.masks.validate_mask_params(get_params(["solvent_radius=0",
    "shrink_truncation_radius=0", "Fmask_res_high=0",
    "mean_shift_for_mask_update=0"]))
  # the mask manager checks
  xrs = get_xrs()
  sorry = None
  try:
    f_mask(xrs, ["n_radial_shells=2", "use_asu_masks=False"])
  except Sorry as e:
    sorry = str(e)
  assert sorry is not None and sorry.count("mask.n_radial_shells") > 0

def exercise_gridding():
  """
  n_real and resolution-based gridding are used with and without asu masks.
  """
  xrs = get_xrs()
  ref = f_mask(xrs, [])
  assert flex.mean(flex.abs(ref.data())) > 1
  uc = xrs.unit_cell().parameters()[:3]
  # coarse grids: all different from the default grid
  n_real = [int(x/1.5)+1 for x in uc]
  n_real = ",".join([str(n+n%2) for n in n_real])
  result = []
  for asu in ["True", "False"]:
    u = ["use_asu_masks=%s"%asu]
    f_step  = f_mask(xrs, u)
    f_nreal = f_mask(xrs, u+["step=None", "n_real=%s"%n_real])
    f_res3  = f_mask(xrs, u+["use_resolution_based_gridding=True",
      "grid_step_factor=3"])
    f_res5  = f_mask(xrs, u+["use_resolution_based_gridding=True",
      "grid_step_factor=5"])
    # asu or whole unit cell: same result (P212121)
    assert mean_diff(f_step, ref) < 1.e-6, (asu, mean_diff(f_step, ref))
    for f in [f_nreal, f_res3]:
      assert mean_diff(f, f_step) > 1, (asu, mean_diff(f, f_step))
    assert mean_diff(f_res3, f_res5) > 1
    result.append([f_nreal, f_res3, f_res5])
  for a, b in zip(result[0], result[1]):
    assert mean_diff(a, b) < 1.e-6, mean_diff(a, b)

def exercise_asu_and_whole_cell_agree():
  """
  Same F_mask with use_asu_masks True and False, whatever defines the grid,
  in space groups where the grid depends on more than the space group
  symmetry alone.
  """
  for symbol, d_min in [("P21", 1.5), ("F222", 2.5), ("P41", 3.3),
                        ("C2", 3.3)]:
    xrs = get_xrs(symbol=symbol, elements=("O","N","C")*10, seed=1,
      volume_per_atom=150)
    p = get_params(["use_resolution_based_gridding=True"])
    n_real_asu = mmtbx.masks.asu_mask(xray_structure=xrs, d_min=d_min,
      mask_params=p).mask_data_whole_uc().focus()
    assert tuple(n_real_asu) == tuple(mmtbx.masks.gridding_n_real(
      mask_params=p, xray_structure=xrs, d_min=d_min)), symbol
    for edits in [[], ["use_resolution_based_gridding=True"],
                  ["use_resolution_based_gridding=True",
                   "grid_step_factor=3"]]:
      edits = edits + ["Fmask_res_high=0"]
      a = f_mask(xrs, edits+["use_asu_masks=True"],  d_min=d_min)
      b = f_mask(xrs, edits+["use_asu_masks=False"], d_min=d_min)
      assert flex.mean(flex.abs(a.data())) > 0.1
      assert mean_diff(a, b) < 1.e-6, (symbol, edits, mean_diff(a, b))

def exercise_fmask_res_high():
  """
  Fmask_res_high is the only resolution cutoff of F_mask.
  """
  xrs = get_xrs()
  for asu in ["True", "False"]:
    u = ["use_asu_masks=%s"%asu]
    for cutoff in [3.0, 2.5, 0.0]:
      f = f_mask(xrs, u+["Fmask_res_high=%s"%cutoff], d_min=2.0)
      d = f.d_spacings().data()
      a = flex.abs(f.data())
      if(cutoff > 0):
        assert flex.max(a.select(d < cutoff)) == 0
      inside = a.select(d >= cutoff)
      assert (inside > 0).count(True) > 0.9*inside.size(), (asu, cutoff)
      if(cutoff == 2.5):
        between = a.select((d >= 2.5) & (d < 3.0))
        assert between.size() > 10
        assert (between > 0).count(True) > 0.9*between.size(), asu

def exercise_number_of_mask_calculations():
  xrs = get_xrs()
  ma = xrs.structure_factors(d_min=3.0).f_calc()
  for asu in ["True", "False"]:
    result = []
    for threshold in [0.5, 0.0]:
      mm = mmtbx.masks.manager(
        miller_array   = ma,
        xray_structure = xrs,
        mask_params    = get_params(["use_asu_masks=%s"%asu,
          "mean_shift_for_mask_update=%s"%threshold]))
      n0 = mmtbx.masks.number_of_mask_calculations
      xrs_shaken = xrs.deep_copy_scatterers()
      xrs_shaken.shake_sites_in_place(mean_distance=0.05)
      mm.shell_f_masks(xray_structure = xrs_shaken)
      result.append(mmtbx.masks.number_of_mask_calculations - n0)
    assert result == [0, 1], (asu, result)

def exercise_log():
  """
  The scaling log tells which mask parameters are used.
  """
  def get_log(edits, hydrogens=False, neutron=False):
    elements = ("O","N","C")*30
    if(hydrogens): elements = elements + ("H",)*30
    xrs = get_xrs(elements=elements)
    if(neutron): xrs.switch_to_neutron_scattering_dictionary()
    f_obs = abs(xrs.structure_factors(d_min=3.0).f_calc())
    flags = f_obs.generate_r_free_flags(fraction=0.1, max_free=99999999)
    data_type = "xray"
    if(neutron): data_type = "neutron"
    fmodel = mmtbx.f_model.manager(
      f_obs          = f_obs,
      r_free_flags   = flags,
      xray_structure = xrs,
      data_type      = data_type,
      mask_params    = get_params(edits))
    log = StringIO()
    fmodel.update_all_scales(log=log, remove_outliers=False)
    return log.getvalue()
  text = get_log([])
  for line in [
    "  mask: solvent_radius=0.95 shrink_truncation_radius=0.70 "
      "use_asu_masks=True",
    "  mask: gridding from step=0.60",
    "  mask: ignore_zero_occupancy_atoms=True ignore_hydrogens=True",
    "  mask: Fmask_res_high=3.00 n_radial_shells=1"]:
    assert text.count(line) == 1, (line, text)
  text = get_log(["step=None", "n_real=20,24,30", "use_asu_masks=False",
    "solvent_radius=0.8", "Fmask_res_high=2", "ignore_hydrogens=False",
    "ignore_zero_occupancy_atoms=False"])
  for line in [
    "  mask: solvent_radius=0.80 shrink_truncation_radius=0.70 "
      "use_asu_masks=False",
    "  mask: gridding from n_real=20,24,30",
    "  mask: ignore_zero_occupancy_atoms=False ignore_hydrogens=False",
    "  mask: Fmask_res_high=2.00 n_radial_shells=1"]:
    assert text.count(line) == 1, (line, text)
  text = get_log(["use_resolution_based_gridding=True", "grid_step_factor=5"])
  assert text.count(
    "  mask: gridding from d_min/grid_step_factor=3.00/5.00") == 1, text
  # neutron data: H are always in the mask, and the log says so
  text = get_log([], hydrogens=True, neutron=True)
  assert text.count("  mask: ignore_zero_occupancy_atoms=True "
    "ignore_hydrogens=False (always False for neutron data)") == 1, text

def exercise_neutron_does_not_change_callers_params():
  xrs = get_xrs(elements=("O","N","C")*30+("H",)*30)
  xrs.switch_to_neutron_scattering_dictionary()
  f_obs = abs(xrs.structure_factors(d_min=3.0).f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.1, max_free=99999999)
  p = get_params()
  fmodel = mmtbx.f_model.manager(
    f_obs          = f_obs,
    r_free_flags   = flags,
    xray_structure = xrs,
    data_type      = "neutron",
    mask_params    = p)
  assert fmodel.data_type() == "neutron"
  assert fmodel.mask_params.ignore_hydrogens is False
  assert p.ignore_hydrogens is True

def run():
  exercise_scope()
  exercise_validate()
  exercise_gridding()
  exercise_asu_and_whole_cell_agree()
  exercise_fmask_res_high()
  exercise_number_of_mask_calculations()
  exercise_log()
  exercise_neutron_does_not_change_callers_params()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")

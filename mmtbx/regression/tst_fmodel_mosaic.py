from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from mmtbx.bulk_solvent import mosaic
from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
from mmtbx.bulk_solvent import tst_mosaic as tm
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import mmtbx.f_model
import random, time

"""
Mosaic bulk-solvent model inside mmtbx.f_model.manager:
update_all_scales(mosaic=True), reuse of the stored correction, maps.
"""

def scale(fm, **kwds):
  """update_all_scales with the settings used throughout this test."""
  args = dict(params = slow_params(), remove_outliers = False,
    apply_scale_k1_to_f_obs = False)
  args.update(kwds)
  return fm.update_all_scales(**args)

def slow_params():
  params = bss.master_params.extract()
  params.mode = "slow"
  return params

def new_fmodel(o):
  return mmtbx.f_model.manager(
    f_obs          = o.f_obs,
    r_free_flags   = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers(),
    mask_params    = o.mask_params)

def f_model_expected(fm, correction):
  """k_total*(F_calc + k_mask*(F_mask + correction)), no H in these models."""
  k_total = fm.k_isotropic()*fm.k_anisotropic()
  return k_total*(fm.f_calc().data() +
    fm.k_masks()[0]*(fm.f_masks()[0].data() + correction))

def is_zero(a):
  return flex.abs(a).all_eq(0)

def exercise_off(o):
  """
  Without mosaic=True nothing changes, whatever the parameters say.
  """
  fm = new_fmodel(o)
  assert fm.mosaic_correction is None and fm.mosaic_result() is None
  n0 = mosaic.number_of_runs
  scale(fm)
  r_off = fm.r_work()
  assert is_zero(fm.f_part1().data())
  assert fm.mosaic_correction is None
  assert approx_equal(fm.f_model().data(),
    f_model_expected(fm, flex.complex_double(fm.f_obs().data().size(), 0)))
  # explicit False, non-default mosaic settings
  fm2 = new_fmodel(o)
  params = slow_params()
  params.mosaic.volume_min = 10
  params.mosaic.algorithm = "alg2"
  scale(fm2, params = params, mosaic = False)
  assert approx_equal(fm2.r_work(), r_off, 1.e-9)
  assert approx_equal(fm2.f_model().data(), fm.f_model().data(), 1.e-9)
  assert is_zero(fm2.f_part1().data())
  assert mosaic.number_of_runs == n0
  # remove_mosaic with nothing stored
  fm2.remove_mosaic()
  scale(fm2)
  assert approx_equal(fm2.r_work(), r_off, 1.e-9)
  return r_off

def exercise_on(o, r_off):
  fm = new_fmodel(o)
  n0 = mosaic.number_of_runs
  sio = StringIO()
  scale(fm, mosaic = True, log = sio)
  assert mosaic.number_of_runs == n0+1
  assert "mosaic: regions 3, small 0, fitted 3, <q> " in sio.getvalue()
  assert fm.r_work() < r_off - 0.005, [fm.r_work(), r_off]
  res = fm.mosaic_result()
  assert res.n_fitted == 3
  assert approx_equal(tm.occupancies(res), o.q, 0.08), tm.occupancies(res)
  c = fm.mosaic_correction
  assert c.data().size() == fm.f_obs().data().size()
  assert c.indices().all_eq(fm.f_obs().indices())
  assert approx_equal(c.data(), res.correction)
  assert fm.mosaic_correction_is_current()
  # F_model is what it is documented to be
  assert approx_equal(fm.f_model().data(), f_model_expected(fm, c.data()),
    1.e-6)
  assert approx_equal(fm.f_part1().data(), fm.k_masks()[0]*c.data(), 1.e-6)
  assert not is_zero(fm.f_part1().data())
  # f_model_no_scales is F_calc plus the flat bulk solvent, as it always was
  assert approx_equal((fm.f_model_no_scales().data() + fm.f_part1().data())*
    fm.k_isotropic()*fm.k_anisotropic(), fm.f_model().data(), 1.e-6)
  # F_mask itself is not changed
  ref = new_fmodel(o)
  assert approx_equal(fm.f_masks()[0].data(), ref.f_masks()[0].data())
  # settings come from params.mosaic
  fm2 = new_fmodel(o)
  params = slow_params()
  params.mosaic.volume_min = 1.e6
  scale(fm2, params = params, mosaic = True)
  assert fm2.mosaic_result().n_small == 3
  assert fm2.mosaic_result().n_fitted == 0
  params.mosaic.volume_min = -1
  try: scale(fm2, params = params, mosaic = True)
  except Exception as e: assert "volume_min" in str(e)
  else: raise RuntimeError("Exception expected")
  # no parameters at all: defaults (fast scaling)
  fm3 = new_fmodel(o)
  fm3.update_all_scales(remove_outliers = False,
    apply_scale_k1_to_f_obs = False)
  r3 = fm3.r_work()
  fm3.update_all_scales(mosaic = True, remove_outliers = False,
    apply_scale_k1_to_f_obs = False)
  assert fm3.mosaic_result().n_fitted == 3
  assert fm3.r_work() < r3
  return fm

def exercise_reuse(o, fm_on):
  """
  Calls without mosaic=True keep the stored correction and do no region
  analysis; remove_mosaic() goes back to the flat model.
  """
  r_on = fm_on.r_work()
  c = fm_on.mosaic_correction.data().deep_copy()
  n0 = mosaic.number_of_runs
  for fm in [fm_on.deep_copy(), fm_on]:
    assert fm.mosaic_correction is not None
    assert fm.mosaic_result() is fm_on.mosaic_result()
    for i in range(2):
      scale(fm)
      assert approx_equal(fm.r_work(), r_on, 1.e-4)
      assert approx_equal(fm.mosaic_correction.data(), c)
      assert approx_equal(fm.f_model().data(), f_model_expected(fm, c), 1.e-6)
  # default parameters, fast scaling
  fm = fm_on.deep_copy()
  fm.update_all_scales(remove_outliers = False,
    apply_scale_k1_to_f_obs = False)
  assert approx_equal(fm.mosaic_correction.data(), c)
  assert approx_equal(fm.f_model().data(), f_model_expected(fm, c), 1.e-6)
  assert mosaic.number_of_runs == n0
  # The correction belongs to the mask it was computed for. Atoms moved
  # and the mask is the same (not updated): kept.
  fm = fm_on.deep_copy()
  xrs = fm.xray_structure.deep_copy_scatterers()
  xrs.shake_sites_in_place(mean_distance = 0.02)
  fm.update_xray_structure(xray_structure = xrs, update_f_calc = True,
    update_f_mask = False)
  assert fm.mosaic_correction_is_current()
  scale(fm)
  assert approx_equal(fm.mosaic_correction.data(), c)
  assert approx_equal(fm.f_model().data(), f_model_expected(fm, c), 1.e-6)
  assert mosaic.number_of_runs == n0
  # Atoms moved (same number) and the mask computed again: the regions are
  # not those of the correction any more, it is dropped.
  fm = fm_on.deep_copy()
  xrs = fm.xray_structure.deep_copy_scatterers()
  xrs.shake_sites_in_place(mean_distance = 0.5)
  f_mask_before = fm.f_masks()[0].data().deep_copy()
  fm.update_xray_structure(xray_structure = xrs, update_f_calc = True,
    update_f_mask = True)
  assert flex.max(flex.abs(fm.f_masks()[0].data() - f_mask_before)) > 1
  assert not fm.mosaic_correction_is_current()
  sio = StringIO()
  scale(fm, log = sio)
  assert "mosaic: stored correction dropped" in sio.getvalue()
  assert fm.mosaic_correction is None and is_zero(fm.f_part1().data())
  assert mosaic.number_of_runs == n0
  # a selection of reflections of a current correction is current
  fm = fm_on.deep_copy()
  sel = flex.random_bool(fm.f_obs().data().size(), 0.8)
  assert fm.select(selection = sel).mosaic_correction_is_current()
  # back to the flat model
  fm = fm_on.deep_copy()
  fm.remove_mosaic()
  assert fm.mosaic_correction is None and fm.mosaic_result() is None
  scale(fm)
  assert is_zero(fm.f_part1().data())
  ref = new_fmodel(o)
  scale(ref)
  assert approx_equal(fm.r_work(), ref.r_work(), 1.e-9)
  assert approx_equal(fm.f_model().data(), ref.f_model().data(), 1.e-9)
  assert mosaic.number_of_runs == n0
  # fm_on itself is still mosaic
  assert approx_equal(fm_on.r_work(), r_on, 1.e-4)

def exercise_select(o, fm_on):
  """
  Selections of reflections carry the correction.
  """
  c = fm_on.mosaic_correction
  sel = flex.random_bool(c.data().size(), 0.7)
  fm = fm_on.select(selection = sel)
  assert fm.mosaic_correction.data().size() == sel.count(True)
  assert approx_equal(fm.mosaic_correction.data(), c.data().select(sel))
  assert fm.mosaic_correction.indices().all_eq(fm.f_obs().indices())
  assert approx_equal(fm.f_model().data(),
    fm_on.f_model().data().select(sel))
  scale(fm)
  assert approx_equal(fm.f_model().data(),
    f_model_expected(fm, fm.mosaic_correction.data()), 1.e-6)
  fm = fm_on.resolution_filter(d_min = 4.5)
  n = fm.f_obs().data().size()
  assert n < c.data().size()
  assert fm.mosaic_correction.data().size() == n
  scale(fm)
  assert 0 < fm.r_work() < 0.1
  # in place, as outlier rejection does
  fm = fm_on.deep_copy()
  fm.select(selection = sel, in_place = True)
  assert fm.mosaic_correction.data().size() == sel.count(True)
  # outlier rejection during scaling
  fm = new_fmodel(o)
  data = o.f_obs.data().deep_copy()
  i_max = flex.max_index(data)
  data[i_max] = data[i_max]*50 # an outlier
  fm = mmtbx.f_model.manager(
    f_obs          = o.f_obs.array(data = data),
    r_free_flags   = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers(),
    mask_params    = o.mask_params)
  scale(fm, mosaic = True, remove_outliers = True)
  n = fm.f_obs().data().size()
  assert n < data.size()
  assert fm.mosaic_correction.data().size() == n
  assert fm.f_part1().data().size() == n
  assert approx_equal(fm.f_model().data(),
    f_model_expected(fm, fm.mosaic_correction.data()), 1.e-6)

def exercise_atoms_change(o, fm_on):
  """
  Number of atoms changed since the correction was computed: the regions it
  belongs to are gone, so it is dropped.
  """
  n0 = mosaic.number_of_runs
  fm = fm_on.deep_copy()
  xrs = fm.xray_structure
  sel = flex.bool(xrs.scatterers().size(), True)
  for i in range(5): sel[i] = False
  fm.update_xray_structure(xray_structure = xrs.select(sel),
    update_f_calc = True, update_f_mask = True)
  sio = StringIO()
  scale(fm, log = sio)
  assert fm.mosaic_correction is None
  assert fm.mosaic_result() is None
  assert is_zero(fm.f_part1().data())
  assert "mosaic: stored correction dropped" in sio.getvalue()
  assert mosaic.number_of_runs == n0
  # and computed again on request
  scale(fm, mosaic = True)
  assert fm.mosaic_correction is not None
  assert fm.mosaic_correction_is_current()
  assert mosaic.number_of_runs == n0+1

def difference_density(fm, sites_frac):
  m = fm.electron_density_map().fft_map(
    resolution_factor = 0.25, map_type = "mFobs-DFmodel")
  m.apply_sigma_scaling()
  md = m.real_map_unpadded()
  return flex.mean(flex.double(
    [md.eight_point_interpolation(s) for s in sites_frac]))

def exercise_maps():
  """
  What the mosaic model does to the difference map.
  """
  # Empty cavity: the flat model fills it with solvent, which shows as
  # negative density; with the mosaic model it is gone.
  o = tm.synthetic(q = (0., 0.5, 0.8))
  nx,ny,nz = o.regions.n_real
  sites = flex.vec3_double()
  for p in list(o.regions.secondary[0].points)[::7]:
    sites.append(((p//(ny*nz))/nx, ((p//nz)%ny)/ny, (p%nz)/nz))
  assert sites.size() > 20
  fm_off = new_fmodel(o)
  scale(fm_off)
  d_off = difference_density(fm_off, sites)
  assert d_off < -1.5, d_off
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  d_on = difference_density(fm, sites)
  assert abs(d_on) < 0.3, d_on
  assert not hasattr(fm, "for_maps")
  # Atoms that are in the crystal but not in the model, in a cavity without
  # solvent: the region gets occupancy one (the limit), which is what the
  # flat model has there. The atoms show in the difference map as they do
  # with the flat model.
  center, radius = tm.SHELLS[2]
  extra = tm.sites_in_sphere(center = center, radius = 3.6, spacing = 1.6)
  o = tm.synthetic(q = (1., 1., 0.), extra = extra)
  uc = o.xray_structure.unit_cell()
  sites = flex.vec3_double([uc.fractionalize(s) for s in extra])
  fm_off = new_fmodel(o)
  scale(fm_off)
  d_off = difference_density(fm_off, sites)
  assert d_off > 2, d_off
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  q = tm.occupancies(fm.mosaic_result())
  assert q[2] == 1 and max(q) <= 1, q
  d_on = difference_density(fm, sites)
  assert d_on > d_off - 0.3, [d_on, d_off]
  # Without the limit the solvent model takes the atoms for dense solvent
  # and removes them from the map (why the limit is there).
  fm_raw = fm.deep_copy()
  sel = mosaic.work_selection(fm_off)
  k = (fm_off.k_isotropic()*fm_off.k_anisotropic()*fm_off.k_masks()[0])
  G = [(f*k).select(sel) for f in o.F]
  whole = (o.f_whole*k).select(sel)
  q_raw = mosaic.fit_occupancies(
    f_obs = fm_off.f_obs().data().select(sel),
    atoms = fm_off.f_model().data().select(sel) - whole,
    rest  = whole - G[0] - G[1] - G[2], G = G, q_max = None)
  assert q_raw[2] > 1.2, list(q_raw)
  c = flex.complex_double(o.f_whole.size(), 0)
  for qn, f in zip(q_raw, o.F): c += f*(qn-1.)
  fm_raw._set_mosaic_arrays(correction = fm.f_obs().array(data = c),
    signature = fm.mosaic_signature(), result = fm.mosaic_result())
  scale(fm_raw)
  d_raw = difference_density(fm_raw, sites)
  assert d_raw < d_off - 0.5, [d_raw, d_off]

def exercise_skip(o):
  """
  Twinned data and neutron data: skipped with a message, scales as without
  mosaic.
  """
  def new(twin_law = None, neutron = False):
    xrs = o.xray_structure.deep_copy_scatterers()
    if(neutron): xrs.switch_to_neutron_scattering_dictionary()
    return mmtbx.f_model.manager(
      f_obs          = o.f_obs,
      r_free_flags   = o.r_free_flags,
      xray_structure = xrs,
      twin_law       = twin_law,
      mask_params    = o.mask_params)
  for kwds, text in [(dict(neutron = True), "mosaic: skipped (neutron data)"),
                     (dict(twin_law = "k,h,-l"),
                      "mosaic: skipped (twinned data)")]:
    n0 = mosaic.number_of_runs
    fm_off, fm_on = new(**kwds), new(**kwds)
    args = dict(remove_outliers = False, apply_scale_k1_to_f_obs = False)
    fm_off.update_all_scales(**args)
    sio = StringIO()
    fm_on.update_all_scales(mosaic = True, log = sio, **args)
    assert text in sio.getvalue(), sio.getvalue()
    assert fm_on.mosaic_correction is None and fm_on.mosaic_result() is None
    assert approx_equal(fm_on.r_work(), fm_off.r_work(), 1.e-6)
    assert approx_equal(fm_on.f_model().data(), fm_off.f_model().data(), 1.e-6)
    assert mosaic.number_of_runs == n0

def with_riding_h(o):
  """
  The model of o with one H at zero occupancy (riding H) next to each atom.
  """
  from cctbx import xray, adptbx
  xrs = o.xray_structure.deep_copy_scatterers()
  sites = xrs.sites_cart()
  n = xrs.scatterers().size()
  uc = xrs.unit_cell()
  for i in range(n):
    s = sites[i]
    xrs.add_scatterer(xray.scatterer(label = "H%d"%i,
      site = uc.fractionalize((s[0]+0.6, s[1]+0.6, s[2]+0.5)),
      u = adptbx.b_as_u(20), scattering_type = "H", occupancy = 0))
  assert xrs.hd_selection().count(True) == n
  return xrs

def check_parts(fm):
  """
  F_model = k_total*(F_calc + F_H + k_mask*(F_mask + correction)), the
  mosaic part being f_part1 and the H part f_part2.
  """
  c = fm.mosaic_correction.data()
  k_total = fm.k_isotropic()*fm.k_anisotropic()
  expected = k_total*(fm.f_calc().data() + fm.f_part2().data() +
    fm.k_masks()[0]*(fm.f_masks()[0].data() + c))
  assert approx_equal(fm.f_model().data(), expected, 1.e-6)
  assert approx_equal(fm.f_part1().data(), fm.k_masks()[0]*c, 1.e-6)

def fast_params():
  params = bss.master_params.extract()
  params.mode = "fast"
  return params

def exercise_hydrogens():
  """
  Model with riding H at zero occupancy (H contribution goes to f_part2):
  both contributions are kept apart, F_model is consistent. Both scaling
  modes, and for several free sets: the H search may or may not find a
  contribution.
  """
  for seed in range(4):
    random.seed(seed)
    flex.set_random_seed(seed)
    o = tm.synthetic(q = (0., 0.5, 0.8))
    xrs = with_riding_h(o)
    for params in [slow_params(), fast_params()]:
      fm = mmtbx.f_model.manager(f_obs = o.f_obs,
        r_free_flags = o.r_free_flags, xray_structure = xrs,
        mask_params = o.mask_params)
      scale(fm, params = params)
      r_off = fm.r_work()
      scale(fm, params = params, mosaic = True)
      assert fm.mosaic_result().n_fitted == 3
      assert fm.r_work() < r_off - 0.005, [params.mode, fm.r_work(), r_off]
      check_parts(fm)
      c = fm.mosaic_correction.data().deep_copy()
      scale(fm, params = params)
      assert approx_equal(fm.mosaic_correction.data(), c)
      check_parts(fm)
  # H contribution present in the data: found, together with the mosaic part
  random.seed(0)
  flex.set_random_seed(0)
  o = tm.synthetic(q = (0., 0.5, 0.8))
  xrs = with_riding_h(o)
  xrs_h = xrs.deep_copy_scatterers()
  xrs_h.set_occupancies(value = 1)
  f_h = o.f_obs.structure_factors_from_scatterers(
    xray_structure = xrs_h.select(xrs_h.hd_selection())).f_calc().data()
  f_atoms = o.f_calc.data()
  f_mask = o.f_whole.deep_copy()
  for qn, f in zip(o.q, o.F): f_mask += f*(qn-1.)
  f_obs = o.f_obs.array(data = flex.abs(f_atoms + f_h + o.k_mask*f_mask))
  for params in [slow_params(), fast_params()]:
    fm = mmtbx.f_model.manager(f_obs = f_obs, r_free_flags = o.r_free_flags,
      xray_structure = xrs, mask_params = o.mask_params)
    scale(fm, params = params, mosaic = True)
    assert flex.max(flex.abs(fm.f_part2().data())) > 0, params.mode
    assert fm.mosaic_result().n_fitted == 3
    check_parts(fm)
    if(params.mode == "slow"):
      assert approx_equal(tm.occupancies(fm.mosaic_result()), o.q, 0.1)
      assert fm.r_work() < 0.03, fm.r_work()

def exercise_back_trace():
  """
  Overall B moved to the atoms (apply_back_trace) with the mosaic model:
  F_model does not change, the mosaic part stays.
  """
  o = tm.synthetic(q = (0., 0.5, 0.8))
  # data with an overall B that the atoms do not have
  f_obs = o.f_obs.array(data = o.f_obs.data()*flex.exp(
    -15./4*o.f_obs.d_star_sq().data()))
  for xrs, with_h in [(o.xray_structure, False), (with_riding_h(o), True)]:
    for params in [slow_params(), fast_params()]:
      def new():
        return mmtbx.f_model.manager(f_obs = f_obs,
          r_free_flags = o.r_free_flags,
          xray_structure = xrs.deep_copy_scatterers(),
          mask_params = o.mask_params)
      fm1 = new()
      scale(fm1, params = params, mosaic = True, apply_back_trace = False)
      fm2 = new()
      scale(fm2, params = params, mosaic = True, apply_back_trace = True)
      b1 = fm1.xray_structure.extract_u_iso_or_u_equiv()
      b2 = fm2.xray_structure.extract_u_iso_or_u_equiv()
      heavy = ~fm2.xray_structure.hd_selection()
      shift = flex.mean((b2-b1).select(heavy))
      assert shift > 0.05, [params.mode, with_h, shift] # 15 A^2 is 0.19 in U
      check_parts(fm2)
      assert fm2.mosaic_result().n_fitted == 3
      # per-bin scales (fast) depend on where the overall B is: looser
      tol_r, tol_q = 5.e-3, 0.05
      if(params.mode == "fast"): tol_r, tol_q = 1.5e-2, 0.15
      assert approx_equal(fm2.r_work(), fm1.r_work(), tol_r), [
        params.mode, with_h, fm2.r_work(), fm1.r_work()]
      assert approx_equal(tm.occupancies(fm2.mosaic_result()),
        tm.occupancies(fm1.mosaic_result()), tol_q), [params.mode, with_h,
        tm.occupancies(fm2.mosaic_result()),
        tm.occupancies(fm1.mosaic_result())]
      # stored correction survives the next call with back trace
      c = fm2.mosaic_correction.data().deep_copy()
      scale(fm2, params = params, apply_back_trace = True)
      assert approx_equal(fm2.mosaic_correction.data(), c)
      check_parts(fm2)

def exercise_anomalous():
  """
  Anomalous data (Friedel mates are separate reflections) through
  update_all_scales.
  """
  o = tm.synthetic(q = (0., 0.5, 0.8))
  f_obs = o.f_obs.generate_bijvoet_mates()
  flags = o.r_free_flags.generate_bijvoet_mates()
  assert f_obs.anomalous_flag()
  def new():
    return mmtbx.f_model.manager(f_obs = f_obs, r_free_flags = flags,
      xray_structure = o.xray_structure.deep_copy_scatterers(),
      mask_params = o.mask_params)
  fm0 = new()
  scale(fm0)
  fm = new()
  scale(fm, mosaic = True)
  assert fm.mosaic_correction.data().size() == f_obs.data().size()
  assert fm.r_work() < fm0.r_work() - 0.005
  assert approx_equal(tm.occupancies(fm.mosaic_result()), o.q, 0.1)
  assert fm.r_work() < 0.01, fm.r_work()
  assert approx_equal(fm.f_model().data(),
    f_model_expected(fm, fm.mosaic_correction.data()), 1.e-6)
  # no anomalous scatterers: the correction obeys Friedel's law
  c = fm.mosaic_correction
  cp, cm = c.hemispheres_acentrics()
  assert cp.indices().size() > 100
  assert approx_equal(cp.data(), flex.conj(cm.data()), 1.e-6)
  # nearly the same occupancies as with merged data (acentric reflections
  # count twice here)
  fm_merged = new_fmodel(o)
  scale(fm_merged, mosaic = True)
  assert approx_equal(tm.occupancies(fm.mosaic_result()),
    tm.occupancies(fm_merged.mosaic_result()), 0.05)

def exercise_optimize_mask(o):
  """
  Mask optimization in the same call changes the mask: regions and
  occupancies must belong to the final mask.
  """
  fm = new_fmodel(o)
  mp = fm.mask_params
  mp.solvent_radius = 1.3 # not the best one: optimization moves it
  mp.shrink_truncation_radius = 1.1
  fm.update(mask_params = mp)
  radii = (fm.mask_params.solvent_radius,
           fm.mask_params.shrink_truncation_radius)
  n0 = mosaic.number_of_runs
  scale(fm, mosaic = True, optimize_mask = True)
  assert (fm.mask_params.solvent_radius,
          fm.mask_params.shrink_truncation_radius) != radii
  assert mosaic.number_of_runs == n0+1
  res = fm.mosaic_result()
  assert res is not None and res.n_fitted == 3
  c = fm.mosaic_correction.data()
  assert approx_equal(fm.f_model().data(), f_model_expected(fm, c), 1.e-6)
  # the correction is the one of the final mask
  assert approx_equal(mosaic.run(fmodel = _flat_copy(fm)).correction, c, 1.e-3)
  # a stored correction does not survive a change of the mask
  fm2 = new_fmodel(o)
  fm2.update(mask_params = mp)
  scale(fm2, mosaic = True)
  assert fm2.mosaic_correction is not None
  scale(fm2, optimize_mask = True)
  assert fm2.mosaic_correction is None
  assert is_zero(fm2.f_part1().data())

def _flat_copy(fm):
  result = fm.deep_copy()
  result.remove_mosaic()
  scale(result)
  return result

def exercise_review_findings(o):
  """
  Cases found by code review.
  """
  # Bulk solvent switched off: nothing to distribute over regions
  fm = new_fmodel(o)
  params = slow_params()
  params.bulk_solvent = False
  n0 = mosaic.number_of_runs
  sio = StringIO()
  scale(fm, params = params, mosaic = True, log = sio)
  assert "mosaic: skipped (no bulk-solvent contribution)" in sio.getvalue()
  assert fm.mosaic_correction is None and fm.mosaic_result() is None
  assert mosaic.number_of_runs == n0
  # Mask parameters changed after the correction was computed: dropped
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  assert fm.mosaic_correction_is_current()
  mp = fm.mask_params
  mp.solvent_radius = 1.4
  fm.update(mask_params = mp)
  assert not fm.mosaic_correction_is_current()
  scale(fm)
  assert fm.mosaic_correction is None and is_zero(fm.f_part1().data())
  # ... and atoms that stop contributing to the mask (occupancy zero)
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  xrs = fm.xray_structure.deep_copy_scatterers()
  sel = flex.bool(xrs.scatterers().size(), False)
  for i in range(5): sel[i] = True
  xrs.set_occupancies(value = 0, selection = sel)
  fm.update_xray_structure(xray_structure = xrs, update_f_calc = True,
    update_f_mask = True)
  assert not fm.mosaic_correction_is_current()
  # Pickle
  from libtbx import easy_pickle
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  fm2 = easy_pickle.loads(easy_pickle.dumps(fm))
  assert approx_equal(fm2.mosaic_correction.data(), fm.mosaic_correction.data())
  assert approx_equal(fm2.f_model().data(), fm.f_model().data())
  assert fm2.mosaic_correction_is_current()
  # object from before the mosaic attributes existed
  state = dict(fm.__dict__)
  for k in ["mosaic_correction", "mosaic_correction_signature",
            "_mosaic_result"]:
    del state[k]
  old = mmtbx.f_model.manager.__new__(mmtbx.f_model.manager)
  old.__setstate__((state,))
  assert old.mosaic_correction is None and old.mosaic_result() is None
  old.select(selection = flex.bool(old.f_obs().data().size(), True))
  # Kick maps are made with the mosaic model (kick_fmodel checks that its
  # copy of the fmodel has the same R as the fmodel)
  from mmtbx.maps import kick
  from cctbx import maptbx
  cg = maptbx.crystal_gridding(unit_cell = fm.f_obs().unit_cell(),
    space_group_info = fm.f_obs().space_group_info(), d_min = 3.5)
  flex.set_random_seed(1); random.seed(1)
  m = kick.kick_fmodel(fmodel = fm, map_type = "mFo-DFc",
    crystal_gridding = cg, number_of_kicks = 1, macro_cycles = 1)
  assert m is not None

def exercise_bulk_and_export(o):
  """
  The mosaic part belongs to the bulk-solvent contribution: f_bulk() has
  it, and the arrays written by export() reproduce F_model.
  """
  fm = new_fmodel(o)
  scale(fm, mosaic = True)
  c = fm.mosaic_correction.data()
  assert flex.max(flex.abs(c)) > 1
  k_mask = fm.k_masks()[0]
  f_mask = fm.f_masks()[0].data()
  assert approx_equal(fm.f_bulk().data(), k_mask*(f_mask + c), 1.e-6)
  k_total = fm.k_isotropic()*fm.k_anisotropic()
  assert approx_equal(fm.f_model().data(),
    k_total*(fm.f_calc().data() + fm.f_bulk().data()), 1.e-6)
  # flat model: as before
  fm0 = new_fmodel(o)
  scale(fm0)
  assert approx_equal(fm0.f_bulk().data(),
    fm0.k_masks()[0]*fm0.f_masks()[0].data(), 1.e-6)
  # exported file
  from iotbx import reflection_file_reader
  for fmodel, name in [(fm, "tst_fmodel_mosaic_export.mtz"),
                       (fm0, "tst_fmodel_mosaic_export_flat.mtz")]:
    with open(name, "w") as f: fmodel.export(out = f)
    arrays = {}
    for ma in reflection_file_reader.any_reflection_file(
        file_name = name).as_miller_arrays():
      arrays[ma.info().labels[0]] = ma
    ref = arrays["FMODEL"]
    def get(label):
      a, b = ref.common_sets(arrays[label])
      assert a.indices().size() == ref.indices().size()
      return b.data()
    recomputed = get("K_ISOTROPIC")*get("K_ANISOTROPIC")*(
      get("FCALC") + get("K_MASK")*get("FMASK"))
    d1, d2 = flex.abs(recomputed), flex.abs(ref.data())
    k = flex.sum(d1*d2)/flex.sum(d1*d1) # FMODEL is on the scale of F_obs
    r = flex.sum(flex.abs(d1*k-d2))/flex.sum(d2)
    assert r < 1.e-3, r

def exercise_fmodels(o):
  """
  mmtbx.fmodels passes the keyword to the X-ray fmodel.
  """
  import mmtbx
  from libtbx.utils import null_out
  fmodels = mmtbx.fmodels(fmodel_xray = new_fmodel(o), log = null_out())
  n0 = mosaic.number_of_runs
  fmodels.update_all_scales(update_f_part1 = False, remove_outliers = False,
    params = slow_params())
  assert fmodels.fmodel_xray().mosaic_correction is None
  r_off = fmodels.fmodel_xray().r_work()
  assert mosaic.number_of_runs == n0
  fmodels.update_all_scales(update_f_part1 = False, remove_outliers = False,
    params = slow_params(), mosaic = True)
  assert mosaic.number_of_runs == n0+1
  assert fmodels.fmodel_xray().mosaic_result().n_fitted == 3
  assert fmodels.fmodel_xray().r_work() < r_off - 0.005
  # kept by the next call, which recomputes F_calc and F_mask first
  fmodels.update_all_scales(update_f_part1 = False, remove_outliers = False,
    params = slow_params())
  assert mosaic.number_of_runs == n0+1
  assert fmodels.fmodel_xray().mosaic_correction is not None
  assert fmodels.fmodel_xray().r_work() < r_off - 0.005

def run():
  random.seed(0)
  flex.set_random_seed(0)
  o = tm.synthetic(q = (0., 0.5, 0.8))
  r_off = exercise_off(o)
  fm_on = exercise_on(o, r_off)
  exercise_reuse(o, fm_on)
  exercise_select(o, fm_on)
  exercise_atoms_change(o, fm_on)
  exercise_maps()
  exercise_skip(o)
  exercise_bulk_and_export(o)
  exercise_fmodels(o)
  exercise_review_findings(o)
  exercise_optimize_mask(o)
  exercise_hydrogens()
  exercise_back_trace()
  exercise_anomalous()

if (__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %.2f"%(time.time()-t0))
  print("OK")

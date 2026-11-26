from cctbx.array_family import flex
from mmtbx.bulk_solvent import squeeze
from mmtbx.bulk_solvent import bulk_solvent_and_scaling as bss
from mmtbx.bulk_solvent import tst_squeeze as ts
from libtbx.test_utils import approx_equal
from six.moves import cStringIO as StringIO
import mmtbx.f_model
import random, time

"""
Squeeze inside mmtbx.f_model.manager: update_all_scales(squeeze=True).
Simulated data: man-made blobs of density of different size in the solvent
region that the model does not have must disappear from the difference map.
"""

def scale(fm, **kwds):
  args = dict(remove_outliers = False, apply_scale_k1_to_f_obs = False)
  args.update(kwds)
  return fm.update_all_scales(**args)

def is_zero(a):
  return flex.abs(a).all_eq(0)

def f_model_expected(fm):
  """k_total*(F_calc + k_mask*F_mask + F_mosaic + F_squeeze), no H here."""
  k_total = fm.k_isotropic()*fm.k_anisotropic()
  return k_total*(fm.f_calc().data() +
    fm.k_masks()[0]*fm.f_masks()[0].data() +
    fm.f_mosaic().data() + fm.f_squeeze().data())

def difference_map(fm):
  gridding = squeeze.crystal_gridding(fmodel = fm,
    params = squeeze.default_params())
  return squeeze.difference_map(fmodel = fm, gridding = gridding)

def exercise_off(o):
  """
  Without squeeze=True nothing changes.
  """
  fm = ts.new_fmodel(o, scale = False)
  assert fm.squeeze_result() is None
  assert is_zero(fm.f_squeeze().data()) and is_zero(fm.f_mosaic().data())
  scale(fm)
  assert fm.squeeze_result() is None
  assert is_zero(fm.f_squeeze().data()) and is_zero(fm.f_part1().data())
  assert fm.f_squeeze().indices().all_eq(fm.f_obs().indices())
  assert approx_equal(fm.f_model().data(), f_model_expected(fm))
  fm2 = ts.new_fmodel(o, scale = False)
  scale(fm2, squeeze = False, water_selection = o.water_selection)
  assert approx_equal(fm2.f_model().data(), fm.f_model().data(), 1.e-9)
  fm2.remove_squeeze()
  scale(fm2)
  assert approx_equal(fm2.f_model().data(), fm.f_model().data(), 1.e-9)
  return fm

def exercise_blobs_disappear(o, fm_off):
  """
  The point of it all: blobs of different size are in the difference map
  before and gone after; the map elsewhere is cleaner; R-factors drop.
  """
  m0 = difference_map(fm_off)
  sigma0 = m0.sample_standard_deviation()
  before = [ts.density_at(m0, s)/sigma0 for s in o.blob_sites_frac]
  for d in before: assert d > 5, before
  fm = ts.new_fmodel(o, scale = False)
  sio = StringIO()
  scale(fm, squeeze = True, water_selection = o.water_selection, log = sio)
  assert "squeeze: density above 2.5 sigma" in sio.getvalue()
  res = fm.squeeze_result()
  assert res is not None and res.k > 1 and res.water_excluded
  assert fm.squeeze_is_current()
  m1 = difference_map(fm)
  # in units of the sigma of the map before
  after = [ts.density_at(m1, s)/sigma0 for s in o.blob_sites_frac]
  for b, a in zip(before, after):
    assert a < 4 and a < b/5, [before, after]
  # noise of the map goes down, also away from the blobs
  assert m1.sample_standard_deviation() < 0.8*sigma0
  def rms_at(m, sites):
    v = flex.double([m.tricubic_interpolation(s) for s in sites])
    return flex.mean_sq(v)**0.5
  assert rms_at(m1, o.control_sites_frac) < rms_at(m0, o.control_sites_frac)
  # R-factors
  assert fm.r_work() < 0.8*fm_off.r_work(), [fm.r_work(), fm_off.r_work()]
  assert fm.r_free() < 0.8*fm_off.r_free(), [fm.r_free(), fm_off.r_free()]
  # F_model is what it is documented to be
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
  assert approx_equal(fm.f_squeeze().data(), res.scaled(ss = fm.ss))
  assert approx_equal(fm.f_part1().data(), fm.f_squeeze().data())
  assert is_zero(fm.f_mosaic().data())
  assert ts.correlation(fm.f_squeeze().data(), o.f_blobs) > 0.85
  # bulk solvent does not include it
  assert approx_equal(fm.f_bulk().data(),
    fm.k_masks()[0]*fm.f_masks()[0].data())
  return fm

def exercise_reuse(o, fm_off, fm_on):
  """
  Calls without squeeze=True keep the stored contribution; remove_squeeze()
  goes back.
  """
  r_on = fm_on.r_work()
  f = fm_on.squeeze_result().f.deep_copy()
  for fm in [fm_on.deep_copy(), fm_on]:
    assert fm.squeeze_result() is not None
    for i in range(2):
      scale(fm)
      assert approx_equal(fm.r_work(), r_on, 1.e-4)
      assert approx_equal(fm.squeeze_result().f, f)
      assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
      assert not is_zero(fm.f_squeeze().data())
  # slow scaling
  fm = fm_on.deep_copy()
  params = bss.master_params.extract()
  params.mode = "slow"
  scale(fm, params = params)
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
  assert fm.r_work() < 0.8*fm_off.r_work()
  # overall B moved to atoms: F_model is the same, B of squeeze follows
  fm = fm_on.deep_copy()
  b0 = fm.squeeze_result().b
  b_atoms0 = flex.mean(fm.xray_structure.extract_u_iso_or_u_equiv())
  scale(fm, apply_back_trace = True)
  b_shift = (flex.mean(fm.xray_structure.extract_u_iso_or_u_equiv()) -
    b_atoms0)*79.
  assert approx_equal(fm.squeeze_result().b - b0, b_shift, 0.01)
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
  assert approx_equal(fm.r_work(), r_on, 1.e-3)
  # computed again on request: about the same
  fm = fm_on.deep_copy()
  scale(fm, squeeze = True, water_selection = o.water_selection)
  assert approx_equal(fm.r_work(), r_on, 1.e-3)
  assert approx_equal(fm.squeeze_result().k, fm_on.squeeze_result().k, 0.05)
  # atoms moved and the mask computed again: dropped
  fm = fm_on.deep_copy()
  xrs = fm.xray_structure.deep_copy_scatterers()
  xrs.shake_sites_in_place(mean_distance = 0.3)
  fm.update_xray_structure(xray_structure = xrs, update_f_calc = True,
    update_f_mask = True)
  assert not fm.squeeze_is_current()
  sio = StringIO()
  scale(fm, log = sio)
  assert "squeeze: stored contribution dropped" in sio.getvalue()
  assert fm.squeeze_result() is None and is_zero(fm.f_squeeze().data())
  assert is_zero(fm.f_part1().data())
  # back to the model without it
  fm = fm_on.deep_copy()
  fm.remove_squeeze()
  assert fm.squeeze_result() is None
  scale(fm)
  assert is_zero(fm.f_part1().data())
  assert approx_equal(fm.r_work(), fm_off.r_work(), 1.e-6)
  assert approx_equal(fm.f_model().data(), fm_off.f_model().data(), 1.e-6)
  # fm_on itself still has it
  assert approx_equal(fm_on.r_work(), r_on, 1.e-4)

def exercise_select(o, fm_on):
  f = fm_on.f_squeeze()
  sel = flex.random_bool(f.data().size(), 0.7)
  fm = fm_on.select(selection = sel)
  assert fm.squeeze_result().f.size() == sel.count(True)
  assert approx_equal(fm.f_squeeze().data(), f.data().select(sel))
  assert fm.f_squeeze().indices().all_eq(fm.f_obs().indices())
  assert approx_equal(fm.f_model().data(),
    fm_on.f_model().data().select(sel))
  assert fm.squeeze_is_current()
  scale(fm)
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
  fm = fm_on.deep_copy()
  fm.select(selection = sel, in_place = True)
  assert fm.squeeze_result().f.size() == sel.count(True)
  # outlier rejection during scaling
  data = o.f_obs.data().deep_copy()
  i_max = flex.max_index(data)
  data[i_max] = data[i_max]*50
  fm = mmtbx.f_model.manager(
    f_obs          = o.f_obs.array(data = data),
    r_free_flags   = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers())
  scale(fm, squeeze = True, water_selection = o.water_selection,
    remove_outliers = True)
  n = fm.f_obs().data().size()
  assert n < data.size()
  assert fm.squeeze_result().f.size() == n
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)

def exercise_with_mosaic(o, fm_on):
  """
  Mosaic bulk solvent first, then squeeze; both are in F_model.
  """
  fm = ts.new_fmodel(o, scale = False)
  sio = StringIO()
  scale(fm, mosaic = True, squeeze = True,
    water_selection = o.water_selection, log = sio)
  assert fm.squeeze_result() is not None
  assert approx_equal(fm.f_model().data(), f_model_expected(fm), 1.e-6)
  assert approx_equal(fm.f_part1().data(),
    fm.f_mosaic().data() + fm.f_squeeze().data())
  assert approx_equal(fm.f_bulk().data(),
    fm.k_masks()[0]*fm.f_masks()[0].data() + fm.f_mosaic().data())
  assert approx_equal(fm.r_work(), fm_on.r_work(), 0.005)

def exercise_skip(o):
  """
  Not applicable: skipped with a message, F_model as without squeeze.
  """
  # resolution
  fm = ts.new_fmodel(o, scale = False).resolution_filter(d_min = 1.5)
  ref = fm.deep_copy()
  scale(ref)
  sio = StringIO()
  scale(fm, squeeze = True, log = sio)
  assert "squeeze: skipped (data resolution 1.50 A" in sio.getvalue()
  assert fm.squeeze_result() is None
  assert approx_equal(fm.f_model().data(), ref.f_model().data(), 1.e-9)
  # ... unless allowed
  params = bss.master_params.extract()
  params.squeeze.d_min_max = None
  scale(fm, squeeze = True, params = params)
  assert fm.squeeze_result() is not None
  assert fm.r_work() < ref.r_work()
  # settings come from params.squeeze
  params.squeeze.threshold = 1000
  sio = StringIO()
  scale(fm, squeeze = True, params = params, log = sio)
  assert "no density above" in sio.getvalue()
  assert fm.squeeze_result() is None
  assert approx_equal(fm.f_model().data(), ref.f_model().data(), 1.e-6)
  # twinned
  fm = mmtbx.f_model.manager(
    f_obs          = o.f_obs,
    r_free_flags   = o.r_free_flags,
    xray_structure = o.xray_structure.deep_copy_scatterers(),
    twin_law       = "k,h,-l")
  sio = StringIO()
  fm.update_all_scales(squeeze = True, remove_outliers = False, log = sio)
  assert "squeeze: skipped (twinned data)" in sio.getvalue()
  assert fm.squeeze_result() is None

def run():
  random.seed(0)
  flex.set_random_seed(0)
  o = ts.synthetic()
  fm_off = exercise_off(o)
  fm_on = exercise_blobs_disappear(o, fm_off)
  exercise_reuse(o, fm_off, fm_on)
  exercise_select(o, fm_on)
  exercise_with_mosaic(o, fm_on)
  exercise_skip(o)

if(__name__ == "__main__"):
  t0 = time.time()
  run()
  print("Time: %6.2f"%(time.time()-t0))
  print("OK")

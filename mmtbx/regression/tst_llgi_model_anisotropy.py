from __future__ import absolute_import, division, print_function
import math
import numpy as np
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_model_anisotropy as llgi_model_anisotropy
import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac

def miller_set(cell, symbol, d_min):
  cs = crystal.symmetry(unit_cell=cell, space_group_symbol=symbol)
  return cs.build_miller_set(anomalous_flag=False, d_min=d_min)

def exercise_tensor_basis():
  # trace-free tensors invariant under the point group: 5 (triclinic),
  # 3 (monoclinic), 2 (orthorhombic), 1 (tetragonal, trigonal, hexagonal),
  # 0 (cubic)
  for cell, symbol, n in [
      ((40, 50, 60, 80, 95, 100), "P 1", 5),
      ((51, 62, 73, 90, 97, 90), "C 2", 3),
      ((41, 47, 283, 90, 90, 90), "P 2 2 21", 2),
      ((60, 60, 90, 90, 90, 90), "P 41 21 2", 1),
      ((172, 172, 277, 90, 90, 120), "P 63 2 2", 1),
      ((80, 80, 80, 90, 90, 90), "P 21 3", 0)]:
    ms = miller_set(cell, symbol, 4.0)
    basis = llgi_model_anisotropy.tensor_basis(ms.unit_cell(), ms.space_group())
    assert len(basis) == n, (symbol, len(basis))
    for t in basis:
      assert abs(np.trace(t)) < 1e-12 and np.allclose(t, t.T)
    # the factor is the same for symmetry-equivalent reflections
    if(n > 0):
      f = llgi_model_anisotropy.anisotropic_factor(ms.unit_cell(), ms.space_group(),
        "resolution_dependent", d_star_sq_max=flex.max(ms.d_star_sq().data()),
        b=np.linspace(-20, 30, 2 * n))
      p1 = ms.expand_to_p1()
      a_p1 = f.evaluate(p1.indices())
      a_asu = f.evaluate(p1.map_to_asu().indices())
      assert approx_equal(a_p1, a_asu, eps=1e-10)

def wilson_intensities(ms, true_factor, seed=0):
  # intensities drawn from the Wilson distributions about
  # Sigma = iso(s) * A(h) (epsilon included in the mean)
  rng = np.random.RandomState(seed)
  dss = ms.d_star_sq().data().as_numpy_array()
  iso = np.exp(-20 * dss / 4) * (1 + 0.3 * np.cos(6 * dss))
  a = np.array(true_factor.evaluate(ms.indices())) if true_factor else np.ones(ms.size())
  mean = iso * a
  cen = ms.centric_flags().data().as_numpy_array()
  i = np.where(cen, mean * rng.normal(size=ms.size()) ** 2, rng.exponential(mean))
  return flex.double(i)

def binned_iso_fn(ms, n_bins=20):
  dss = ms.d_star_sq().data().as_numpy_array()
  order = np.argsort(dss)
  bins = np.array_split(order, n_bins)
  d_mean = np.array([dss[c].mean() for c in bins])
  def iso_fn(v):
    v = np.array(v)
    return flex.double(np.interp(dss, d_mean, np.array([v[c].mean() for c in bins])))
  return iso_fn

def exercise_fit_recovers_factor():
  # a planted resolution-dependent anisotropy is recovered from model-like
  # intensities; isotropic intensities give a factor close to 1; off gives 1
  ms = miller_set((51, 62, 73, 90, 97, 90), "C 2", 2.2)
  dss_max = flex.max(ms.d_star_sq().data())
  true = llgi_model_anisotropy.anisotropic_factor(ms.unit_cell(), ms.space_group(),
    "resolution_dependent", d_star_sq_max=dss_max, b=[30, -20, 10, -40, 25, 15])
  i = wilson_intensities(ms, true, seed=1)
  r = llgi_model_anisotropy.fit(i, ms.indices(), ms.centric_flags().data(),
    ms.unit_cell(), ms.space_group(), binned_iso_fn(ms), "resolution_dependent")
  assert np.allclose(r.factor.b, true.b, atol=4.0), r.factor.b
  rel = np.array(r.a) / np.array(true.evaluate(ms.indices()))
  assert np.std(np.log(rel)) < 0.03, np.std(np.log(rel))
  # normalised intensities: unit mean
  assert approx_equal(flex.mean(i / r.sigma), 1, eps=0.02)
  i_iso = wilson_intensities(ms, None, seed=2)
  r = llgi_model_anisotropy.fit(i_iso, ms.indices(), ms.centric_flags().data(),
    ms.unit_cell(), ms.space_group(), binned_iso_fn(ms), "ellipsoid")
  assert np.max(np.abs(np.log(np.array(r.a)))) < 0.05
  r = llgi_model_anisotropy.fit(i_iso, ms.indices(), ms.centric_flags().data(),
    ms.unit_cell(), ms.space_group(), binned_iso_fn(ms), "off")
  assert r.factor is None and r.a.all_eq(1)
  # cubic: nothing to fit
  mc = miller_set((80, 80, 80, 90, 90, 90), "P 21 3", 3.0)
  r = llgi_model_anisotropy.fit(wilson_intensities(mc, None), mc.indices(),
    mc.centric_flags().data(), mc.unit_cell(), mc.space_group(), binned_iso_fn(mc),
    "resolution_dependent")
  assert r.factor.n_parameters() == 0 and r.a.all_eq(1)

def exercise_build_e_model_and_scatfrac():
  # off is the isotropic normalisation as before; with model_anisotropy the
  # model's E (and the ScatFrac-normalised intensities) have unit mean
  # square in each direction
  ms = miller_set((41, 47, 283, 90, 90, 90), "P 2 2 21", 2.5)
  true = llgi_model_anisotropy.anisotropic_factor(ms.unit_cell(), ms.space_group(),
    "ellipsoid", d_star_sq_max=flex.max(ms.d_star_sq().data()), b=[60, -35])
  rng = np.random.RandomState(4)
  i = np.array(wilson_intensities(ms, true, seed=3))
  eps = ms.epsilons().data().as_double()
  f = flex.polar(flex.double(np.sqrt(i * np.array(eps))),
    flex.double(rng.uniform(0, 2 * math.pi, ms.size())))
  dss = ms.d_star_sq().data()
  old = llgi_e_sigmaa.build_sigma_p(f, dss)
  em0 = llgi_e_sigmaa.build_e_model(f, eps, dss)
  assert approx_equal(em0.sigma_p, old, eps=0) and em0.factor is None
  em1 = llgi_e_sigmaa.build_e_model(f, eps, dss, miller_set=ms,
    model_anisotropy="off")
  assert approx_equal(em1.e_model, em0.e_model, eps=0)
  em = llgi_e_sigmaa.build_e_model(f, eps, dss, miller_set=ms,
    model_anisotropy="ellipsoid")
  a_true = np.array(true.evaluate(ms.indices()))
  e2 = np.abs(np.array(em.e_model)) ** 2
  e2_iso = np.abs(np.array(em0.e_model)) ** 2
  weak, strong = a_true < np.percentile(a_true, 20), a_true > np.percentile(a_true, 80)
  assert abs(e2[weak].mean() - 1) < 0.05 and abs(e2[strong].mean() - 1) < 0.05
  assert e2_iso[weak].mean() < 0.8 and e2_iso[strong].mean() > 1.2
  # sigma_p_at evaluates the isotropic part elsewhere; the factor too
  assert approx_equal(em.sigma_p_at(dss), em.sigma_p, eps=1e-10)
  assert approx_equal(em.factor.evaluate(ms.indices()), em.anisotropy, eps=1e-10)
  # ScatFrac: off as before; anisotropic, x/ScatFrac has unit mean by direction
  teps, resn = flex.double(ms.size(), 1), flex.double(ms.size(), 1)
  s0 = llgi_scatfrac.scatfrac_moments(f, teps, resn, dss)
  s1 = llgi_scatfrac.scatfrac_moments(f, teps, resn, dss, miller_set=ms,
    model_anisotropy="off")
  assert approx_equal(s0, s1, eps=0)
  s = llgi_scatfrac.scatfrac_moments(f, teps, resn, dss, miller_set=ms,
    model_anisotropy="ellipsoid")
  x = np.abs(np.array(f)) ** 2 / np.array(s)
  assert abs(x[weak].mean() / x.mean() - 1) < 0.05 and abs(x[strong].mean() / x.mean() - 1) < 0.05

def exercise_resolve():
  # auto takes the family of nacelle's ANISOBETA: resolution_dependent if
  # it changes with resolution, otherwise (also without it) ellipsoid
  from libtbx import group_args
  ms = miller_set((51, 62, 73, 90, 97, 90), "C 2", 2.5)
  dss_max = flex.max(ms.d_star_sq().data())
  def llgi_data(b, mode="ellipsoid", **kw):
    f = llgi_model_anisotropy.anisotropic_factor(ms.unit_cell(),
      ms.space_group(), mode, d_star_sq_max=dss_max, b=b)
    return group_args(anisobeta=ms.array(data=f.evaluate(ms.indices())), **kw)
  ell = llgi_data([30, -20, 10])
  rd = llgi_data([30, -20, 10, 3, -2, 1], "resolution_dependent")
  flat = group_args(anisobeta=ms.array(data=flex.double(ms.size(), 1)))
  resolve = llgi_model_anisotropy.resolve
  assert resolve("auto", ell) == "ellipsoid"
  assert resolve("auto", rd) == "resolution_dependent"
  assert resolve("auto", flat) == "ellipsoid"
  assert resolve("auto", None) == "ellipsoid"
  assert resolve("auto", group_args(anisobeta=None)) == "ellipsoid"
  for mode in llgi_model_anisotropy.MODES:
    assert resolve(mode, rd) == mode
  # mode_for: the sigmaA fit's setting (the default, auto, if e_scale is
  # None), resolved; off without a fit
  assert llgi_model_anisotropy.mode_for(rd) == "off"
  fp = group_args(scatfrac=None, e_scale=None)
  assert llgi_model_anisotropy.mode_for(llgi_data([3, -2, 1, 3, -2, 1],
    "resolution_dependent", fit_params=fp)) == "resolution_dependent"
  fp = group_args(scatfrac=None, e_scale=group_args(model_anisotropy="off"))
  assert llgi_model_anisotropy.mode_for(llgi_data([3, -2, 1, 3, -2, 1],
    "resolution_dependent", fit_params=fp)) == "off"

def run():
  exercise_tensor_basis()
  exercise_fit_recovers_factor()
  exercise_build_e_model_and_scatfrac()
  exercise_resolve()
  print("OK")

if (__name__ == "__main__"):
  run()

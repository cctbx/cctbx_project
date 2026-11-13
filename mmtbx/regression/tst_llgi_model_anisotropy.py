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
      ((60, 60, 60, 80, 80, 80), "R 3 :R", 1),
      ((75, 75, 140, 90, 90, 120), "R 3 2 :H", 1),
      ((80, 80, 80, 90, 90, 90), "P 21 3", 0)]:
    ms = miller_set(cell, symbol, 4.0)
    basis = llgi_model_anisotropy.tensor_basis(ms.unit_cell(), ms.space_group())
    assert len(basis) == n, (symbol, len(basis))
    # the symmetry-allowed ADPs, less the isotropic one
    assert n == ms.space_group().adp_constraints().n_independent_params() - 1
    for i, t in enumerate(basis):
      assert abs(np.trace(t)) < 1e-12 and np.allclose(t, t.T)
      for j, t2 in enumerate(basis):
        assert abs(np.sum(t * t2) - (i == j)) < 1e-12
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
  # off is the isotropic normalisation; with model_anisotropy the model's E
  # (and the ScatFrac-normalised intensities) have unit mean square in each
  # direction
  ms = miller_set((41, 47, 283, 90, 90, 90), "P 2 2 21", 2.5)
  true = llgi_model_anisotropy.anisotropic_factor(ms.unit_cell(), ms.space_group(),
    "ellipsoid", d_star_sq_max=flex.max(ms.d_star_sq().data()), b=[60, -35])
  rng = np.random.RandomState(4)
  i = np.array(wilson_intensities(ms, true, seed=3))
  eps = ms.epsilons().data().as_double()
  f = flex.polar(flex.double(np.sqrt(i * np.array(eps))),
    flex.double(rng.uniform(0, 2 * math.pi, ms.size())))
  dss = ms.d_star_sq().data()
  one = flex.double(ms.size(), 1)
  def e_model(mode):   # RESN = sqrt(EPS), TEPS = k_anisotropic = 1
    return llgi_e_sigmaa.build_e_model(f, one, ms, flex.sqrt(eps), one,
      model_anisotropy=mode)
  em0 = e_model("off")
  em = e_model("ellipsoid")
  a_true = np.array(true.evaluate(ms.indices()))
  e2 = np.abs(np.array(em.e_model)) ** 2
  e2_iso = np.abs(np.array(em0.e_model)) ** 2
  weak, strong = a_true < np.percentile(a_true, 20), a_true > np.percentile(a_true, 80)
  assert abs(e2[weak].mean() - 1) < 0.05 and abs(e2[strong].mean() - 1) < 0.05
  assert e2_iso[weak].mean() < 0.8 and e2_iso[strong].mean() > 1.2
  # e_abs_for evaluates the normalisation for other reflections (here the
  # same ones)
  assert approx_equal(em.e_abs_for(ms, f, flex.sqrt(eps), one, one),
    flex.abs(em.e_model), eps=1e-10)
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

def exercise_safeguards():
  # too few reflections per parameter: A = 1; "auto" must be resolved first
  ms = miller_set((51, 62, 73, 90, 97, 90), "C 2", 8.0)
  assert ms.size() < llgi_model_anisotropy.MIN_REFLECTIONS_PER_PARAMETER * 3
  r = llgi_model_anisotropy.fit(wilson_intensities(ms, None), ms.indices(),
    ms.centric_flags().data(), ms.unit_cell(), ms.space_group(),
    binned_iso_fn(ms), "ellipsoid")
  assert r.a.all_eq(1) and np.all(r.factor.b == 0)
  f = flex.complex_double(ms.size(), 1)
  for call in (
      lambda: llgi_e_sigmaa.build_e_model(f, flex.double(ms.size(), 1), ms,
        flex.double(ms.size(), 1), flex.double(ms.size(), 1),
        model_anisotropy="auto"),
      lambda: llgi_scatfrac.scatfrac_moments(f, flex.double(ms.size(), 1),
        flex.double(ms.size(), 1), ms.d_star_sq().data(), miller_set=ms,
        model_anisotropy="auto")):
    try: call()
    except ValueError: pass
    else: raise AssertionError("auto accepted")
  # ScatFrac's likelihood form stays isotropic: a warning
  import mmtbx.refinement.llgi_scatfrac as ls
  p = ls.llgi_scatfrac_params.extract()
  assert llgi_e_sigmaa.isotropic_scatfrac_warnings(p, "ellipsoid") == []
  p.method = "likelihood"
  assert len(llgi_e_sigmaa.isotropic_scatfrac_warnings(p, "ellipsoid")) == 1
  assert llgi_e_sigmaa.isotropic_scatfrac_warnings(p, "off") == []

def exercise_missing_reflections():
  # the maps' fill normalises the model as the map coefficients (and the
  # sigmaA fit) do: the same |Emodel| for the observed reflections
  from mmtbx import map_tools as mt
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  fmodel = build_llgi_fmodel(n_atoms=40, d_min=2.2, seed=3)
  llgi_data = fmodel.llgi_data()
  assert llgi_model_anisotropy.mode_for(llgi_data) == "ellipsoid"
  mro = mt.model_missing_reflections_llgi(fmodel=fmodel)
  mro.get_missing()
  m = mro.e_scale_missing
  f_obs = fmodel.f_obs()
  em = llgi_e_sigmaa.build_e_model_as_fitted(fmodel)
  from cctbx import miller as cmiller
  pairs = cmiller.match_indices(m.miller_set.indices(), f_obs.indices()).pairs()
  assert pairs.size() == f_obs.size()
  assert approx_equal(m.e_model_abs.select(pairs.column(0)),
    flex.abs(em.e_model).select(pairs.column(1)), eps=1e-8)

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
  # fill values (1) for reflections without nacelle's value, and values
  # <= 0: nacelle's own values (nacelle_all) are used if present, values
  # <= 0 are ignored
  filled = ell.anisobeta.deep_copy()
  filled.data().set_selected(flex.size_t(range(0, filled.size(), 50)), 1.0)
  assert resolve("auto", group_args(anisobeta=filled)) == "resolution_dependent"
  assert resolve("auto", group_args(anisobeta=filled,
    nacelle_all=group_args(anisobeta=ell.anisobeta))) == "ellipsoid"
  zeros = ell.anisobeta.deep_copy()
  zeros.data().set_selected(flex.size_t(range(0, zeros.size(), 50)), 0.0)
  assert resolve("auto", group_args(anisobeta=zeros)) == "ellipsoid"
  # kept on llgi_data while the array is the same
  assert resolve("auto", rd) == "resolution_dependent"
  assert rd._model_anisotropy_auto[0] is rd.anisobeta
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
  exercise_safeguards()
  exercise_missing_reflections()
  exercise_resolve()
  print("OK")

if (__name__ == "__main__"):
  run()

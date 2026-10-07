from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from cctbx.xray import ext
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_hybrid as llgi_hybrid
import numpy as np
import random

from mmtbx.regression.tst_llgi_map_coefficients import (
  build_fmodel, synthetic_llgi_data)

def add_hybrid_data(fmodel, llgi_data, rice_kappa, seed=7):
  """ Synthetic intensities consistent with the synthetic Feff/RESN:
  E_obs^2 = Eeff^2 + noise, sigma(E_obs^2) in [0.2, 2]. """
  rnd = random.Random(seed)
  f_obs = fmodel.f_obs()
  eeff = llgi_data.feff.data() / llgi_data.resn.data()
  n = eeff.size()
  sig = flex.double([0.2 + 1.8 * rnd.random() for i in range(n)])
  e2 = eeff * eeff + sig * flex.double([rnd.gauss(0, 1) for i in range(n)])
  params = llgi_hybrid.llgi_hybrid_params.extract()
  params.rice_kappa = rice_kappa
  flags = llgi_hybrid.force_exact_flags(e2, sig,
    f_obs.centric_flags().data())
  return llgi_hybrid.replace_llgi_data(llgi_data,
    e_obs_sq=f_obs.array(data=e2), sig_e_obs_sq=f_obs.array(data=sig),
    force_exact=f_obs.array(data=flags.force_exact), hybrid_params=params)

def build_hybrid_fmodel(rice_kappa, n_atoms=40, d_min=2.2, seed=3):
  fmodel = build_fmodel(n_atoms=n_atoms, d_min=d_min, seed=seed)
  llgi_data = synthetic_llgi_data(fmodel, seed=seed + 100)
  llgi_data = add_hybrid_data(fmodel, llgi_data, rice_kappa)
  fmodel.set_llgi_data(llgi_data)
  fmodel.set_target_name("llgi")
  fmodel.update_llgi_sigmaa_scatfrac()
  return fmodel

def exercise_helpers_and_select():
  fmodel = build_hybrid_fmodel(rice_kappa=0.0)
  llgi_data = fmodel.llgi_data()
  h = llgi_hybrid.get_hybrid(llgi_data)
  assert h is not None and h.size() == fmodel.f_obs().size()
  assert llgi_hybrid.get_hybrid(llgi_data) is h  # cached
  # replace keeps the cache, select rebuilds on the new index set
  same = llgi_hybrid.replace_llgi_data(llgi_data, sigmaa=llgi_data.sigmaa)
  assert llgi_hybrid.get_hybrid(same) is h
  sel = flex.bool([i % 3 != 0 for i in range(fmodel.f_obs().size())])
  sub = fmodel.select(sel)
  hs = llgi_hybrid.get_hybrid(sub.llgi_data())
  assert hs.size() == sel.count(True)
  assert approx_equal(hs.e_obs_sq, h.e_obs_sq.select(sel))
  assert approx_equal(hs.null_log_z, h.null_log_z.select(sel))
  # disabled hybrid
  p = llgi_data.hybrid_params
  p.enabled = False
  assert llgi_hybrid.get_hybrid(
    llgi_hybrid.replace_llgi_data(llgi_data, hybrid_params=p)) is None
  p.enabled = True

def exercise_map_coefficients_exact_branch():
  """ With every reflection exact, the difference coefficient
  (m*F - D*Fc on the E scale) must be <E> - sigmaA*Emodel, proportional
  to the exact LLGI gradient with respect to Emodel. """
  import mmtbx.refinement.llgi_e_bulk_solvent as llgi_e_bs
  fmodel = build_hybrid_fmodel(rice_kappa=0.0)
  mch = fmodel.map_calculation_helper_llgi()
  assert mch.n_exact > 0
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  fmnas = llgi_e_bs.f_model_no_aniso_scale(fmodel)
  em = llgi_e_bs.build_e_model(fmnas.data(),
    f_obs.epsilons().data().as_double(), f_obs.d_star_sq().data())
  emodel_abs = flex.abs(em.e_model)
  sa = llgi_data.sigmaa.data()
  resn = llgi_data.resn.data()
  centric = f_obs.centric_flags().data()
  h = llgi_hybrid.get_hybrid(llgi_data)
  work = ~fmodel.r_free_flags().data()
  r = ext.llgi_e_emodel_target_and_gradients(
    e_eff=llgi_data.feff.data() / resn, selection=work, e_model=emodel_abs,
    dobs=llgi_data.dobs.data(), sigmaa=sa, centric_flags=centric, hybrid=h)
  assert r.n_exact() == work.count(True)
  g = r.d_target_by_demodel()
  n_work = work.count(True)
  inv = 1 / flex.sqrt(em.sigma_p * f_obs.epsilons().data().as_double())
  n_checked = 0
  for i in range(f_obs.size()):
    if(not work[i] or sa[i] <= 0): continue
    e_expected = mch.fom[i] * mch.f_obs.data()[i] / resn[i]
    d_emodel = mch.alpha.data()[i] / (resn[i] * inv[i])
    assert approx_equal(d_emodel, sa[i], eps=1e-10)
    b = 0.5 if centric[i] else 1.0
    s = 1 - sa[i]**2
    expected = -(2 * b * sa[i] / s) * (e_expected - sa[i] * emodel_abs[i])
    # d/dEmodel of the exact LLGI: d/dEc with Ec = Emodel
    assert approx_equal(g[i] * n_work, expected, eps=1e-6), (i,)
    assert 0 <= mch.fom[i] < 1
    n_checked += 1
  assert n_checked > 0
  # Phase errors remain well defined
  pe = fmodel.phase_errors_llgi(mch)
  assert flex.min(pe) >= 0 and flex.max(pe) <= 90.0001

def exercise_map_coefficients_exact_without_hybrid():
  """ The map coefficients use the exact posterior for every measured
  reflection whether or not the hybrid target is enabled: with it
  disabled, they match those of an all-exact hybrid. """
  fmodel = build_hybrid_fmodel(rice_kappa=0.0)
  m_on = fmodel.map_calculation_helper_llgi()
  llgi_data = fmodel.llgi_data()
  p = llgi_data.hybrid_params
  p.enabled = False
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(llgi_data,
    hybrid_params=p))
  assert llgi_hybrid.get_hybrid(fmodel.llgi_data()) is None
  m_off = fmodel.map_calculation_helper_llgi()
  p.enabled = True
  assert m_off.n_exact > 0 and m_off.n_exact == m_on.n_exact
  assert approx_equal(m_off.fom, m_on.fom)
  assert approx_equal(m_off.f_obs.data(), m_on.f_obs.data())
  assert approx_equal(m_off.alpha.data(), m_on.alpha.data())

def exercise_d_model_target_with_hybrid():
  import mmtbx.refinement.llgi_e_dmodel_target as target
  import mmtbx.refinement.llgi_e_dmodel_fit as fit
  rnd = np.random.default_rng(4)
  n = 120
  s2 = rnd.uniform(0.002, 0.08, n)
  e_eff = rnd.uniform(0.2, 2.5, n)
  e_c = rnd.uniform(0.1, 2.5, n)
  dobs = rnd.uniform(0.4, 0.95, n)
  centric = rnd.random(n) < 0.2
  e2 = e_eff**2 + rnd.normal(0, 0.5, n)
  sig = rnd.uniform(0.2, 2.0, n)
  # half the reflections forced exact, the rest by the rice_kappa rule
  h = ext.llgi_hybrid(e_obs_sq=flex.double(e2), sig_e_obs_sq=flex.double(sig),
    force_exact=flex.bool((np.arange(n) % 2 == 0).tolist()),
    centric_flags=flex.bool(centric.tolist()), rice_kappa=0.1)
  b_k_grid = fit.default_b_k_grid(2, s2)
  theta = np.array([0.6, 0.4, 0.1, 40.0])
  ll, grad = target.total_ll_and_gradient(
    theta, s2, e_eff, e_c, dobs, centric, b_k_grid, hybrid=h)
  ll0, _ = target.total_ll_and_gradient(
    theta, s2, e_eff, e_c, dobs, centric, b_k_grid)
  assert abs(ll - ll0) > 1e-6  # the exact reflections changed the target
  step = 1e-6
  for i in range(theta.size):
    tp = theta.copy(); tp[i] += step
    tm = theta.copy(); tm[i] -= step
    lp, _ = target.total_ll_and_gradient(
      tp, s2, e_eff, e_c, dobs, centric, b_k_grid, hybrid=h)
    lm, _ = target.total_ll_and_gradient(
      tm, s2, e_eff, e_c, dobs, centric, b_k_grid, hybrid=h)
    scale = max(1, abs(grad[i]))
    assert abs((lp - lm) / (2 * step) - grad[i]) < 1e-5 * scale, i

def exercise_k1_scale_carries_feff_and_resn():
  """ apply_scale_k1_to_f_obs must rescale llgi_data's FEFF and RESN with
  f_obs, so the F-scale target sees model and data on one scale; the E
  scale quantities do not change. """
  fmodel = build_fmodel(n_atoms=40, d_min=2.2, seed=5)
  llgi_data = synthetic_llgi_data(fmodel, seed=105)
  f_obs = fmodel.f_obs()
  scaled = f_obs.customized_copy(data=f_obs.data() * 4)
  fmodel.update(f_obs=scaled)
  llgi_data = llgi_hybrid.replace_llgi_data(llgi_data,
    feff=llgi_data.feff.customized_copy(data=llgi_data.feff.data() * 4),
    resn=llgi_data.resn.customized_copy(data=llgi_data.resn.data() * 4))
  fmodel.set_llgi_data(llgi_data)
  eeff0 = llgi_data.feff.data() / llgi_data.resn.data()
  ratio0 = flex.mean(llgi_data.feff.data()) / flex.mean(fmodel.f_obs().data())
  fmodel.apply_scale_k1_to_f_obs()
  new = fmodel.llgi_data()
  assert approx_equal(flex.mean(fmodel.f_obs().data()),
    flex.mean(f_obs.data()), eps=0.2 * flex.mean(f_obs.data()))
  assert approx_equal(
    flex.mean(new.feff.data()) / flex.mean(fmodel.f_obs().data()), ratio0)
  assert approx_equal(new.feff.data() / new.resn.data(), eeff0)

def exercise_intensities_from_amplitudes():
  """ cctbx French-Wilson amplitudes from synthetic intensities, then the
  reconstruction: close to the original intensities for most reflections;
  plain amplitudes are recognised as not French-Wilson. """
  from cctbx import french_wilson
  import io
  fmodel = build_fmodel(n_atoms=60, d_min=1.8, seed=7)
  f_calc = abs(fmodel.f_model())
  rnd = random.Random(3)
  i_true = flex.pow2(f_calc.data())
  mean_i = flex.mean(i_true)
  sig = flex.double([0.05 * x + 0.02 * mean_i * (1 + rnd.random())
    for x in i_true])
  i_obs = f_calc.customized_copy(
    data=i_true + sig * flex.double([rnd.gauss(0, 1) for x in i_true]),
    sigmas=sig).set_observation_type_xray_intensity()
  f_fw = french_wilson.french_wilson_scale(miller_array=i_obs,
    log=io.StringIO())
  log = io.StringIO()
  r = llgi_hybrid.intensities_from_amplitudes(f_obs=f_fw, log=log)
  assert r.french_wilson
  assert "WARNING: the LLGI data were derived from amplitudes" in log.getvalue()
  i_c, rec = i_obs.common_sets(r.i_obs)
  z = flex.abs(rec.data() - i_c.data()) / i_c.sigmas()
  zs = sorted(z)
  assert zs[len(zs)//2] < 0.05, zs[len(zs)//2]
  assert zs[int(0.95*len(zs))] < 0.5, zs[int(0.95*len(zs))]
  ratio = rec.sigmas() / i_c.sigmas()
  assert abs(sorted(ratio)[len(ratio)//2] - 1) < 0.01
  # plain amplitudes (not French-Wilson): I = F^2, SIGI = 2 F SIGF
  f_plain = f_calc.customized_copy(data=f_calc.data(),
    sigmas=f_calc.data() * 0.9)  # SIGF/F above the French-Wilson limit
  r2 = llgi_hybrid.intensities_from_amplitudes(f_obs=f_plain)
  assert not r2.french_wilson
  assert approx_equal(r2.i_obs.data(), flex.pow2(f_calc.data()))

def run():
  exercise_helpers_and_select()
  exercise_map_coefficients_exact_branch()
  exercise_map_coefficients_exact_without_hybrid()
  exercise_d_model_target_with_hybrid()
  exercise_k1_scale_carries_feff_and_resn()
  exercise_intensities_from_amplitudes()
  print("OK")

if (__name__ == "__main__"):
  run()

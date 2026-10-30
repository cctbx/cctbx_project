from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_maps as llgi_maps
import mmtbx.refinement.llgi_phi_table as llgi_phi_table
import math

# Reference values from the map-coefficient handoff, computed independently
# by numerical quadrature of the posterior:
# (E_obs^2, sigma, sigmaA, Ec, centric, phi, <E>, M)
handoff = [
  (5.0, 2.0, 0.6, 1.5, False, 0.047302, 1.32072349, 0.70513705),
  (5.0, 2.0, 0.9, 2.0, False, 0.033725, 1.89729238, 0.72703561),
  (46.5, 30.0, 0.6, 1.5, False, 0.000272, 0.92928558, 0.04613415),
  (46.5, 30.0, 0.9, 2.0, False, 0.000195, 1.81641909, 0.08819034),
  (5.0, 2.0, 0.6, 1.5, True, 0.147178, 1.57175276, 1.16050712),
  (5.0, 2.0, 0.3, 1.0, True, 0.062550, 0.72673069, 0.48654580),
  (1.0, 1.0, 0.6, 1.5, True, 0.261818, 0.66594126, 0.00214839),
  (1.0, 1.0, 0.3, 1.0, True, 0.102838, 0.20391978, -0.07095832)]

def exercise_handoff_fixtures():
  col = lambda k: [h[k] for h in handoff]
  e = llgi_maps.posterior_mean_e(
    e_obs_sq=flex.double(col(0)), sig_e_obs_sq=flex.double(col(1)),
    e_calc=flex.double(col(3)), sigmaa=flex.double(col(2)),
    centric_flags=flex.bool(col(4)))
  assert approx_equal(e, col(6), eps=1e-7)
  # bias-reduced coefficient with the handoff's own phi
  c = llgi_maps.e_scale_coefficients(e_expected=flex.double(col(6)),
    e_calc=flex.double(col(3)), sigmaa=flex.double(col(2)),
    phi=flex.double(col(5)), centric_flags=flex.bool(col(4)))
  assert approx_equal(c.bias_reduced, col(7), eps=1e-7)
  assert approx_equal(c.difference,
    [h[6] - h[2]*h[3] for h in handoff], eps=1e-12)
  # fill: 2<E> - sigmaA*Ec (acentric), <E> (centric)
  assert approx_equal(c.fill,
    [h[6] if h[4] else 2*h[6] - h[2]*h[3] for h in handoff], eps=1e-12)
  # tabulated phi (phi_ref's acentric values at large sigma are slightly
  # high; the table follows the large-sigma limit)
  for h in handoff:
    p = llgi_phi_table.phi_one(h[2], h[1], h[4])
    if(h[1] < 10): assert approx_equal(p, h[5], eps=3e-4), h
    else: assert abs(p/(h[2]**2*(1-h[2]**2)/h[1]**2) - 1) < 0.02, h

def exercise_phi_table_limits():
  # Monte Carlo values (phi_quad.py), sigma = 0, and sigmaA -> 0
  for sa, c, ref in [(0.3, False, 0.074), (0.9, False, 0.405),
                     (0.99, False, 0.485), (0.6, True, 0.389),
                     (0.99, True, 0.911)]:
    assert abs(llgi_phi_table.phi_one(sa, 0, c) - ref) < 0.003
  assert llgi_phi_table.phi_one(0, 1.0, False) == 0
  # beyond the table: 1/sigma^2 with the large-sigma constant
  for c, k in [(False, 1), (True, 4)]:
    p = llgi_phi_table.phi_one(0.7, 100.0, c)
    assert abs(p/(k*0.49*0.51/100.0**2) - 1) < 0.02
  # monotone decreasing in sigma
  ps = [llgi_phi_table.phi_one(0.8, s, False) for s in [0, 0.3, 1, 3, 10, 30]]
  assert all(a > b for a, b in zip(ps, ps[1:]))

def exercise_no_measurement_error():
  # sigma = 0: |E| known, <E> = |E| I1/I0(X) (acentric), |E| tanh(X/2)
  e = llgi_maps.posterior_mean_e(
    e_obs_sq=flex.double([4.0, 4.0, -0.1]), sig_e_obs_sq=flex.double(3, 0),
    e_calc=flex.double([1.5, 1.5, 1.0]), sigmaa=flex.double(3, 0.6),
    centric_flags=flex.bool([False, True, False]))
  x = 2*0.6*1.5*2/(1 - 0.36)
  import scitbx.math
  assert approx_equal(e[0], 2*scitbx.math.bessel_i1_over_i0(x), eps=1e-9)
  assert approx_equal(e[1], 2*math.tanh(x/2), eps=1e-12)
  assert e[2] == 0

def exercise_fit_anisobeta():
  # an exact quadratic form is recovered; evaluation at other indices
  idx = flex.miller_index([(h, k, l) for h in range(-4, 5) for k in range(0, 5)
    for l in range(1, 6)])
  b = [-0.010, -0.004, 0.006, 0.002, -0.003, 0.001]
  f = lambda h: math.exp(b[0]*h[0]**2 + b[1]*h[1]**2 + b[2]*h[2]**2 +
    b[3]*h[0]*h[1] + b[4]*h[0]*h[2] + b[5]*h[1]*h[2])
  fit = llgi_maps.fit_anisobeta(idx, flex.double([f(h) for h in idx]))
  assert approx_equal(list(fit.coefficients), b, eps=1e-10)
  assert fit.rms_log_residual < 1e-10
  other = flex.miller_index([(7, 3, -2), (0, 9, 4)])
  assert approx_equal(fit.evaluate(other), [f(h) for h in other], eps=1e-9)

def exercise_fit_anisobeta_resolution_dependent():
  # nacelle's resolution-dependent anisotropy, exp(h^T B0 h + d*^2 h^T B1 h),
  # is recovered exactly with the unit cell (twelve coefficients); a plain
  # quadratic form gives B1 = 0; without the unit cell the fit is poorer
  from cctbx import uctbx
  uc = uctbx.unit_cell((51, 62, 73, 90, 97, 90))
  idx = flex.miller_index([(h, k, l) for h in range(-6, 7) for k in range(0, 7)
    for l in range(1, 8)])
  b0 = [-0.010, -0.004, 0.006, 0.002, -0.003, 0.001]
  b1 = [0.20, -0.10, 0.05, 0.0, 0.08, -0.03]
  m = lambda h: [h[0]**2, h[1]**2, h[2]**2, h[0]*h[1], h[0]*h[2], h[1]*h[2]]
  f = lambda h, c1: math.exp(sum(c*x for c, x in zip(b0, m(h))) +
    uc.d_star_sq(h) * sum(c*x for c, x in zip(c1, m(h))))
  data = flex.double([f(h, b1) for h in idx])
  fit = llgi_maps.fit_anisobeta(idx, data, unit_cell=uc)
  assert approx_equal(list(fit.coefficients), b0 + b1, eps=1e-8)
  assert fit.rms_log_residual < 1e-10
  other = flex.miller_index([(7, 3, -2), (0, 9, 4)])
  assert approx_equal(fit.evaluate(other), [f(h, b1) for h in other], eps=1e-8)
  fit = llgi_maps.fit_anisobeta(idx, flex.double([f(h, [0]*6) for h in idx]),
    unit_cell=uc)
  assert approx_equal(list(fit.coefficients), b0 + [0]*6, eps=1e-8)
  assert llgi_maps.fit_anisobeta(idx, data).rms_log_residual > 1e-3

def exercise_compute_on_fmodel():
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  fmodel = build_llgi_fmodel(40, 2.2, seed=3, rice_kappa=0.1)
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  # anisotropy as nacelle writes it, exp(h^T B h)
  quad = lambda h: math.exp(-0.002*h[0]**2 + 0.001*h[1]**2 + 0.0005*h[2]**2
    + 0.0003*h[0]*h[1])
  aniso = f_obs.array(data=flex.double([quad(h) for h in f_obs.indices()]))
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(llgi_data,
    anisobeta=aniso))
  # default: the maps with nacelle's anisotropy removed
  r = llgi_maps.compute(fmodel)
  assert [a[0] for a in r.arrays] == ["2FOFCWT", "2FOFCWT_fill", "FOFCWT"]
  # with anisotropic_maps: the _aniso versions as well
  p = llgi_maps.llgi_map_params.extract()
  p.anisotropic_maps = True
  r = llgi_maps.compute(fmodel, params=p)
  labels = [a[0] for a in r.arrays]
  assert labels == ["2FOFCWT", "2FOFCWT_fill", "FOFCWT", "2FOFCWT_aniso", "2FOFCWT_fill_aniso",
    "FOFCWT_aniso"], labels
  arrays = dict(r.arrays)
  n = f_obs.size()
  for k in ["2FOFCWT", "FOFCWT", "2FOFCWT_aniso", "FOFCWT_aniso"]:
    assert arrays[k].size() == n
    assert arrays[k].indices().all_eq(f_obs.indices())
  for k in ["2FOFCWT_fill", "2FOFCWT_fill_aniso"]:
    assert arrays[k].size() == n + r.n_filled
  ce = r.coefficients_e
  resn = fmodel.llgi_data().resn.data()   # TEPS = 1 here
  # default F scale RESN/sqrt(ANISOBETA); _aniso: RESN
  assert approx_equal(r.f_scale, resn / flex.sqrt(aniso.data()), eps=1e-12)
  assert approx_equal(r.f_scale_aniso, resn, eps=1e-12)
  # signed coefficients along the model phase
  phase = flex.polar(flex.double(n, 1),
    flex.arg(fmodel.f_model_no_aniso_scale().data()))
  def along_phase(v): return flex.complex_double(v, flex.double(n, 0)) * phase
  for lab, coef, scale in [("2FOFCWT", ce.bias_reduced, r.f_scale),
      ("FOFCWT", ce.difference, r.f_scale),
      ("2FOFCWT_aniso", ce.bias_reduced, resn),
      ("FOFCWT_aniso", ce.difference, resn)]:
    assert approx_equal(arrays[lab].data(), along_phase(coef * scale),
      eps=1e-9), lab
  # the fill maps on the reflections used: the fill coefficient
  for lab, scale in [("2FOFCWT_fill", r.f_scale), ("2FOFCWT_fill_aniso", resn)]:
    obs = arrays[lab].common_set(f_obs)
    assert approx_equal(obs.data(), along_phase(ce.fill * scale),
      eps=1e-9), lab
  # added reflections (sigmaA*Ec): the two fill maps differ by
  # sqrt(ANISOBETA(h))
  assert r.n_filled > 0
  a = arrays["2FOFCWT_fill_aniso"].lone_set(f_obs)
  b = arrays["2FOFCWT_fill"].common_set(a)
  ratio = flex.abs(a.data()) / flex.abs(b.data())
  assert approx_equal(ratio,
    flex.sqrt(flex.double([quad(h) for h in a.indices()])), eps=1e-6)
  # the fill (map_tools.model_missing_reflections_llgi: Ec from the whole
  # model on the complete set, the fitted sigmaA curve, resn_model's RESN)
  # evaluated for the observed reflections: sigmaA*Ec on their F scale
  from mmtbx.map_tools import model_missing_reflections_llgi
  mro = model_missing_reflections_llgi(fmodel=fmodel)
  mro.get_missing(isotropize=False)
  mis = mro.e_scale_missing
  rm = llgi_maps.resn_model(fmodel.llgi_data())
  sa = mis.sigmaa.deep_copy()
  sa.set_selected(sa > llgi_maps.SIGMAA_MAX, llgi_maps.SIGMAA_MAX)
  fill = mis.miller_set.array(data=sa * mis.e_model_abs * mis.resn /
    flex.sqrt(rm.anisobeta(mis.miller_set))).common_set(f_obs)
  assert fill.indices().all_eq(f_obs.indices())
  expected = r.sigmaa * r.e_calc * r.f_scale
  assert flex.max(flex.abs(fill.data() - expected) / expected) < 0.01
  assert flex.max(r.sigmaa) <= llgi_maps.SIGMAA_MAX
  # without ANISOBETA: the maps keep the anisotropy (scale RESN) and no
  # _aniso maps are written
  fmodel.set_llgi_data(llgi_data)
  r = llgi_maps.compute(fmodel, params=p)
  assert [a[0] for a in r.arrays] == ["2FOFCWT", "2FOFCWT_fill", "FOFCWT"]
  assert approx_equal(r.f_scale, resn, eps=1e-12)

def exercise_resn_model():
  # RESN as nacelle computes it, epsilon*G(s^2)*ANISOBETA(h) (G smooth in
  # resolution, with a step between bins; ANISOBETA a quadratic form):
  # reflections held out are recovered from the rest; nacelle's own values
  # are returned as they are; nacelle_all takes precedence over the refined
  # subset; without ANISOBETA the model is isotropic
  from cctbx import crystal
  from libtbx import group_args
  cs = crystal.symmetry(unit_cell=(51, 62, 73, 90, 97, 90), space_group_symbol="C 2")
  ms = cs.build_miller_set(anomalous_flag=False, d_min=2.5)
  dss = ms.d_star_sq().data()
  g = flex.exp(-12 * dss) * (1 + 0.3 * flex.cos(5 * dss))
  g.set_selected(dss > 0.08, g.select(dss > 0.08) * 1.05)    # bin step
  quad = lambda h: math.exp(0.004*h[0]**2 - 0.002*h[1]**2 + 0.001*h[2]**2
    - 0.0007*h[0]*h[2])
  ab = flex.double([quad(h) for h in ms.indices()])
  resn = ms.array(data=flex.sqrt(ms.epsilons().data().as_double() * g * ab))
  anisobeta = ms.array(data=ab)
  # the weakest directions held out, as phenix's sigma(E^2) cutoff does
  held = ab < 0.6
  assert 0 < held.count(True) < ms.size() // 2
  m = llgi_maps.resn_model(group_args(resn=resn.select(~held),
    anisobeta=anisobeta.select(~held)))
  assert m.fit.rms_log_residual < 1e-8
  pred, exact = m.resn_for(resn.select(held))
  assert exact.count(True) == 0
  rel = flex.abs(pred / resn.data().select(held) - 1)
  assert flex.max(rel) < 0.06 and flex.mean(rel) < 2e-3, (flex.max(rel), flex.mean(rel))
  # nacelle's values where it has them, in the order asked
  sub = resn.select(~held)
  perm = flex.random_permutation(sub.size())
  pred, exact = m.resn_for(sub.select(perm))
  assert exact.all_eq(True)
  assert approx_equal(pred, sub.data().select(perm), eps=0)
  # nacelle_all: every reflection, though only some are refined
  m = llgi_maps.resn_model(group_args(resn=sub, anisobeta=anisobeta.select(~held),
    nacelle_all=group_args(resn=resn, anisobeta=anisobeta)))
  pred, exact = m.resn_for(resn.select(held))
  assert exact.all_eq(True)
  # put on the scale of the refined RESN (rescaled with f_obs by
  # apply_scale_k1_to_f_obs)
  m = llgi_maps.resn_model(group_args(resn=sub.customized_copy(
    data=sub.data() * 0.1), anisobeta=anisobeta.select(~held),
    nacelle_all=group_args(resn=resn, anisobeta=anisobeta)))
  pred, exact = m.resn_for(resn.select(held))
  assert approx_equal(pred, resn.data().select(held) * 0.1, eps=1e-9)
  # resolution-dependent anisotropy (phasertng.macaniso.resolution_dependent):
  # held-out reflections recovered as well
  quartic = lambda h: quad(h) * math.exp(cs.unit_cell().d_star_sq(h) * (
    0.05*h[0]**2 - 0.03*h[1]**2 + 0.02*h[2]**2 + 0.01*h[0]*h[2]))
  ab4 = flex.double([quartic(h) for h in ms.indices()])
  resn4 = ms.array(data=flex.sqrt(ms.epsilons().data().as_double() * g * ab4))
  held4 = ab4 < 0.6
  m = llgi_maps.resn_model(group_args(resn=resn4.select(~held4),
    anisobeta=ms.array(data=ab4).select(~held4)))
  assert m.fit.rms_log_residual < 1e-8
  pred, exact = m.resn_for(resn4.select(held4))
  rel = flex.abs(pred / resn4.data().select(held4) - 1)
  assert flex.max(rel) < 0.06 and flex.mean(rel) < 2e-3, (flex.max(rel), flex.mean(rel))
  # no ANISOBETA: isotropic, G = RESN^2/epsilon
  iso = ms.array(data=flex.sqrt(ms.epsilons().data().as_double() * g))
  m = llgi_maps.resn_model(group_args(resn=iso.select(~held), anisobeta=None))
  assert not m.anisotropic
  pred, exact = m.resn_for(iso.select(held))
  assert flex.max(flex.abs(pred / iso.data().select(held) - 1)) < 0.06

def exercise_fit_k_anisotropic():
  # both forms of fmodel's anisotropic scale (mmtbx.bulk_solvent.scaler) are
  # recovered exactly and evaluated for other reflections
  from cctbx import crystal, adptbx
  from mmtbx.f_model import ext as f_model_ext
  cs = crystal.symmetry(unit_cell=(51, 62, 73, 90, 97, 90), space_group_symbol="C 2")
  ms = cs.build_miller_set(anomalous_flag=False, d_min=2.5)
  sel = flex.random_bool(ms.size(), 0.8)
  fit_on, other = ms.select(sel), ms.select(~sel)
  u_star = adptbx.u_cart_as_u_star(cs.unit_cell(), adptbx.b_as_u((4, -3, 2, 0, 1.5, 0)))
  a = flex.double([0.002, -0.3, -0.001, 0.2, 0.0005, 0.1, 0, 0, 0.0003, -0.05, 0, 0])
  for form, k in [("exp", lambda m: f_model_ext.k_anisotropic(m.indices(), u_star)),
                  ("polynomial", lambda m: f_model_ext.k_anisotropic(m.indices(), a, cs.unit_cell()))]:
    f = llgi_maps.fit_k_anisotropic(fit_on.array(data=k(fit_on)))
    assert f.form == form, (f.form, form)
    assert f.rms_log_residual < 1e-8
    assert approx_equal(f.evaluate(other), k(other), eps=1e-6)

def exercise_model_building_maps():
  # fmodel.electron_density_map_for_model_building with real_space_map =
  # bias_reduced: the 2mFo-DFc-type maps (whatever their spelling) are the
  # bias-reduced map,
  # from the working reflections when asked (exclude_free_r_reflections,
  # use_all_data=False); other maps come from the LLGI helper; with
  # posterior_2mfo_dfc, or another target, the helper itself
  from cctbx import miller
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  from mmtbx.regression.llgi_test_utils import build_llgi_fmodel
  fmodel = build_llgi_fmodel(30, 2.2, seed=3, rice_kappa=0.1)
  ld = fmodel.llgi_data()
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld,
    real_space_map="bias_reduced"))
  emap = fmodel.electron_density_map_for_model_building()
  assert isinstance(emap, llgi_maps.bias_reduced_map_server)
  iwt = dict(llgi_maps.compute(fmodel).arrays)[llgi_maps.LABELS.bias_reduced]
  work = fmodel.arrays.work_sel
  assert 0 < work.count(False) < work.size()
  for map_type in ["2mFo-DFc", "2mFobs-DFmodel", "2mFo-DFc_filled"]:
    mc = emap.map_coefficients(map_type=map_type, fill_missing=True,
      isotropize=True)
    assert mc.indices().all_eq(iwt.indices())
    assert approx_equal(mc.data(), iwt.data())
  mc = emap.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True)
  assert mc.indices().all_eq(iwt.indices().select(work))
  assert approx_equal(mc.data(), iwt.data().select(work))
  helper = fmodel.electron_density_map_for_target()
  for map_type in ["mFo-DFc", "mFobs-DFmodel", "Fmodel", "mFo"]:
    assert approx_equal(emap.map_coefficients(map_type=map_type).data(),
      helper.map_coefficients(map_type=map_type).data())
  # the map from the working set, on the grid of another map
  m_all = emap.fft_map(map_type="2mFo-DFc", resolution_factor=0.25)
  m_work = emap.fft_map(map_type="2mFo-DFc", resolution_factor=0.25,
    use_all_data=False, other_fft_map=m_all)
  ref = miller.fft_map(crystal_gridding=m_all,
    fourier_coefficients=iwt.select(work))
  assert approx_equal(m_work.real_map_unpadded(), ref.real_map_unpadded())
  assert flex.max(flex.abs(m_work.real_map_unpadded() -
    m_all.real_map_unpadded())) > 0
  md = emap.fft_map(map_type="mFo-DFc", resolution_factor=0.25,
    use_all_data=False, other_fft_map=m_all)
  ref = miller.fft_map(crystal_gridding=m_all, fourier_coefficients=
    helper.map_coefficients(map_type="mFo-DFc").select(work))
  assert approx_equal(md.real_map_unpadded(), ref.real_map_unpadded())
  # posterior_filled: for the real-space steps (fft_map, working
  # reflections), the posterior map of the working reflections with D*Fc
  # (the LLGI fill) for the free reflections and those without data, as
  # the helper's map_coefficients gives with exclude_free_r_reflections and
  # fill_missing; without fill_missing (water picking), the working
  # reflections alone
  from mmtbx.map_tools import model_missing_reflections_llgi
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld,
    real_space_map="posterior_filled"))
  emap = fmodel.electron_density_map_for_model_building()
  assert isinstance(emap, llgi_maps.filled_posterior_map_server)
  fill = model_missing_reflections_llgi(fmodel=fmodel).get_missing()
  post = helper.map_coefficients(map_type="2mFo-DFc").select(work)
  assert approx_equal(emap.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True).data(), post.data())
  for map_type in ["2mFo-DFc", "2mFobs-DFmodel"]:
    mc = emap.map_coefficients(map_type=map_type,
      exclude_free_r_reflections=True, fill_missing=True)
    assert mc.size() == fill.size()
    a, b = mc.common_sets(post)
    assert a.size() == post.size()
    assert approx_equal(a.data(), b.data())
    a, b = mc.lone_set(post).common_sets(fill)
    assert a.size() == fill.size() - post.size()
    assert approx_equal(a.data(), b.data())
  m_work = emap.fft_map(map_type="2mFo-DFc", resolution_factor=0.25,
    use_all_data=False, other_fft_map=m_all)
  ref = miller.fft_map(crystal_gridding=m_all,
    fourier_coefficients=emap.map_coefficients(map_type="2mFo-DFc",
      exclude_free_r_reflections=True, fill_missing=True))
  assert approx_equal(m_work.real_map_unpadded(), ref.real_map_unpadded())
  for kw in [dict(), dict(fill_missing=False)]:
    assert approx_equal(emap.map_coefficients(map_type="2mFo-DFc", **kw).data(),
      helper.map_coefficients(map_type="2mFo-DFc", **kw).data())
  md = emap.fft_map(map_type="mFo-DFc", resolution_factor=0.25,
    use_all_data=False, other_fft_map=m_all)
  ref = miller.fft_map(crystal_gridding=m_all, fourier_coefficients=
    helper.map_coefficients(map_type="mFo-DFc").select(work))
  assert approx_equal(md.real_map_unpadded(), ref.real_map_unpadded())
  # posterior_2mfo_dfc: the LLGI helper; target ml: the ordinary map
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld,
    real_space_map="posterior_2mfo_dfc"))
  emap = fmodel.electron_density_map_for_model_building()
  assert not isinstance(emap, llgi_maps.bias_reduced_map_server)
  assert approx_equal(emap.map_coefficients(map_type="2mFo-DFc").data(),
    helper.map_coefficients(map_type="2mFo-DFc").data())
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(ld,
    real_space_map="bias_reduced"))
  fmodel.set_target_name("ml")
  emap = fmodel.electron_density_map_for_model_building()
  assert not isinstance(emap, llgi_maps.bias_reduced_map_server)
  assert approx_equal(emap.map_coefficients(map_type="2mFo-DFc").data(),
    fmodel.electron_density_map().map_coefficients(
      map_type="2mFo-DFc").data())

def run():
  exercise_handoff_fixtures()
  exercise_phi_table_limits()
  exercise_no_measurement_error()
  exercise_fit_anisobeta()
  exercise_fit_anisobeta_resolution_dependent()
  exercise_compute_on_fmodel()
  exercise_resn_model()
  exercise_fit_k_anisotropic()
  exercise_model_building_maps()
  print("OK")

if (__name__ == "__main__"):
  run()

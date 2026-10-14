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
  assert [a[0] for a in r.arrays] == ["IWT", "IWT_fill", "dIWT"]
  # with anisotropic_maps: the _aniso versions as well
  p = llgi_maps.llgi_map_params.extract()
  p.anisotropic_maps = True
  r = llgi_maps.compute(fmodel, params=p)
  labels = [a[0] for a in r.arrays]
  assert labels == ["IWT", "IWT_fill", "dIWT", "IWT_aniso", "IWT_fill_aniso",
    "dIWT_aniso"], labels
  arrays = dict(r.arrays)
  n = f_obs.size()
  for k in ["IWT", "dIWT", "IWT_aniso", "dIWT_aniso"]:
    assert arrays[k].size() == n
    assert arrays[k].indices().all_eq(f_obs.indices())
  for k in ["IWT_fill", "IWT_fill_aniso"]:
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
  for lab, coef, scale in [("IWT", ce.bias_reduced, r.f_scale),
      ("dIWT", ce.difference, r.f_scale),
      ("IWT_aniso", ce.bias_reduced, resn),
      ("dIWT_aniso", ce.difference, resn)]:
    assert approx_equal(arrays[lab].data(), along_phase(coef * scale),
      eps=1e-9), lab
  # the fill maps on the reflections used: the fill coefficient
  for lab, scale in [("IWT_fill", r.f_scale), ("IWT_fill_aniso", resn)]:
    obs = arrays[lab].common_set(f_obs)
    assert approx_equal(obs.data(), along_phase(ce.fill * scale),
      eps=1e-9), lab
  # added reflections (sigmaA*Ec): on the level of the observed ones at
  # similar resolution, and the two fill maps differ by sqrt(ANISOBETA(h))
  assert r.n_filled > 0
  a = arrays["IWT_fill_aniso"].lone_set(f_obs)
  b = arrays["IWT_fill"].common_set(a)
  ratio = flex.abs(a.data()) / flex.abs(b.data())
  assert approx_equal(ratio,
    flex.sqrt(flex.double([quad(h) for h in a.indices()])), eps=1e-6)
  obs = arrays["IWT_fill"].common_set(f_obs)
  ss_obs = obs.d_star_sq().data()
  for amp, ss in zip(flex.abs(b.data()), b.d_star_sq().data()):
    order = flex.sort_permutation(flex.abs(ss_obs - ss))[:50]
    near = flex.mean(flex.abs(obs.data()).select(order))
    assert 0.05 * near < amp < 20 * near, (amp, near)
  assert flex.max(r.sigmaa) <= llgi_maps.SIGMAA_MAX
  # without ANISOBETA: the maps keep the anisotropy (scale RESN) and no
  # _aniso maps are written
  fmodel.set_llgi_data(llgi_data)
  r = llgi_maps.compute(fmodel, params=p)
  assert [a[0] for a in r.arrays] == ["IWT", "IWT_fill", "dIWT"]
  assert approx_equal(r.f_scale, resn, eps=1e-12)

def run():
  exercise_handoff_fixtures()
  exercise_phi_table_limits()
  exercise_no_measurement_error()
  exercise_fit_anisobeta()
  exercise_compute_on_fmodel()
  print("OK")

if (__name__ == "__main__"):
  run()

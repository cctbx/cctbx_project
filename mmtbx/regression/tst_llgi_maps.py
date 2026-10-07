from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import mmtbx.refinement.llgi_maps as llgi_maps
import mmtbx.refinement.llgi_phi_table as llgi_phi_table
import math

# Map-coefficient handoff fixtures (fixtures.py / mfix.py, Oct 2026):
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
    phi=flex.double(col(5)))
  assert approx_equal(c.bias_reduced, col(7), eps=1e-7)
  assert approx_equal(c.difference,
    [h[6] - h[2]*h[3] for h in handoff], eps=1e-12)
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

def exercise_compute_on_fmodel():
  from mmtbx.regression.tst_llgi_hybrid import build_hybrid_fmodel
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  fmodel = build_hybrid_fmodel(rice_kappa=0.1)
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  aniso = f_obs.array(data=flex.double(
    [0.8 + 0.4*((i*37) % 11)/10. for i in range(f_obs.size())]))
  fmodel.set_llgi_data(llgi_hybrid.replace_llgi_data(llgi_data,
    anisobeta=aniso))
  r = llgi_maps.compute(fmodel)
  labels = [a[0] for a in r.arrays]
  assert labels == ["LLGI_FILLED", "LLGI_BIASRED", "LLGI_DIFF",
    "LLGI_BIASRED_E", "LLGI_DIFF_E"], labels
  arrays = dict(r.arrays)
  n = f_obs.size()
  assert arrays["LLGI_FILLED"].size() == n + r.n_filled
  for k in labels[1:]:
    assert arrays[k].size() == n
    assert arrays[k].indices().all_eq(f_obs.indices())
  # F scale = RESN/sqrt(ANISOBETA)
  resn = fmodel.llgi_data().resn.data()
  ratio = flex.abs(arrays["LLGI_DIFF"].data()) / flex.abs(
    arrays["LLGI_DIFF_E"].data()).set_selected(
      flex.abs(arrays["LLGI_DIFF_E"].data()) == 0, 1)
  sel = flex.abs(arrays["LLGI_DIFF_E"].data()) > 1e-6
  assert approx_equal(ratio.select(sel),
    (resn / flex.sqrt(aniso.data())).select(sel), eps=1e-9)
  # sigmaA clamp, and the phase of the coefficients is the model phase
  assert flex.max(r.sigmaa) <= llgi_maps.SIGMAA_MAX
  # without anisotropy removal
  p = llgi_maps.llgi_map_params.extract()
  p.remove_anisotropy = False
  p.e_scale = False
  r2 = llgi_maps.compute(fmodel, params=p)
  assert [a[0] for a in r2.arrays] == labels[:3]
  d2 = dict(r2.arrays)["LLGI_DIFF"]
  assert approx_equal(flex.abs(d2.data()).select(sel),
    (flex.abs(arrays["LLGI_DIFF_E"].data()) * resn).select(sel), eps=1e-9)

def run():
  exercise_handoff_fixtures()
  exercise_phi_table_limits()
  exercise_no_measurement_error()
  exercise_compute_on_fmodel()
  print("OK")

if (__name__ == "__main__"):
  run()

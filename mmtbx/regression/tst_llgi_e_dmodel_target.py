from __future__ import absolute_import, division, print_function
import numpy as np
import mmtbx.refinement.llgi_e_dmodel_target as target

def _build_reflections(rnd, n=30):
  s2 = rnd.uniform(0.0, 1.2, size=n)
  e_eff = rnd.uniform(0.3, 2.2, size=n)
  e_c = rnd.uniform(0.3, 2.2, size=n)
  dobs = rnd.uniform(0.3, 0.95, size=n)
  centric_flags = rnd.uniform(size=n) < 0.25
  return s2, e_eff, e_c, dobs, centric_flags

def _random_theta_and_grid(rnd, k):
  a = rnd.uniform(0.05, 0.6, size=k)
  b_k_grid = np.sort(rnd.uniform(1.0, 100.0, size=k))
  b = rnd.uniform(0.02, 0.3)
  b_defect = rnd.uniform(5.0, 150.0)
  theta = np.empty(k + 2, dtype=float)
  theta[:k] = a
  theta[-2] = b
  theta[-1] = b_defect
  return theta, b_k_grid

def _ll_only(theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid):
  ll, _ = target.total_ll_and_gradient(
    theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  return ll

def exercise_gradient_matches_finite_difference_acentric_only():
  rnd = np.random.RandomState(10)
  s2, e_eff, e_c, dobs, centric_flags = _build_reflections(rnd)
  centric_flags = np.zeros_like(centric_flags, dtype=bool)  # force all acentric
  theta, b_k_grid = _random_theta_and_grid(rnd, 2)
  ll, grad = target.total_ll_and_gradient(
    theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  h = 1.e-6
  worst = 0.0
  for i in range(theta.size):
    tp = theta.copy(); tp[i] += h
    tm = theta.copy(); tm[i] -= h
    fd = (_ll_only(tp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
          - _ll_only(tm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)) / (2*h)
    worst = max(worst, abs(grad[i] - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-4, worst

def exercise_gradient_matches_finite_difference_mixed_centric():
  rnd = np.random.RandomState(11)
  s2, e_eff, e_c, dobs, centric_flags = _build_reflections(rnd)
  theta, b_k_grid = _random_theta_and_grid(rnd, 2)
  ll, grad = target.total_ll_and_gradient(
    theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  h = 1.e-6
  worst = 0.0
  for i in range(theta.size):
    tp = theta.copy(); tp[i] += h
    tm = theta.copy(); tm[i] -= h
    fd = (_ll_only(tp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
          - _ll_only(tm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)) / (2*h)
    worst = max(worst, abs(grad[i] - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-4, worst

def exercise_single_reflection_matches_llgi_e_likelihood_directly():
  # A one-reflection "sum" must reduce exactly to
  # llgi_e_likelihood's own l evaluated at D=dobs*D_model(s2).
  import mmtbx.refinement.llgi_e_likelihood as lik
  rnd = np.random.RandomState(13)
  theta, b_k_grid = _random_theta_and_grid(rnd, 1)
  s2 = np.array([0.3])
  e_eff = np.array([1.1])
  e_c = np.array([0.9])
  dobs = np.array([0.6])
  for centric in [False, True]:
    centric_flags = np.array([centric])
    ll, grad = target.total_ll_and_gradient(
      theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
    import mmtbx.refinement.llgi_e_dmodel as dmodel
    D = float(dobs[0] * dmodel.d_model(s2, theta, b_k_grid)[0])
    l_fn = lik.centric_l if centric else lik.acentric_l
    expected_ll = float(l_fn(D, e_eff[0], e_c[0]))
    assert abs(ll - expected_ll) < 1.e-10, (centric, ll, expected_ll)

def run():
  exercise_gradient_matches_finite_difference_acentric_only()
  exercise_gradient_matches_finite_difference_mixed_centric()
  exercise_single_reflection_matches_llgi_e_likelihood_directly()
  print("OK")

if (__name__ == "__main__"):
  run()

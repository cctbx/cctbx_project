from __future__ import absolute_import, division, print_function
import numpy as np
import mmtbx.refinement.llgi_e_dmodel_reparam as rep
import mmtbx.refinement.llgi_e_dmodel_target as target

def _build_case(rnd, n=15, k=2):
  s2 = rnd.uniform(0.0, 1.0, size=n)
  e_eff = rnd.uniform(0.4, 2.0, size=n)
  e_c = rnd.uniform(0.4, 2.0, size=n)
  dobs = rnd.uniform(0.3, 0.9, size=n)
  centric_flags = rnd.uniform(size=n) < 0.3
  a = rnd.uniform(0.1, 0.6, size=k)
  b_k_grid = np.sort(rnd.uniform(5.0, 80.0, size=k))
  b = rnd.uniform(0.05, 0.3)
  b_defect = rnd.uniform(10.0, 100.0)
  theta0 = np.empty(k + 2)
  theta0[:k] = a
  theta0[-2] = b
  theta0[-1] = b_defect
  return s2, e_eff, e_c, dobs, centric_flags, theta0, b_k_grid

def exercise_theta_q_roundtrip():
  rnd = np.random.RandomState(0)
  theta = rnd.uniform(0.01, 100.0, size=8)
  q = rep.q_from_theta(theta)
  theta_back = rep.theta_from_q(q)
  assert np.max(np.abs(theta - theta_back)) < 1.e-10

def exercise_theta_from_q_always_positive():
  q = np.array([-50.0, 0.0, 50.0, -1.e-3, 1.e-3])
  theta = rep.theta_from_q(q)
  assert np.all(theta > 0.0)

def _ll_of_q(q, s2, e_eff, e_c, dobs, centric_flags, b_k_grid):
  theta = rep.theta_from_q(q)
  ll, _, _ = target.total_ll_gradient_hessian(
    theta, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  return ll

def exercise_reparametrized_gradient_matches_finite_difference():
  rnd = np.random.RandomState(1)
  s2, e_eff, e_c, dobs, centric_flags, theta0, b_k_grid = _build_case(rnd)
  q0 = rep.q_from_theta(theta0)
  _, grad_p, _ = target.total_ll_gradient_hessian(
    theta0, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  grad_q = rep.reparametrize_gradient(theta0, grad_p)
  h = 1.e-6
  worst = 0.0
  for i in range(theta0.size):
    qp = q0.copy(); qp[i] += h
    qm = q0.copy(); qm[i] -= h
    fd = (_ll_of_q(qp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
          - _ll_of_q(qm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)) / (2*h)
    worst = max(worst, abs(grad_q[i] - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-5, worst

def exercise_reparametrized_hessian_matches_finite_difference():
  rnd = np.random.RandomState(2)
  s2, e_eff, e_c, dobs, centric_flags, theta0, b_k_grid = _build_case(rnd)
  q0 = rep.q_from_theta(theta0)
  _, grad_p, hess_p = target.total_ll_gradient_hessian(
    theta0, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  hess_q = rep.reparametrize_hessian(theta0, grad_p, hess_p)
  h = 1.e-4
  n = theta0.size
  worst = 0.0
  for i in range(n):
    for j in range(n):
      qpp = q0.copy(); qpp[i] += h; qpp[j] += h
      qpm = q0.copy(); qpm[i] += h; qpm[j] -= h
      qmp = q0.copy(); qmp[i] -= h; qmp[j] += h
      qmm = q0.copy(); qmm[i] -= h; qmm[j] -= h
      fd = (_ll_of_q(qpp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
            - _ll_of_q(qpm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
            - _ll_of_q(qmp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
            + _ll_of_q(qmm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)) / (4*h*h)
      worst = max(worst, abs(hess_q[i, j] - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-3, worst

def exercise_naive_hessian_without_extra_term_is_wrong():
  # Documents/guards the actual bug this module's docstring warns
  # about: dropping the dF/dp*p diagonal term is NOT a negligible
  # simplification away from convergence. The extra term's size is
  # dF/dp_i * p_i -- it is large wherever the ORIGINAL-parameter
  # gradient dF/dp_i is not small, and can be negligible at other
  # (i, j) or other random seeds where dF/dp_i happens to be small at
  # theta0 -- so this checks the DIAGONAL entry with the largest
  # |dF/dp_i * p_i| (guaranteed, by construction, to be where the
  # naive/correct reparametrizations disagree most), rather than a
  # fixed index that could pass or fail depending on the random seed.
  rnd = np.random.RandomState(3)
  s2, e_eff, e_c, dobs, centric_flags, theta0, b_k_grid = _build_case(rnd)
  q0 = rep.q_from_theta(theta0)
  _, grad_p, hess_p = target.total_ll_gradient_hessian(
    theta0, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  hess_q_correct = rep.reparametrize_hessian(theta0, grad_p, hess_p)
  hess_q_naive = hess_p * np.outer(theta0, theta0)
  extra_term = grad_p * theta0
  i = int(np.argmax(np.abs(extra_term)))
  j = i
  h = 1.e-4
  qpp = q0.copy(); qpp[i] += h; qpp[j] += h
  qpm = q0.copy(); qpm[i] += h; qpm[j] -= h
  qmp = q0.copy(); qmp[i] -= h; qmp[j] += h
  qmm = q0.copy(); qmm[i] -= h; qmm[j] -= h
  fd = (_ll_of_q(qpp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
        - _ll_of_q(qpm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
        - _ll_of_q(qmp, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
        + _ll_of_q(qmm, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)) / (4*h*h)
  correct_rel_err = abs(hess_q_correct[i, j] - fd) / max(1.0, abs(fd))
  naive_rel_err = abs(hess_q_naive[i, j] - fd) / max(1.0, abs(fd))
  assert correct_rel_err < 1.e-3, correct_rel_err
  assert naive_rel_err > 5 * correct_rel_err, (naive_rel_err, correct_rel_err)

def exercise_hessian_diagonal_helper_matches_full_hessian_diagonal():
  rnd = np.random.RandomState(4)
  s2, e_eff, e_c, dobs, centric_flags, theta0, b_k_grid = _build_case(rnd)
  _, grad_p, hess_p = target.total_ll_gradient_hessian(
    theta0, s2, e_eff, e_c, dobs, centric_flags, b_k_grid)
  hess_q_full = rep.reparametrize_hessian(theta0, grad_p, hess_p)
  hess_p_diag = np.diag(hess_p)
  diag_helper = rep.reparametrize_hessian_diagonal(
    theta0, grad_p, hess_p_diag)
  assert np.max(np.abs(diag_helper - np.diag(hess_q_full))) < 1.e-12

def run():
  exercise_theta_q_roundtrip()
  exercise_theta_from_q_always_positive()
  exercise_reparametrized_gradient_matches_finite_difference()
  exercise_reparametrized_hessian_matches_finite_difference()
  exercise_naive_hessian_without_extra_term_is_wrong()
  exercise_hessian_diagonal_helper_matches_full_hessian_diagonal()
  print("OK")

if (__name__ == "__main__"):
  run()

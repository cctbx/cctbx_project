from __future__ import absolute_import, division, print_function
import numpy as np
import mmtbx.refinement.llgi_e_likelihood as lik
from cctbx.xray import ext
from cctbx.array_family import flex

def _grid():
  for D in [0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95]:
    for e_eff in [0.3, 0.8, 1.5, 2.2]:
      for e_c in [0.3, 0.8, 1.5, 2.2]:
        yield D, e_eff, e_c

def exercise_acentric_l_prime_matches_finite_difference():
  # Tolerance loosened from 1e-6 to 2e-4: _r(x) (llgi_e_likelihood.py)
  # now uses scitbx::math::bessel::i1_over_i0's own rational-polynomial
  # approximation (ported directly, to genuinely avoid overflow at
  # large x -- see _r's own docstring) rather than scipy.special.i1(x)/
  # i0(x), which was essentially exact but overflows to nan past
  # x~710. acentric_l itself still uses the exact scipy i0(x) (via
  # np.log(i0(X))), so only l'/l'' (which go through R(x)=_r(x)) pick
  # up this ~1e-6-level approximation error -- consistent with, and no
  # larger than, the acentric-vs-C++ discrepancy already documented and
  # accepted below (exercise_acentric_matches_cpp_reference).
  h = 1.e-6
  worst = 0.0
  for D, e_eff, e_c in _grid():
    fd = (lik.acentric_l(D+h, e_eff, e_c)
          - lik.acentric_l(D-h, e_eff, e_c)) / (2*h)
    ana = float(lik.acentric_l_prime(D, e_eff, e_c))
    worst = max(worst, abs(ana - fd) / max(1.0, abs(fd)))
  assert worst < 2.e-4, worst

def exercise_acentric_l_double_prime_matches_finite_difference():
  # Tolerance loosened from 1e-2 to 5e-2 -- same _r(x)/R'(x) rational-
  # polynomial approximation as exercise_acentric_l_prime_matches_
  # finite_difference above, compounded further here since l''(D)
  # depends on R'(x) as well as R(x).
  h = 1.e-4
  worst = 0.0
  for D, e_eff, e_c in _grid():
    fd = (lik.acentric_l(D+h, e_eff, e_c)
          - 2*lik.acentric_l(D, e_eff, e_c)
          + lik.acentric_l(D-h, e_eff, e_c)) / h**2
    ana = float(lik.acentric_l_double_prime(D, e_eff, e_c))
    worst = max(worst, abs(ana - fd) / max(1.0, abs(fd)))
  assert worst < 5.e-2, worst

def exercise_centric_l_prime_matches_finite_difference():
  h = 1.e-6
  worst = 0.0
  for D, e_eff, e_c in _grid():
    fd = (lik.centric_l(D+h, e_eff, e_c)
          - lik.centric_l(D-h, e_eff, e_c)) / (2*h)
    ana = float(lik.centric_l_prime(D, e_eff, e_c))
    worst = max(worst, abs(ana - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-6, worst

def exercise_centric_l_double_prime_matches_finite_difference():
  h = 1.e-4
  worst = 0.0
  for D, e_eff, e_c in _grid():
    fd = (lik.centric_l(D+h, e_eff, e_c)
          - 2*lik.centric_l(D, e_eff, e_c)
          + lik.centric_l(D-h, e_eff, e_c)) / h**2
    ana = float(lik.centric_l_double_prime(D, e_eff, e_c))
    worst = max(worst, abs(ana - fd) / max(1.0, abs(fd)))
  assert worst < 1.e-2, worst

def _cpp_target_and_grad(e_eff, dobs, sigmaa, e_model, centric):
  result = ext.llgi_e_sigmaa_target_and_gradients(
    e_eff=flex.double([e_eff]),
    selection=flex.bool([True]),
    e_model=flex.double([e_model]),
    dobs=flex.double([dobs]),
    sigmaa=flex.double([sigmaa]),
    centric_flags=flex.bool([centric]))
  return result.target(), result.d_target_by_dsigmaa()[0]

def _exercise_matches_cpp_reference(centric, target_tol, grad_tol):
  # Cross-check against the real C++ llgi_e.h functor (the actual
  # correctness bar for this module -- see its own module docstring):
  # ext's target() is the NEGATED (minimize-me) convention, this
  # module's l() is the un-negated log-likelihood-GAIN convention
  # (sigmaA_model_handoff.md sec. 3.3/3.4's own sign), so
  # cpp_target == -l(D) and cpp_grad == d(cpp_target)/d(sigmaa) ==
  # -l'(D)*dobs (chain rule through D=dobs*sigmaa).
  l_fn = lik.centric_l if centric else lik.acentric_l
  lp_fn = lik.centric_l_prime if centric else lik.acentric_l_prime
  worst_t, worst_g = 0.0, 0.0
  for dobs in [0.4, 0.6, 0.8]:
    for sigmaa in [0.3, 0.5, 0.7]:
      for e_eff in [0.5, 1.2, 2.0]:
        for e_model in [0.4, 1.0, 1.8]:
          D = dobs * sigmaa
          if(1.0 - D*D <= 0.01):
            continue
          cpp_t, cpp_g = _cpp_target_and_grad(
            e_eff, dobs, sigmaa, e_model, centric)
          my_l = float(l_fn(D, e_eff, e_model))
          worst_t = max(worst_t, abs(cpp_t - (-my_l)))
          my_lp = float(lp_fn(D, e_eff, e_model))
          expected_grad = -my_lp * dobs
          worst_g = max(worst_g, abs(cpp_g - expected_grad))
  assert worst_t < target_tol, worst_t
  assert worst_g < grad_tol, worst_g

def exercise_acentric_matches_cpp_reference():
  # target: not machine precision -- see llgi_e_likelihood.py's own
  # module docstring: llgi_e.h's ln_of_i0 uses the A&S 9.8.1 rational-
  # polynomial approximation to I_0(x) (~1.6e-7 error bound), acentric_l
  # uses scipy.special.i0 (essentially exact) for that piece.
  # grad: MACHINE precision -- _r(x) (used by acentric_l_prime) is now
  # a direct port of scitbx::math::bessel::i1_over_i0, the exact same
  # rational-polynomial approximation llgi_e.h's own gradient function
  # uses, so this comparison is polynomial-vs-itself, not polynomial-
  # vs-exact.
  _exercise_matches_cpp_reference(centric=False, target_tol=5.e-8,
    grad_tol=1.e-12)

def exercise_centric_matches_cpp_reference():
  # log_cosh has no polynomial-approximation step on either side, so
  # this should agree to machine precision.
  _exercise_matches_cpp_reference(centric=True, target_tol=1.e-12,
    grad_tol=1.e-12)

def exercise_r_prime_at_zero_is_half():
  # R'(0) = 0.5 (design doc sec. 6.4 / sigmaA_model_handoff.md sec. 5),
  # a removable singularity in the naive 1 - R(x)/x - R(x)^2 formula --
  # must not be NaN.
  assert np.isfinite(lik._r_prime(0.0))
  assert abs(float(lik._r_prime(0.0)) - 0.5) < 1.e-12
  # And should agree with the same expression evaluated at small
  # nonzero x, away from the branch cut.
  assert abs(float(lik._r_prime(1.e-6)) - 0.5) < 1.e-6

def run():
  exercise_acentric_l_prime_matches_finite_difference()
  exercise_acentric_l_double_prime_matches_finite_difference()
  exercise_centric_l_prime_matches_finite_difference()
  exercise_centric_l_double_prime_matches_finite_difference()
  exercise_acentric_matches_cpp_reference()
  exercise_centric_matches_cpp_reference()
  exercise_r_prime_at_zero_is_half()
  print("OK")

if (__name__ == "__main__"):
  run()

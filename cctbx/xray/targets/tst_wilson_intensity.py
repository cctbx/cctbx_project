from __future__ import absolute_import, division, print_function
from cctbx.xray import ext as xray_ext
from cctbx.array_family import flex
import math
import random
from scipy import integrate
from scipy.special import erfc

# Oracle: ln of int_0^inf p(J) N(e; J, s^2) dJ, J = I/Sigma.
def ln_quad(e, s, centric):
  def log_integrand(j):
    if centric:
      lp = -j / 2. - 0.5 * math.log(2 * math.pi * j)
    else:
      lp = -j
    return lp - 0.5 * ((e - j) / s) ** 2 - math.log(math.sqrt(2 * math.pi) * s)
  # Break points around the data and the Wilson scale. The integrand is
  # scaled by its maximum on a dense grid so that it does not underflow.
  centre = max(e, 0.)
  top = centre + 40 * s + 60.
  pts = sorted(set([0., min(1e-3, centre + 1.), 1., centre,
    max(centre - 8 * s, 0.), centre + 8 * s, top]))
  grid = [1e-6 + (top - 1e-6) * k / 20000. for k in range(20001)]
  grid += [max(centre + s * k / 10., 1e-6) for k in range(-100, 101)]
  c = max(log_integrand(j) for j in grid)
  total = 0.
  for a, b in zip(pts[:-1], pts[1:]):
    if b > a:
      total += integrate.quad(lambda j: math.exp(log_integrand(j) - c),
        a, b, epsabs=0, epsrel=1e-13, limit=500)[0]
  return c + math.log(total)

grid_e = (-2, -0.3, 0, 0.5, 1, 3, 10, 40)
grid_s = (0.01, 0.1, 1, 5, 30)
grid_sigma = (1., 3.7)

def evaluate(i, si, sigma, centric):
  return xray_ext.wilson_intensity(
    i_obs=flex.double(i), sig_i_obs=flex.double(si),
    sigma=flex.double(sigma), centric_flags=flex.bool(centric))

def check_against_quadrature(centric, tol):
  for sigma in grid_sigma:
    for e in grid_e:
      for s in grid_s:
        r = evaluate([e * sigma], [s * sigma], [sigma], [centric])
        ref = ln_quad(e, s, centric) - math.log(sigma)
        assert abs(r.log_likelihood[0] - ref) < tol, (
          centric, e, s, sigma, r.log_likelihood[0], ref)

def exercise_acentric_matches_quadrature():
  check_against_quadrature(False, 1e-8)

def exercise_centric_matches_quadrature():
  check_against_quadrature(True, 1e-6)

def exercise_derivatives_fd():
  h = 1e-5
  for centric in (False, True):
    for sigma in grid_sigma:
      for e in grid_e:
        for s in grid_s:
          i, si = e * sigma, s * sigma
          def f(ln_sig):
            return evaluate([i], [si], [math.exp(ln_sig)], [centric])
          x0 = math.log(sigma)
          r0 = f(x0)
          rp, rm = f(x0 + h), f(x0 - h)
          fd1 = (rp.log_likelihood[0] - rm.log_likelihood[0]) / (2 * h)
          a1 = r0.d_ll_d_ln_sigma[0]
          assert abs(fd1 - a1) < 1e-6 * max(1., abs(a1)), (centric, e, s, fd1, a1)
          fd2 = (rp.d_ll_d_ln_sigma[0] - rm.d_ll_d_ln_sigma[0]) / (2 * h)
          a2 = r0.d2_ll_d_ln_sigma2[0]
          assert abs(fd2 - a2) < 1e-5 * max(1., abs(a2)), (centric, e, s, fd2, a2)

def exercise_tail_probability():
  s = 1e-4
  for e in (1., 5., 20.):
    for centric, expected in ((False, math.exp(-e)),
                              (True, erfc(math.sqrt(e / 2.)))):
      p = xray_ext.wilson_intensity_tail_probability(
        i_obs=flex.double([e]), sig_i_obs=flex.double([s]),
        sigma=flex.double([1.]), centric_flags=flex.bool([centric]))[0]
      assert abs(p - expected) < 1e-4 * expected, (centric, e, p, expected)
  es = [0., 0.5, 1., 2., 5., 10., 20.]
  for centric in (False, True):
    p = list(xray_ext.wilson_intensity_tail_probability(
      i_obs=flex.double(es), sig_i_obs=flex.double([0.3] * len(es)),
      sigma=flex.double([1.] * len(es)),
      centric_flags=flex.bool([centric] * len(es))))
    for a, b in zip(p[:-1], p[1:]):
      assert a > b, (centric, p)
  # Negative intensities give a small-side probability in [0, 1].
  p = xray_ext.wilson_intensity_tail_probability(
    i_obs=flex.double([-0.5, -3., -3.]), sig_i_obs=flex.double([1., 1., 1.]),
    sigma=flex.double([1.] * 3), centric_flags=flex.bool([False, False, True]))
  assert 0 <= p[1] < p[0] <= 1 and 0 <= p[2] <= 1

def exercise_vectorised():
  rnd = random.Random(1)
  n = 1000
  i = [rnd.uniform(-3, 30) for _ in range(n)]
  si = [math.exp(rnd.uniform(-3, 2)) for _ in range(n)]
  sigma = [math.exp(rnd.uniform(-1, 3)) for _ in range(n)]
  cen = [rnd.random() < 0.3 for _ in range(n)]
  r = evaluate(i, si, sigma, cen)
  for k in range(n):
    one = evaluate([i[k]], [si[k]], [sigma[k]], [cen[k]])
    assert r.log_likelihood[k] == one.log_likelihood[0]
    assert r.d_ll_d_ln_sigma[k] == one.d_ll_d_ln_sigma[0]
    assert r.d2_ll_d_ln_sigma2[k] == one.d2_ll_d_ln_sigma2[0]
  t = xray_ext.wilson_intensity_tail_probability(
    flex.double(i[:50]), flex.double(si[:50]), flex.double(sigma[:50]),
    flex.bool(cen[:50]))
  for k in range(50):
    one = xray_ext.wilson_intensity_tail_probability(
      flex.double([i[k]]), flex.double([si[k]]), flex.double([sigma[k]]),
      flex.bool([cen[k]]))
    assert t[k] == one[0]

def run():
  exercise_acentric_matches_quadrature()
  exercise_centric_matches_quadrature()
  exercise_derivatives_fd()
  exercise_tail_probability()
  exercise_vectorised()
  print("OK")

if __name__ == "__main__":
  run()

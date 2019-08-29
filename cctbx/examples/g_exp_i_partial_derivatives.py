from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
import cmath
import math
from six.moves import range

def empirical_proof(g, ffp, fdp, alpha):
  # Mathematica: f = g (ffp + I fdp) Exp[I alpha]
  c = g * (ffp + 1j*fdp) * cmath.exp(1j*alpha)
  a = g * ffp * math.cos(alpha) - g * fdp * math.sin(alpha)
  b = g * ffp * math.sin(alpha) + g * fdp * math.cos(alpha)
  assert approx_equal(a, c.real)
  assert approx_equal(b, c.imag)
  #
  # Mathematica: D[f,alpha]
  d_c_d_alpha = g * (ffp + 1j*fdp) * 1j*cmath.exp(1j*alpha)
  d_a_d_alpha = g * ffp * -math.sin(alpha) - g * fdp * math.cos(alpha)
  d_b_d_alpha = g * ffp *  math.cos(alpha) - g * fdp * math.sin(alpha)
  assert approx_equal(d_a_d_alpha, d_c_d_alpha.real)
  assert approx_equal(d_b_d_alpha, d_c_d_alpha.imag)
  #
  # Mathematica: D[f,g]
  d_c_d_g = (ffp + 1j*fdp) * cmath.exp(1j*alpha)
  d_a_d_g = ffp * math.cos(alpha) - fdp * math.sin(alpha)
  d_b_d_g = ffp * math.sin(alpha) + fdp * math.cos(alpha)
  assert approx_equal(d_a_d_g, d_c_d_g.real)
  assert approx_equal(d_b_d_g, d_c_d_g.imag)
  #
  # Mathematica: D[f,ffp]
  d_c_d_ffp = g * cmath.exp(1j*alpha)
  d_a_d_ffp = g * math.cos(alpha)
  d_b_d_ffp = g * math.sin(alpha)
  assert approx_equal(d_a_d_ffp, d_c_d_ffp.real)
  assert approx_equal(d_b_d_ffp, d_c_d_ffp.imag)
  #
  # Mathematica: D[f,fdp]
  d_c_d_fdp = g * 1j * cmath.exp(1j*alpha)
  d_a_d_fdp = -g * math.sin(alpha)
  d_b_d_fdp =  g * math.cos(alpha)
  assert approx_equal(d_a_d_fdp, d_c_d_fdp.real)
  assert approx_equal(d_b_d_fdp, d_c_d_fdp.imag)

def exercise():
  for g in range(-3,4):
    for ffp in range(-3,4):
      for fdp in range(-3,4):
        for alpha_deg in range(0,360,15):
          empirical_proof(g, ffp, fdp, alpha_deg*math.pi/180)
  print("OK")

if (__name__ == "__main__"):
  exercise()

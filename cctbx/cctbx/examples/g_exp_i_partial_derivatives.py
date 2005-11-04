from libtbx.test_utils import approx_equal
import cmath
import math

def empirical_proof(g, fp, fdp, alpha):
  # Mathematica: f = g (fp + I fdp) Exp[I alpha]
  c = g * (fp + 1j*fdp) * cmath.exp(1j*alpha)
  a = g * fp * math.cos(alpha) - g * fdp * math.sin(alpha)
  b = g * fp * math.sin(alpha) + g * fdp * math.cos(alpha)
  assert approx_equal(a, c.real)
  assert approx_equal(b, c.imag)
  #
  # Mathematica: D[f,alpha]
  d_c_d_alpha = g * (fp + 1j*fdp) * 1j*cmath.exp(1j*alpha)
  d_a_d_alpha = g * fp * -math.sin(alpha) - g * fdp * math.cos(alpha)
  d_b_d_alpha = g * fp *  math.cos(alpha) - g * fdp * math.sin(alpha)
  assert approx_equal(d_a_d_alpha, d_c_d_alpha.real)
  assert approx_equal(d_b_d_alpha, d_c_d_alpha.imag)
  #
  # Mathematica: D[f,g]
  d_c_d_g = (fp + 1j*fdp) * cmath.exp(1j*alpha)
  d_a_d_g = fp * math.cos(alpha) - fdp * math.sin(alpha)
  d_b_d_g = fp * math.sin(alpha) + fdp * math.cos(alpha)
  assert approx_equal(d_a_d_g, d_c_d_g.real)
  assert approx_equal(d_b_d_g, d_c_d_g.imag)
  #
  # Mathematica: D[f,fp]
  d_c_d_fp = g * cmath.exp(1j*alpha)
  d_a_d_fp = g * math.cos(alpha)
  d_b_d_fp = g * math.sin(alpha)
  assert approx_equal(d_a_d_fp, d_c_d_fp.real)
  assert approx_equal(d_b_d_fp, d_c_d_fp.imag)
  #
  # Mathematica: D[f,fdp]
  d_c_d_fdp = g * 1j * cmath.exp(1j*alpha)
  d_a_d_fdp = -g * math.sin(alpha)
  d_b_d_fdp =  g * math.cos(alpha)
  assert approx_equal(d_a_d_fdp, d_c_d_fdp.real)
  assert approx_equal(d_b_d_fdp, d_c_d_fdp.imag)

def exercise():
  for g in xrange(-3,4):
    for fp in xrange(-3,4):
      for fdp in xrange(-3,4):
        for alpha_deg in xrange(0,360,15):
          empirical_proof(g, fp, fdp, alpha_deg/180*math.pi)
  print "OK"

if (__name__ == "__main__"):
  exercise()

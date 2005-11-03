from libtbx.test_utils import approx_equal
import cmath
import math

def empirical_proof(g, alpha):
  c = g * cmath.exp(1j*alpha)
  a = g * math.cos(alpha)
  b = g * math.sin(alpha)
  assert approx_equal(a, c.real)
  assert approx_equal(b, c.imag)
  #                                        Mathematica:
  d_c_d_alpha = g * 1j*cmath.exp(1j*alpha) # D[g Exp[I alpha], alpha]
  d_a_d_alpha = g * (-math.sin(alpha))     # D[g Cos[alpha], alpha]
  d_b_d_alpha = g * ( math.cos(alpha))     # D[g Sin[alpha], alpha]
  # assert D[Re[g Exp[I alpha],alpha]] == Re[D[g Exp[I alpha],alpha]]
  assert approx_equal(d_a_d_alpha, d_c_d_alpha.real)
  # assert D[Im[g Exp[I alpha],alpha]] == Im[D[g Exp[I alpha],alpha]]
  assert approx_equal(d_b_d_alpha, d_c_d_alpha.imag)

def exercise():
  for g in xrange(-10,10):
    for alpha_deg in xrange(0,360):
      empirical_proof(g, alpha_deg/180*math.pi)
  print "OK"

if (__name__ == "__main__"):
  exercise()

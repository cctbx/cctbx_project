from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = "--verbose" in args
  n_refl = 10
  f_calc = flex.polar(
    flex.random_double(n_refl)*10-5,
    flex.random_double(n_refl)*10-5)
  f_obs = flex.abs(f_calc) + (flex.random_double(n_refl)*2-1)
  weights = flex.random_double(n_refl)
  r = xray.targets_least_squares_residual(
    f_obs, weights, f_calc, True, 0)
  scale_factor = r.scale_factor()
  gr_ana = r.derivatives()
  gr_fin = flex.complex_double()
  eps = 1.e-6
  for i_refl in xrange(n_refl):
    gc = []
    for i_part in [0,1]:
      fc0 = f_calc[i_refl]
      ts = []
      for signed_eps in [eps,-eps]:
        if (i_part == 0):
          f_calc[i_refl] = complex(fc0.real + signed_eps, fc0.imag)
        else:
          f_calc[i_refl] = complex(fc0.real, fc0.imag + signed_eps)
        r = xray.targets_least_squares_residual(
          f_obs, weights, f_calc, False, scale_factor)
        ts.append(r.target())
      f_calc[i_refl] = fc0
      gc.append((ts[0]-ts[1])/(2*eps))
    gr_fin.append(complex(*gc))
  if (verbose):
    print "ana:", list(gr_ana)
    print "fin:", list(gr_fin)
  assert approx_equal(gr_fin, gr_ana)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])

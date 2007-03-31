from cctbx import xray
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = "--verbose" in args
  exercise_core_LS(xray.targets_least_squares_residual, verbose)
  exercise_core_LS(xray.targets_least_squares_residual_for_F_square, verbose)
  
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,11,12,85,95,100),
    space_group_symbol="P 1")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=3)
  f_calc = miller_set.array(
    data=flex.polar(
      flex.random_double(miller_set.size())*10-5,
      flex.random_double(miller_set.size())*10-5))

  obs = miller_set.array(
    data=flex.abs(f_calc.data()) + (flex.random_double(miller_set.size())*2-1),
    sigmas=flex.random_double(miller_set.size()))
  obs.set_observation_type_xray_amplitude()
  weighting = xray.weighting_schemes.unit_weighting()
  exercise_py_LS(obs, f_calc, weighting, verbose)
  
  obs = miller_set.array(
    data=flex.pow2(flex.abs(f_calc.data())) + (flex.random_double(miller_set.size())*2-1),
    sigmas=flex.random_double(miller_set.size()))
  obs.set_observation_type_xray_intensity()
  weighting = xray.weighting_schemes.unit_weighting()
  exercise_py_LS(obs, f_calc, weighting, verbose)

  
def exercise_core_LS(target_class, verbose):
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

def exercise_py_LS(obs, f_calc, weighting, verbose):
  r = xray.unified_least_squares_residual(obs, weighting=weighting)
  rt = r(f_calc, compute_derivatives=True)
  if obs.is_xray_amplitude_array():
    assert(isinstance(rt, xray.targets_least_squares_residual))
  elif obs.is_xray_intensity_array():
    assert(isinstance(rt, xray.targets_least_squares_residual_for_F_square))
  scale_factor = rt.scale_factor()
  gr_ana = rt.derivatives()
  gr_fin = flex.complex_double()
  eps = 1.e-6
  for i_refl in xrange(obs.size()):
    gc = []
    for i_part in [0,1]:
      fc0 = f_calc.data()[i_refl]
      ts = []
      for signed_eps in [eps,-eps]:
        if (i_part == 0):
          f_calc.data()[i_refl] = complex(fc0.real + signed_eps, fc0.imag)
        else:
          f_calc.data()[i_refl] = complex(fc0.real, fc0.imag + signed_eps)
        rt = r(f_calc, compute_derivatives=False, scale_factor=scale_factor)
        ts.append(rt.target())
      f_calc.data()[i_refl] = fc0
      gc.append((ts[0]-ts[1])/(2*eps))
    gr_fin.append(complex(*gc))
  if (verbose):
    print "ana:", list(gr_ana)
    print "fin:", list(gr_fin)
  assert approx_equal(gr_fin, gr_ana)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])

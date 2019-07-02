from __future__ import absolute_import, division, print_function
from cctbx import xray
from cctbx import crystal
from cctbx import miller
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import sys
from six.moves import range

def run(args):
  assert args in [[], ["--verbose"]]
  verbose = "--verbose" in args
  exercise_least_squares_residual()
  exercise_core_LS(xray.targets_least_squares_residual, verbose)
  exercise_core_LS(xray.targets_least_squares_residual_for_intensity, verbose)

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
  weighting = xray.weighting_schemes.amplitude_unit_weighting()
  exercise_py_LS(obs, f_calc, weighting, verbose)

  obs = miller_set.array(
    data=flex.norm(f_calc.data()) + (flex.random_double(miller_set.size())*2-1),
    sigmas=flex.random_double(miller_set.size()))
  obs.set_observation_type_xray_intensity()

  weighting = xray.weighting_schemes.intensity_quasi_unit_weighting()
  exercise_py_LS(obs, f_calc, weighting, verbose)

  weighting = xray.weighting_schemes.simple_shelx_weighting(a=100, b=150)
  exercise_py_LS(obs, f_calc, weighting, verbose)

  weighting = xray.weighting_schemes.shelx_weighting(a=100, b=150)
  exercise_py_LS(obs, f_calc, weighting, verbose)

  print("OK")


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
  for i_refl in range(n_refl):
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
    print("ana:", list(gr_ana))
    print("fin:", list(gr_fin))
  assert approx_equal(gr_fin, gr_ana)

def exercise_py_LS(obs, f_calc, weighting, verbose):
  weighting.computing_derivatives_wrt_f_c = True
  r = xray.unified_least_squares_residual(obs, weighting=weighting)
  rt = r(f_calc, compute_derivatives=True)
  if obs.is_xray_amplitude_array():
    assert(isinstance(rt, xray.targets_least_squares_residual))
  elif obs.is_xray_intensity_array():
    assert(isinstance(rt, xray.targets_least_squares_residual_for_intensity))
  scale_factor = rt.scale_factor()
  gr_ana = rt.derivatives()
  K = scale_factor
  w = weighting.weights
  if w is not None: w = w.deep_copy()
  dw_dfc = weighting.derivatives_wrt_f_c
  if dw_dfc is not None: dw_dfc = dw_dfc.deep_copy()

  y_o = obs.data()
  if w is None: w = flex.double(obs.size(), 1)
  sum_w_y_o_sqr = flex.sum(w * y_o * y_o)
  f_c = f_calc.data().deep_copy()
  if obs.is_xray_amplitude_array():
    y_c = flex.abs(f_c)
    der = f_c * (1/y_c)
  elif obs.is_xray_intensity_array():
    y_c = flex.norm(f_c)
    der = 2 * f_c
  gr_explicit = w*2*K*(K*y_c - y_o) * der / sum_w_y_o_sqr
  sum_w_squares = flex.sum(w*flex.pow2(K*y_c - y_o))
  assert approx_equal(gr_ana, gr_explicit)

  gr_fin = flex.complex_double()
  eps = 1.e-6
  for i_refl in range(obs.size()):
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
    print("ana:", list(gr_ana))
    print("fin:", list(gr_fin))
  if dw_dfc is None:
    assert approx_equal(gr_fin, gr_ana)
  else:
    gr_total_ana = ( gr_ana
                     + dw_dfc*(flex.pow2(K*y_c - y_o)/sum_w_y_o_sqr
                        - sum_w_squares*flex.pow2(y_o)/sum_w_y_o_sqr**2) )
    assert approx_equal(gr_fin, gr_total_ana)

def exercise_least_squares_residual():
  crystal_symmetry = crystal.symmetry(
    unit_cell=(6,3,8,90,90,90),
    space_group_symbol="P222")
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=0.7)
  f_obs = miller_set.array(
    data=flex.random_double(miller_set.size()),
    sigmas=flex.random_double(miller_set.size())*0.05)
  ls = xray.least_squares_residual(
    f_obs,
    use_sigmas_as_weights=True,
  )



if (__name__ == "__main__"):
  run(sys.argv[1:])

from __future__ import absolute_import, division, print_function
import os
import random
import time
import iotbx.pdb
import mmtbx.model
import mmtbx.f_model
import libtbx.load_env
from cctbx import xray
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal

random.seed(0)
flex.set_random_seed(0)

def get_fmodel(d_min=1.5, noise_fraction=0.05):
  """
  1yjp (P 1 21 1, 66 atoms) with simulated amplitudes: |Fcalc| from direct
  summation plus Gaussian noise of relative size noise_fraction, and the
  amplitude maximum-likelihood target ("ml").
  """
  pdb_file = libtbx.env.find_in_repositories(
    relative_path=os.path.join("iotbx", "regression", "data", "1yjp.pdb"),
    test=os.path.isfile)
  assert pdb_file is not None
  pdb_inp = iotbx.pdb.input(file_name=pdb_file)
  model = mmtbx.model.manager(model_input=pdb_inp)
  xrs = model.get_xray_structure()
  f_calc = xrs.structure_factors(d_min=d_min, algorithm="direct").f_calc()
  f_obs = abs(f_calc)
  sigmas = f_obs.data() * noise_fraction
  noise = flex.random_double(f_obs.data().size(), 1.0) - 0.5
  noise = noise * (12 ** 0.5)  # unit variance
  f_obs = f_obs.customized_copy(
    data=flex.abs(f_obs.data() + sigmas * noise), sigmas=sigmas)
  r_free_flags = f_obs.generate_r_free_flags(fraction=0.1)
  sfg_params = mmtbx.f_model.sf_and_grads_accuracy_master_params.extract()
  sfg_params.algorithm = "direct"
  fmodel = mmtbx.f_model.manager(
    f_obs=f_obs,
    r_free_flags=r_free_flags,
    xray_structure=xrs,
    target_name="ml",
    sf_and_grads_accuracy_params=sfg_params)
  return model, fmodel

def mlf_arguments(fmodel):
  alpha, beta = fmodel.alpha_beta()
  return dict(
    f_obs=fmodel.f_obs().data(),
    r_free_flags=fmodel.r_free_flags().data(),
    alpha=alpha.data(),
    beta=beta.data(),
    scale_factor=fmodel.scale_ml_wrapper(),
    epsilons=fmodel.epsilons,
    centric_flags=fmodel.f_obs().centric_flags().data())

def exact_i1_over_i0(x):
  """I1(x)/I0(x) from the power series, summed in log space."""
  import math
  if x == 0:
    return 0.0
  n_terms = int(2 * x) + 200
  l = math.log(x / 2)
  s0, s1 = 0.0, 0.0
  for k in range(n_terms):
    s0 += math.exp(2 * k * l - 2 * math.lgamma(k + 1) - x)
    s1 += math.exp((2 * k + 1) * l - math.lgamma(k + 1)
                   - math.lgamma(k + 2) - x)
  return s1 / s0

def exercise_mlf_gradients_match_phenix():
  """
  The radial first derivative of mmtbx.refinement.curvatures reproduces
  the gradients of cctbx.xray.mlf_target_and_gradients (to the accuracy of
  the Bessel-ratio polynomial fit in scitbx/math/bessel.h).
  """
  from mmtbx.refinement import curvatures
  model, fmodel = get_fmodel()
  args = mlf_arguments(fmodel)
  f_calc = fmodel.f_model().data()
  k = args["scale_factor"]
  if k <= 0:
    k = 1.0
  a = args["alpha"] * k
  sigma = args["epsilons"] * args["beta"] * (k * k)
  t1, t2 = curvatures.rice_radial_derivatives(
    f_obs=args["f_obs"], f_calc_abs=flex.abs(f_calc), a=a, sigma=sigma,
    centric_flags=args["centric_flags"])
  work = ~args["r_free_flags"]
  g = curvatures.gradients_from_radial(f_calc=f_calc, t1=t1).select(work) \
    * (1.0 / work.count(True))
  g_phenix = xray.mlf_target_and_gradients(
    f_calc=f_calc, compute_gradients=True, **args).gradients_work()
  assert g.size() == g_phenix.size()
  scale = flex.max(flex.abs(g_phenix))
  assert scale > 0
  assert flex.max(flex.abs(g - g_phenix)) / scale < 1.e-6

def exercise_mlf_radial_second_derivative():
  """
  The radial second derivative t2 of rice_radial_derivatives matches a
  central finite difference of the exact radial first derivative (series
  Bessel ratio), for acentric and centric reflections of 1yjp.
  """
  import math
  from mmtbx.refinement import curvatures
  x_values = []
  for noise_fraction in (0.05, 0.5):
    model, fmodel = get_fmodel(noise_fraction=noise_fraction)
    args = mlf_arguments(fmodel)
    f_calc = fmodel.f_model().data()
    fc_abs = flex.abs(f_calc)
    k = args["scale_factor"]
    if k <= 0:
      k = 1.0
    a = args["alpha"] * k
    sigma = args["epsilons"] * args["beta"] * (k * k)
    centric = args["centric_flags"]
    t1, t2 = curvatures.rice_radial_derivatives(
      f_obs=args["f_obs"], f_calc_abs=fc_abs, a=a, sigma=sigma,
      centric_flags=centric)
    def t1_exact(i, fc):
      fo = args["f_obs"][i]
      if centric[i]:
        x = a[i] * fo * fc / sigma[i]
        return a[i] * a[i] * fc / sigma[i] \
             - a[i] * fo / sigma[i] * math.tanh(x)
      x = 2 * a[i] * fo * fc / sigma[i]
      return 2 * a[i] * a[i] * fc / sigma[i] \
           - 2 * a[i] * fo / sigma[i] * exact_i1_over_i0(x)
    # the Bessel argument, weakest and strongest reflections of each kind
    x = 2 * a * args["f_obs"] * fc_abs / sigma
    picked = []
    for flag in (False, True):
      indices = [i for i in range(x.size()) if centric[i] == flag]
      indices.sort(key=lambda i: x[i])
      assert len(indices) > 10
      picked.extend(indices[:5] + indices[-5:])
    for i in picked:
      fc = fc_abs[i]
      h = 1.e-4 * fc
      fd = (t1_exact(i, fc + h) - t1_exact(i, fc - h)) / (2 * h)
      # t1 is a small difference of two large terms; tolerances scale
      # with the larger term, 2 a^2 fc / sigma. t1 uses the Bessel-ratio
      # polynomial fit of scitbx/math/bessel.h, accurate to about 1e-6.
      term = 2 * a[i] * a[i] * fc / sigma[i]
      assert approx_equal(t1[i], t1_exact(i, fc), eps=5.e-6 * term)
      assert approx_equal(t2[i], fd, eps=1.e-5 * term / fc)
      x_values.append(x[i])
  assert max(x_values) > 100  # strong reflections are represented
  assert min(x_values) < 1    # and weak ones

def exercise_least_squares_limit():
  """
  hessians_from_radial reproduces the second derivatives of the amplitude
  least-squares target in cctbx/xray/targets/least_squares.h, with the
  radial derivatives of w (fo - k |Fc|)^2 / sum(w fo^2) as input.
  """
  from mmtbx.refinement import curvatures
  model, fmodel = get_fmodel()
  f_obs = fmodel.f_obs().data()
  f_calc = fmodel.f_model().data()
  weights = flex.double(f_obs.size(), 1.0)
  k = 1.3
  ls = xray.targets_least_squares(
    compute_scale_using_all_data=False,
    obs_type="F",
    obs=f_obs,
    weights=weights,
    r_free_flags=None,
    f_calc=f_calc,
    derivatives_depth=2,
    scale_factor=k)
  fc_abs = flex.abs(f_calc)
  norm = flex.sum(weights * f_obs * f_obs)
  t1 = -2 * weights * k * (f_obs - k * fc_abs) / norm
  t2 = 2 * weights * k * k / norm
  hessians = curvatures.hessians_from_radial(f_calc=f_calc, t1=t1, t2=t2)
  assert approx_equal(ls.gradients_work(),
    curvatures.gradients_from_radial(f_calc=f_calc, t1=t1), eps=1.e-10)
  assert approx_equal(hessians, ls.hessians_work(), eps=1.e-10)

def run():
  exercise_least_squares_limit()
  exercise_mlf_gradients_match_phenix()
  exercise_mlf_radial_second_derivative()
  print("OK")

if __name__ == "__main__":
  t0 = time.time()
  run()
  print("Time: %6.2f" % (time.time() - t0))

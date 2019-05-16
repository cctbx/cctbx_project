from __future__ import absolute_import, division, print_function
from cctbx.examples.g_exp_i_alpha_derivatives \
  import parameters, gradients, pack_gradients, g_exp_i_alpha_sum
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from libtbx.test_utils import approx_equal
import random
import math
import copy
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

random.seed(0)

def d_target_d_params_finite(obs, params, eps=1.e-8):
  result = []
  params_eps = copy.deepcopy(params)
  for i_param in range(len(params)):
    dx = []
    for ix in range(4):
      ts = []
      for signed_eps in [eps, -eps]:
        pi_eps = params[i_param].as_list()
        pi_eps[ix] += signed_eps
        params_eps[i_param] = parameters(*pi_eps)
        exp_sum = g_exp_i_alpha_sum(params=params_eps)
        target = least_squares(obs=obs, calc=exp_sum.f())
        ts.append(target.f())
      dx.append((ts[0]-ts[1])/(2*eps))
    result.append(gradients(*dx))
    params_eps[i_param] = params[i_param]
  return result

def d2_target_d_params_finite(obs, params, eps=1.e-8):
  result = []
  params_eps = copy.deepcopy(params)
  for i_param in range(len(params)):
    for ix in range(4):
      gs = []
      for signed_eps in [eps, -eps]:
        pi_eps = params[i_param].as_list()
        pi_eps[ix] += signed_eps
        params_eps[i_param] = parameters(*pi_eps)
        exp_sum = g_exp_i_alpha_sum(params=params_eps)
        target = least_squares(obs=obs, calc=exp_sum.f())
        dp = exp_sum.d_target_d_params(target=target)
        gs.append(pack_gradients(dp))
      result.append([(gp-gm)/(2*eps) for gp,gm in zip(gs[0],gs[1])])
    params_eps[i_param] = params[i_param]
  return result

def compare_analytical_and_finite(obs, params, out):
  grads_fin = d_target_d_params_finite(obs=obs, params=params)
  print("grads_fin:", pack_gradients(grads_fin), file=out)
  exp_sum = g_exp_i_alpha_sum(params=params)
  target = least_squares(obs=obs, calc=exp_sum.f())
  grads_ana = exp_sum.d_target_d_params(target=target)
  print("grads_ana:", pack_gradients(grads_ana), file=out)
  assert approx_equal(pack_gradients(grads_ana), pack_gradients(grads_fin))
  curvs_fin = d2_target_d_params_finite(obs=obs, params=params)
  print("curvs_fin:", curvs_fin, file=out)
  curvs_ana = list(exp_sum.d2_target_d_params(target=target))
  print("curvs_ana:", curvs_ana, file=out)
  assert approx_equal(curvs_ana, curvs_fin)
  print(file=out)

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  for n_params in range(2,5):
    for i_trial in range(5):
      params = []
      for i in range(n_params):
        params.append(parameters(
          g=(random.random()-0.5)*2,
          ffp=(random.random()-0.5)*2,
          fdp=(random.random()-0.5)*2,
          alpha=2*math.pi*random.random()))
      exp_sum = g_exp_i_alpha_sum(params=params)
      obs = abs(exp_sum.f())
      compare_analytical_and_finite(
        obs=obs,
        params=params,
        out=out)
      compare_analytical_and_finite(
        obs=obs*(random.random()+0.5),
        params=params,
        out=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

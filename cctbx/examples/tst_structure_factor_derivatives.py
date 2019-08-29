from __future__ import absolute_import, division, print_function
from cctbx.examples.structure_factor_derivatives \
  import parameters, gradients, pack_gradients, structure_factor
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from libtbx.test_utils import approx_equal
import random
import copy
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

random.seed(0)

def d_target_d_params_finite(obs, hkl, d_star_sq, params, eps=1.e-8):
  result = []
  params_eps = copy.deepcopy(params)
  for i_param in range(len(params)):
    dx = []
    for ix in range(7):
      ts = []
      for signed_eps in [eps, -eps]:
        pi_eps = params[i_param].as_list()
        pi_eps[ix] += signed_eps
        params_eps[i_param] = parameters(pi_eps[:3], *pi_eps[3:])
        sf = structure_factor(hkl=hkl, d_star_sq=d_star_sq, params=params_eps)
        target = least_squares(obs=obs, calc=sf.f())
        ts.append(target.f())
      dx.append((ts[0]-ts[1])/(2*eps))
    result.append(gradients(dx[:3], *dx[3:]))
    params_eps[i_param] = params[i_param]
  return result

def d2_target_d_params_finite(obs, hkl, d_star_sq, params, eps=1.e-8):
  result = []
  params_eps = copy.deepcopy(params)
  for i_param in range(len(params)):
    for ix in range(7):
      gs = []
      for signed_eps in [eps, -eps]:
        pi_eps = params[i_param].as_list()
        pi_eps[ix] += signed_eps
        params_eps[i_param] = parameters(pi_eps[:3], *pi_eps[3:])
        sf = structure_factor(hkl=hkl, d_star_sq=d_star_sq, params=params_eps)
        target = least_squares(obs=obs, calc=sf.f())
        dp = sf.d_target_d_params(target=target)
        gs.append(pack_gradients(dp))
      result.append([(gp-gm)/(2*eps) for gp,gm in zip(gs[0],gs[1])])
    params_eps[i_param] = params[i_param]
  return result

def compare_analytical_and_finite(obs, hkl, d_star_sq, params, out):
  grads_fin = d_target_d_params_finite(
    obs=obs, hkl=hkl, d_star_sq=d_star_sq, params=params)
  print("grads_fin:", pack_gradients(grads_fin), file=out)
  sf = structure_factor(hkl=hkl, d_star_sq=d_star_sq, params=params)
  target = least_squares(obs=obs, calc=sf.f())
  grads_ana = sf.d_target_d_params(target=target)
  print("grads_ana:", pack_gradients(grads_ana), file=out)
  assert approx_equal(pack_gradients(grads_ana), pack_gradients(grads_fin))
  curvs_fin = d2_target_d_params_finite(
    obs=obs, hkl=hkl, d_star_sq=d_star_sq, params=params)
  print("curvs_fin:", curvs_fin, file=out)
  curvs_ana = sf.d2_target_d_params(target=target)
  print("curvs_ana:", curvs_ana, file=out)
  assert approx_equal(curvs_ana, curvs_fin, 1.e-5)
  print(file=out)

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  hkl = (1,2,3)
  d_star_sq = 1e-3
  for n_params in range(2,5):
    for i_trial in range(5):
      params = []
      for i in range(n_params):
        params.append(parameters(
          xyz=[random.random() for i in range(3)],
          u=random.random()*0.1,
          w=random.random(),
          fp=(random.random()-0.5)*2,
          fdp=(random.random()-0.5)*2))
      sf = structure_factor(hkl=hkl, d_star_sq=d_star_sq, params=params)
      obs = abs(sf.f())
      compare_analytical_and_finite(
        obs=obs,
        hkl=hkl,
        d_star_sq=d_star_sq,
        params=params,
        out=out)
      compare_analytical_and_finite(
        obs=obs*(random.random()+0.5),
        hkl=hkl,
        d_star_sq=d_star_sq,
        params=params,
        out=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

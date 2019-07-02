from __future__ import absolute_import, division, print_function
from cctbx.examples.exp_i_alpha_derivatives \
  import least_squares, exp_i_alpha_sum
from libtbx.test_utils import approx_equal
import random
import math
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

random.seed(0)

def d_target_d_alphas_finite(obs, alphas, eps=1.e-8):
  result = []
  for i_alpha in range(len(alphas)):
    alphas_eps = list(alphas)
    ts = []
    for signed_eps in [eps, -eps]:
      alphas_eps[i_alpha] = alphas[i_alpha] + signed_eps
      exp_sum = exp_i_alpha_sum(alphas=alphas_eps)
      target = least_squares(obs=obs, calc=exp_sum.f())
      ts.append(target.f())
    result.append((ts[0]-ts[1])/(2*eps))
  return result

def d2_target_d_alphas_finite(obs, alphas, eps=1.e-8):
  result = []
  for i_alpha in range(len(alphas)):
    alphas_eps = list(alphas)
    gs = []
    for signed_eps in [eps, -eps]:
      alphas_eps[i_alpha] = alphas[i_alpha] + signed_eps
      exp_sum = exp_i_alpha_sum(alphas=alphas_eps)
      target = least_squares(obs=obs, calc=exp_sum.f())
      dalphas = exp_sum.d_target_d_alphas(target=target)
      gs.append(dalphas)
    result.append([(gp-gm)/(2*eps) for gp,gm in zip(gs[0],gs[1])])
  return result

def compare_analytical_and_finite(obs, alphas, out):
  grads_fin = d_target_d_alphas_finite(obs=obs, alphas=alphas)
  print("grads_fin:", grads_fin, file=out)
  exp_sum = exp_i_alpha_sum(alphas=alphas)
  target = least_squares(obs=obs, calc=exp_sum.f())
  grads_ana = exp_sum.d_target_d_alphas(target=target)
  print("grads_ana:", grads_ana, file=out)
  assert approx_equal(grads_ana, grads_fin)
  curvs_fin = d2_target_d_alphas_finite(obs=obs, alphas=alphas)
  print("curvs_fin:", curvs_fin, file=out)
  curvs_ana = exp_sum.d2_target_d_alphas(target=target)
  print("curvs_ana:", curvs_ana, file=out)
  assert approx_equal(curvs_ana, curvs_fin)
  print(file=out)

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  for n_alphas in range(2,5):
    for i_trial in range(5):
      alphas = [2*math.pi*random.random() for i in range(n_alphas)]
      exp_sum = exp_i_alpha_sum(alphas=alphas)
      obs = abs(exp_sum.f())
      compare_analytical_and_finite(
        obs=obs,
        alphas=alphas,
        out=out)
      compare_analytical_and_finite(
        obs=obs*(random.random()+0.5),
        alphas=alphas,
        out=out)
  for obs in [0, 0.1]:
    for calc in [0j, 1.e-200j]:
      target = least_squares(obs=obs, calc=calc)
      assert target.f() == obs**2
      assert target.da() == 0
      assert target.db() == 0
      if (obs == 0):
        assert target.daa() == 2
        assert target.dbb() == 2
        assert target.dab() == 0
      else:
        assert target.daa() == -1.e160
        assert target.dbb() == -1.e160
        assert target.dab() == 1.e160
    calc = 1+1j
    while (calc != 0):
      target = least_squares(obs=obs, calc=calc)
      # exercise numerical stability without checking the results
      target.f()
      target.da()
      target.db()
      target.daa()
      target.dbb()
      target.dab()
      calc /= 2
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

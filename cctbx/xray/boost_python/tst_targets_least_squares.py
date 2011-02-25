from cctbx.xray import ext
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.math_utils import iceil
from itertools import count
import sys

class random_inputs(object):

  def __init__(O, mt, n_refl, obs_type):
    O.obs_type = obs_type
    O.obs = mt.random_double(size=n_refl)
    O.weights = mt.random_double(size=n_refl)
    rff = flex.bool(max(1,iceil(n_refl*0.6)), False)
    rff.resize(n_refl, True)
    O.r_free_flags = rff.select(mt.random_permutation(size=n_refl))
    O.scale_factor = 1 + mt.random_double()
    O.a = mt.random_double(size=n_refl)
    O.b = mt.random_double(size=n_refl)

  def get(O, derivatives_depth=0):
    return ext.targets_least_squares(
      compute_scale_using_all_data=False,
      obs_type=O.obs_type,
      obs=O.obs,
      weights=O.weights,
      r_free_flags=O.r_free_flags,
      f_calc=flex.complex_double(O.a, O.b),
      derivatives_depth=derivatives_depth,
      scale_factor=O.scale_factor)

  def gradients_work_fd(O, eps=1.e-6):
    result = flex.complex_double()
    for ih,a,b,f in zip(count(), O.a, O.b, O.r_free_flags):
      if (f): continue
      def fd(x, Ox):
        fs = []
        for signed_eps in [eps, -eps]:
          Ox[ih] = x + signed_eps
          fs.append(O.get(derivatives_depth=0).target_work())
        Ox[ih] = x
        return (fs[0]-fs[1])/(2*eps)
      result.append(complex(fd(a, O.a), fd(b, O.b)))
    return result

  def hessians_work_fd(O, eps=1.e-6):
    result = flex.vec3_double()
    for ih,a,b,f in zip(count(), O.a, O.b, O.r_free_flags):
      if (f): continue
      def fd(x, Ox, ri):
        fs = []
        for signed_eps in [eps, -eps]:
          Ox[ih] = x + signed_eps
          ga = O.get(derivatives_depth=1).gradients_work()
          fs.append(getattr(ga[len(result)], ri))
        Ox[ih] = x
        return (fs[0]-fs[1])/(2*eps)
      daa = fd(a, O.a, "real")
      dbb = fd(b, O.b, "imag")
      dab = fd(a, O.a, "imag")
      dba = fd(b, O.b, "real")
      assert approx_equal(dab, dba)
      result.append((daa, dbb, dab))
    return result

def exercise_random(n_trials=10, n_refl=30):
  mt = flex.mersenne_twister(seed=0)
  for i_trial in xrange(n_trials):
    for obs_type in ["F", "I"]:
      ri = random_inputs(mt=mt, n_refl=n_refl, obs_type=obs_type)
      ls = ri.get(derivatives_depth=2)
      ga = ls.gradients_work()
      gf = ri.gradients_work_fd()
      assert approx_equal(ga, gf)
      ca = ls.hessians_work()
      cf = ri.hessians_work_fd()
      assert approx_equal(ca, cf)

def exercise_singular():
  obs = flex.double([1.234])
  weights = flex.double([2.345])
  r_free_flags = flex.bool([False])
  scale_factor = 3.456
  a = flex.double([0])
  b = flex.double([0])
  for obs_type in ["F", "I"]:
    ls = ext.targets_least_squares(
      compute_scale_using_all_data=False,
      obs_type=obs_type,
      obs=obs,
      weights=weights,
      r_free_flags=r_free_flags,
      f_calc=flex.complex_double(a, b),
      derivatives_depth=2,
      scale_factor=scale_factor)
    assert list(ls.gradients_work()) == [0j]
    assert list(ls.hessians_work()) == [(1,1,1)]

def run(args):
  assert len(args) == 0
  exercise_random()
  exercise_singular()
  print "OK"

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

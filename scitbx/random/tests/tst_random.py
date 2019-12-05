from __future__ import absolute_import, division, print_function
import libtbx.utils
from libtbx.test_utils import approx_equal
import scitbx.random
from scitbx.math import basic_statistics
from scitbx.array_family import flex
import itertools
import math

import libtbx.load_env
from six.moves import range
boost_version = libtbx.env.boost_version

def exercise_distributions():
  from scitbx.random import normal_distribution
  n = normal_distribution()
  assert (n.mean, n.sigma) == (0, 1)
  n = normal_distribution(mean=5, sigma=10)
  assert (n.mean, n.sigma) == (5, 10)

def exercise_variate_generators():
  from scitbx.random \
       import variate, normal_distribution, bernoulli_distribution, \
              gamma_distribution, poisson_distribution
  for i in range(10):
    scitbx.random.set_random_seed(0)
    g = variate(normal_distribution())
    assert approx_equal(g(), -0.917787219374)
    assert approx_equal(g(10),
      (1.21838707856, 1.732426915,
        0.838038157555, -0.296895169923,
        0.246451144946, -0.635474652255,
        -0.0980626986425, 0.36458295417,
        0.534073780268, -0.665073136294))

  stat = basic_statistics(flex.double(itertools.islice(g, 1000000)))
  assert approx_equal(stat.mean,            0, eps=0.005)
  assert approx_equal(stat.biased_variance, 1, eps=0.005)
  assert approx_equal(stat.skew,            0, eps=0.005)
  assert approx_equal(stat.kurtosis,        3, eps=0.005)

  bernoulli_seq = variate(bernoulli_distribution(0.1))
  for b in itertools.islice(bernoulli_seq, 10):
    assert b in (True, False)
  bernoulli_sample = flex.bool(itertools.islice(bernoulli_seq, 10000))
  assert approx_equal(
    bernoulli_sample.count(True)/len(bernoulli_sample),
    0.1,
    eps = 0.01)

  # Boost 1.64 changes the exponential distribution to use Ziggurat algorithm
  scitbx.random.set_random_seed(0)
  g = variate(gamma_distribution())
  if (boost_version < 106400):
    assert approx_equal(g(), 0.79587450456577546)
    assert approx_equal(g(2), (0.89856038848394115, 1.2559307580473893))
  else:
    assert approx_equal(g(), 0.864758191783)
    assert approx_equal(g(2), (1.36660841837, 2.26740986094))
  stat = basic_statistics(flex.double(itertools.islice(g, 1000000)))
  assert approx_equal(stat.mean,            1, eps=0.005)
  assert approx_equal(stat.skew,            2, eps=0.01)
  assert approx_equal(stat.biased_variance, 1, eps=0.005)
  scitbx.random.set_random_seed(0)
  g = variate(gamma_distribution(alpha=2, beta=3))
  assert approx_equal(g(), 16.670850592722729)
  assert approx_equal(g(2), (10.03662877519449, 3.9357158398972873))
  stat = basic_statistics(flex.double(itertools.islice(g, 1000000)))
  assert approx_equal(stat.mean,            6, eps=0.005)
  assert approx_equal(stat.skew,            2/math.sqrt(2), eps=0.05)
  assert approx_equal(stat.biased_variance, 18, eps=0.05)

  mean = 10.0
  pv = variate(poisson_distribution(mean))
  draws = pv(1000000).as_double()
  m = flex.mean(draws)
  v = flex.mean(draws*draws) - m*m
  assert approx_equal(m,mean,eps=0.05)
  assert approx_equal(v,mean,eps=0.05)

def run():
  libtbx.utils.show_times_at_exit()
  exercise_distributions()
  exercise_variate_generators()

if __name__ == '__main__':
  run()

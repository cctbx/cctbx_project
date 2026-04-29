from __future__ import absolute_import, division, print_function
from six.moves import range

#  get_prob_more_than_z_n_tries

# Simple function to calculate probability that a Gaussian-distributed sample
#  with mean zero and SD=1 (such as a Z-score) will be less than Z for all
#  N tries.

# NOTE: Not exact...see exercise_1()

# See info in http://www.panix.com/~kts/Thesis/extreme/extreme2.html
# TT 2019-01-02
# This is based on the Type I Gumbel distribution


# Usage: pp=get_prob_more_than_z_n_tries(z,n)
#  z = best Z-score you obtained
#  n = number of tries
#  pp = result = probability you would have gotten a higher score than this
#   by chance.

# for small values, calculate directly

def p_of_none_greater(z,n):
  from scitbx.math import erf
  x=z/(2**0.5)
  p_one_greater=0.5*(1-erf(x)) # one-tailed prob of Z > z
  p_not_one_greater=1-p_one_greater # prob Z is not > z
  p_all_not_greater=p_not_one_greater**n # prob all N not >z
  p_at_least_one_greater=1-p_all_not_greater # prob at least one >z
  return p_at_least_one_greater

import math

def exp(x):
  return math.exp(max(-20,min(20,x)))

def pdf_of_x(x,a,b):
  z=(x-b)/a
  return exp(-z - exp(-z))/a

def cdf_of_x(x,a,b):
  z=(x-b)/a
  return exp(-exp(-z))

def log(n):
  return math.log(max(1.e-30,min(1.e+30,n)))

def get_an_bn(n):
  bn=max(0.001,2*log(n)-log(log(n))-log(4.*3.14159))**0.5
  an=1/bn
  return an,bn

def get_prob_more_than_z_n_tries(z,n,n_direct=10):
  return 1-get_prob_less_than_z_n_tries(z,n,n_direct=n_direct)

def get_prob_less_than_z_n_tries(z,n,n_direct=10):
  if n <= n_direct:
    return 1-p_of_none_greater(z,n)
  else:
   an,bn=get_an_bn(n)
   return cdf_of_x(z,an,bn)

class sample_extreme_z:
  """Return samples from extreme probability distribution for Gaussian
     variable z (mean 0, variance 1)  with n = n_tries

    Call with number of tries (n_tries) and optional random seed.
    Then sample with the sample() function to return individual samples
      from the extreme probability distribution

    Optionally use the sample_extreme_from_gaussian() function to directly
      sample (much slower).

    For values of n_tries under n_direct, uses the sample_extreme_from_gaussian
      function always. (Gumbel function very inaccurate for small n_tries,
      pretty ok with n_tries 10 or greater, use 100 as this is fast anyhow)

"""

  def __init__(self, n_tries, random_seed = 1779147, n_direct = 100):
    from scitbx.array_family import flex
    import random
    self.n_tries = n_tries
    self.random_seed = random_seed
    self.n_direct = max(1, n_direct)
    flex.set_random_seed(self.random_seed)
    random.seed(self.random_seed)
    import numpy as np
    np.random.seed(self.random_seed)

    # Get Gumbel parameters: an is scale (beta), bn is location (mu)
    self.an,self.bn=get_an_bn(self.n_tries)

  def sample(self, n_sample = 1):
    """ Return n_sample samples from the extreme probability distribution
        of z with n = n_tries.
      This will yield samples that have a distribution that is similar to
      sample_extreme_from_gaussian (but is much faster).

      If n_sample == 1, return a float, otherwise a flex.double array

    """
    from scitbx.array_family import flex

    if self.n_tries<= self.n_direct:
      values = flex.double()
      for i in range(n_sample):
        values.append(self.sample_extreme_from_gaussian())
      if n_sample == 1:
        return values[0]
      else:
        return values

    else:
      import numpy as np
      values = np.random.gumbel(loc=self.bn, scale=self.an, size=n_sample)
      return flex.double(values)

  def sample_extreme_from_gaussian(self):
    """ Return extreme value of z for n_tries by direct sampling"""
    from scitbx.array_family import flex
    import random
    values = flex.double()
    for i in range(self.n_tries):
      values.append(random.gauss(0,1))
    return values.min_max_mean().max

def exercise_sample(random_seed = 1779147):
  for n in (1,10, 20, 30, 40, 50, 100, 1000, 10000):
    s = sample_extreme_z(n, random_seed = random_seed)
    from scitbx.array_family import flex
    gumbel= s.sample(n_sample = 100)
    sample_extreme = flex.double()
    for k in range(100):
      sample_extreme.append(s.sample_extreme_from_gaussian())

    print("N: %s  Gumbel %.2f +/- %.2f  Sample: %.2f +/- %.2f " %(
      n,
      gumbel.min_max_mean().mean,
      gumbel.standard_deviation_of_the_sample(),
      sample_extreme.min_max_mean().mean,
      sample_extreme.standard_deviation_of_the_sample(),))
    d1 = abs(gumbel.min_max_mean().mean - sample_extreme.min_max_mean().mean)
    d2 = abs(gumbel.standard_deviation_of_the_sample() -
         sample_extreme.standard_deviation_of_the_sample())
    assert d1 < 0.2
    assert d2 < 0.2

def exercise():
  from libtbx.test_utils import approx_equal
  for n in [10,100,1000,10000,]:
    an,bn=get_an_bn(n)
    print("\nN=",n," expected maximum Z:",bn)
    last_pdf_u=None
    sum=0
    for i in range(1000000):
      x=i/10000
      pdf_u=pdf_of_x(x,an,bn)
      sum+=pdf_u*0.0001
      if last_pdf_u and pdf_u < last_pdf_u:
        print("N %s Z: %.2f Scale: %.2f Expected: %.2f  " %(
              n,x,an,bn,))
        print("P(Z): %.2f CDF(sum): %.2f  CDF(Z): %.2f " %(
          pdf_u,sum,cdf_of_x(x,an,bn)))
        print("P(z>Z): %.2f  P(z>Z),direct: %.2f" %(
            get_prob_more_than_z_n_tries(x,n),p_of_none_greater(x,n)))
        assert approx_equal(bn,x,eps=0.01)
        if n > 10:
          assert approx_equal(
          1-cdf_of_x(x,an,bn),
          get_prob_more_than_z_n_tries(x,n),eps=0.01)
        if n >=100:
          assert approx_equal(sum,cdf_of_x(x,an,bn),eps=0.01)
          assert approx_equal(get_prob_more_than_z_n_tries(x,n),
            p_of_none_greater(x,n),eps=0.1)
        break
      last_pdf_u=pdf_u

def exercise_1():
  from scitbx.array_family import flex
  a=flex.double()
  b=flex.double()
  for i in range(1,100):
    for xx in range(1,1000):
      x=xx/10
      n=10**i
      if p_of_none_greater(x,n)> 0.001 and p_of_none_greater(x,n)< 0.999:
        a.append(p_of_none_greater(x,n))
        b.append(get_prob_more_than_z_n_tries(x,n))
  print(a.size())
  print(flex.linear_correlation(a,b).coefficient())
  assert flex.linear_correlation(a,b).coefficient() > 0.99


def run(args):
  n=float(args[0])
  z=float(args[1])
  print("N=",n," Z=",z," p( z<=Z) in N tries: ",get_prob_less_than_z_n_tries(z,n))

if __name__=="__main__":
  import sys
  args=sys.argv[1:]
  if 'exercise' in args:
    exercise()
    exercise_1()
    exercise_sample()
  elif len(args) == 2:
    run(args)
  elif len(args) == 1:
    exercise_sample(random_seed=int(args[0]))

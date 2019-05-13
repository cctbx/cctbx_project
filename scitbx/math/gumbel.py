from __future__ import division
from __future__ import print_function

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

def exercise():
  from libtbx.test_utils import approx_equal
  for n in [10,100,1000,10000,]:
    an,bn=get_an_bn(n)
    print("\nN=",n," expected maximum Z:",bn)
    last_pdf_u=None
    sum=0
    for i in xrange(1000000):
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
  for i in xrange(1,100):
    for xx in xrange(1,1000):
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
  else:
    run(args)

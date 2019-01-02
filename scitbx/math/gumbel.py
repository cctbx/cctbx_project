from __future__ import division

#  get_prob_more_than_z_n_tries

# Simple function to calculate probability that a Gaussian-distributed sample
#  with mean zero and SD=1 (such as a Z-score) will be less than Z for all
#  N tries.
# See info in http://www.panix.com/~kts/Thesis/extreme/extreme2.html
# TT 2019-01-02
# This is based on the Type I Gumbel distribution


# Usage: pp=get_prob_more_than_z_n_tries(z,n)
#  z = best Z-score you obtained
#  n = number of tries
#  pp = result = probability you would have gotten a higher score than this
#   by chance.

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
  bn=(2*log(n)-log(log(n))-log(4.*3.14159))**0.5
  an=1/bn
  return an,bn

def get_prob_more_than_z_n_tries(z,n):
  return 1-get_prob_less_than_z_n_tries(z,n)

def get_prob_less_than_z_n_tries(z,n):
 an,bn=get_an_bn(n)
 return cdf_of_x(z,an,bn)

def exercise():
  from libtbx.test_utils import approx_equal
  for n in [100,1000,10000,]:
    an,bn=get_an_bn(n)
    print "N=",n," expected maximum Z:",bn
    last_pdf_u=None
    sum=0
    for i in xrange(1000000):
      x=i/10000
      pdf_u=pdf_of_x(x,an,bn)
      sum+=pdf_u*0.0001
      if last_pdf_u and pdf_u < last_pdf_u:
        print n,x,an,bn,pdf_u,sum,cdf_of_x(x,an,bn),\
           get_prob_less_than_z_n_tries(x,n)
        assert approx_equal(bn,x,eps=0.01)
        assert approx_equal(
          cdf_of_x(x,an,bn),
          get_prob_less_than_z_n_tries(x,n),eps=0.01)
        assert approx_equal(sum,cdf_of_x(x,an,bn),eps=0.01)
        break
      last_pdf_u=pdf_u

def run(args):
  n=float(args[0])
  z=float(args[1])
  print "N=",n," Z=",z," p( z<=Z) in N tries: ",get_prob_less_than_z_n_tries(z,n)

if __name__=="__main__":
  import sys
  args=sys.argv[1:]
  if 'exercise' in args:
    exercise()
  else:
    run(args)

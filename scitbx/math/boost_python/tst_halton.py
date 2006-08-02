import scitbx.math as sm
import math
from libtbx.test_utils import approx_equal

"""
This routine tests the halton sequence that can be used for
numerical integration (or uniformish samping of a hypercube for that manner).
Specifically:
\int_0^{1} f(y) dy = lim(n->infty) (1/n) Sum f(x_n)
where x_n is the n-th halton sequence number.
If we want to extend that range, take into account jacobians!
i.e.:
\int_{-a}^{a} f(y) dy \approx 2*a/n Sum f( (1-2.0*a)*x_n ) 

See also

http://en.wikipedia.org/wiki/Quasi-random_sequence

for more information on quasi random sequences.

"""

def test_halton_sequence_1(n_points, high_limit=1.0):
  h_gen = sm.halton(1)
  result = 0.0
  for ii in xrange(n_points):
    x = h_gen.nth_given_base(5, ii)*high_limit 
    y = h_gen.nth_given_base(7, ii)*high_limit
    result += math.exp( float(-x - y) )
    tmp = ( math.exp(-2.0*high_limit) )*( math.exp(high_limit)-1.0 )**2.0
  result /= float(n_points)
  result *= high_limit*high_limit
  tmp = ( math.exp(-2.0*high_limit) )*( math.exp(high_limit)-1.0 )**2.0
  assert approx_equal( result/tmp, 1.0 ,eps=1e-2) 


def test_halton_sequence_2(n_points, high_limit=5.0):
  norm = math.pi*2.0
  h_gen = sm.halton(1)
  result = 0.0
  for ii in xrange(n_points):
    x = (1.0-2.0*h_gen.nth_given_base(5, ii))*high_limit
    y = (1.0-2.0*h_gen.nth_given_base(7, ii))*high_limit
    result += math.exp( -(x*x+y*y)*0.5 )/norm
  tmp = high_limit*high_limit*4.0*result/(float(n_points))
  assert approx_equal(tmp,1.0,eps=1e-2)

def visualise(n=500):
  h_gen=sm.halton(2)
  for ii in xrange(n):
    print h_gen.nth_given_base(2,ii), h_gen.nth_given_base(3,ii)

def run():
  test_halton_sequence_1(1000,2)
  test_halton_sequence_2(1000,10)

if (__name__ == "__main__"):
  run()
  print "OK"



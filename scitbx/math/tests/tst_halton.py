from __future__ import absolute_import, division, print_function
import scitbx.math as sm
import math
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from six.moves import range

r"""
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
  for ii in range(n_points):
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
  for ii in range(n_points):
    x = (1.0-2.0*h_gen.nth_given_base(5, ii))*high_limit
    y = (1.0-2.0*h_gen.nth_given_base(7, ii))*high_limit
    result += math.exp( -(x*x+y*y)*0.5 )/norm
  tmp = high_limit*high_limit*4.0*result/(float(n_points))
  assert approx_equal(tmp,1.0,eps=1e-2)

def test_cube():
  square = sm.square_halton_sampling(0.1, 0.8,  10.0, 80.0)
  start_values= (0.1, 10.0)
  assert approx_equal( start_values, next(square), eps=1e-4 )
  for ii in range(24):
    next(square)
  assert square.state()==25
  square.set_state(0)
  assert approx_equal( start_values, next(square), eps=1e-4 )

def tst_five_d_cube():
  hcube = sm.halton(5)
  result = flex.double( [0,0,0,0,0] )
  for ii in range(5000):
    vec = flex.double(hcube.nth_all(ii))
    result += vec
  result = result/5000.0
  for ii in result:
    assert approx_equal(ii, 0.5, 0.01)

def run():
  test_halton_sequence_1(1000,2)
  test_halton_sequence_2(1000,10)
  test_cube()
  tst_five_d_cube()

if (__name__ == "__main__"):
  run()
  print("OK")

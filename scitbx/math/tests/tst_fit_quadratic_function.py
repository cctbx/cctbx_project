from __future__ import absolute_import, division, print_function
from scitbx.math import fit_quadratic_function as fqf
from libtbx.test_utils import approx_equal
from six.moves import range

def test_fit():
  x1_obs=[]
  x2_obs=[]
  x1m=1.0
  x2m=3.0
  a=1
  b=2
  c=-0.5
  y_obs=[]

  for ii in range(10):
    for jj in range(10):
      x1_obs.append(ii-5)
      x2_obs.append(jj-5)
      y_obs.append(
        a*(( (ii-5.0)-x1m )**2.0) +
        b*(( (jj-5.0)-x2m )**2.0) +
        2.0*c*( (ii-5.0)-x1m )*( (jj-5.0)-x2m )
        )

  fit =fqf.fit_quadratic_function_2d_data( x1_obs, x2_obs, y_obs  )
  assert approx_equal( fit.x1m, x1m, eps=1e-5)
  assert approx_equal( fit.x2m, x2m, eps=1e-5)
  assert approx_equal( fit.a, a, eps=1e-5)
  assert approx_equal( fit.b, b, eps=1e-5)
  assert approx_equal( fit.c, c, eps=1e-5)
  print('OK')

if (__name__ == "__main__"):
  test_fit()

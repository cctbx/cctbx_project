from __future__ import absolute_import, division, print_function
from scitbx import math
from scitbx.array_family import flex
from six.moves import range

def tst(N=3):
  weights = flex.double(N)*0.0+1.0
  mvo = math.multivariate_moments( weights )
  for ii in range(100000):
    tmp = 0.1*(1.0-2.0*flex.random_double(N))+flex.double(range(N))*0.0+1.0
    mvo.update(tmp)
  vcv = mvo.vcv_upper_triangle_packed()
  mean = mvo.mean()
  var = mvo.variance()
  for m in mean:
    assert abs(m-1.0)<1e-3
  for v in var:
    assert abs(v-0.2*0.2/12)<1e-4
  for c in vcv:
    assert abs(c/0.0033)<1e-3


if __name__ == "__main__":
  tst()
  print("OK")

from __future__ import absolute_import, division, print_function
from scitbx.math import angle_derivative_wrt_vectors
from libtbx.test_utils import approx_equal
from scitbx import matrix
from random import random
from six.moves import range

def finite_difference_approx(u, v):

  delta = 1.e-6
  e1_shift = matrix.col((0.5*delta, 0, 0))
  e2_shift = matrix.col((0, 0.5*delta, 0))
  e3_shift = matrix.col((0, 0, 0.5*delta))

  fd_u1 = ((u + e1_shift).angle(v) - (u - e1_shift).angle(v)) / delta
  fd_u2 = ((u + e2_shift).angle(v) - (u - e2_shift).angle(v)) / delta
  fd_u3 = ((u + e3_shift).angle(v) - (u - e3_shift).angle(v)) / delta

  fd_v1 = ((u).angle(v + e1_shift) - (u).angle(v - e1_shift)) / delta
  fd_v2 = ((u).angle(v + e2_shift) - (u).angle(v - e2_shift)) / delta
  fd_v3 = ((u).angle(v + e3_shift) - (u).angle(v - e3_shift)) / delta

  return ((fd_u1, fd_u2, fd_u3), (fd_v1, fd_v2, fd_v3))

def exercise():

  # find acceptable trial vectors
  while True:
    u = matrix.col.random(3,-1,1)
    v = matrix.col.random(3,-1,1)
    tst = u.accute_angle(v)
    # skip vectors of zero length
    if tst is None: continue
    # skip vectors less than about 1 degrees from collinear
    if tst < 0.017: continue
    # otherwise accept the trial
    break

  dtheta_du, dtheta_dv = angle_derivative_wrt_vectors(u,v)
  fd_dtheta_du, fd_dtheta_dv = finite_difference_approx(u,v)

  try:
    assert approx_equal(dtheta_du, fd_dtheta_du)
    assert approx_equal(dtheta_dv, fd_dtheta_dv)
  except AssertionError:
    print("problem at angle", u.angle(v))

  return


if (__name__ == "__main__"):

  for i in range(10): exercise()
  print("OK")

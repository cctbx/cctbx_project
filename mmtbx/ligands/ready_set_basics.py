from __future__ import absolute_import, division, print_function
import math

from scitbx import matrix

def construct_xyz(ba, bv,
                  aa, av,
                  da, dv,
                  period=3,
                  ):
  assert ba is not None
  assert aa is not None
  assert da is not None
  rn = matrix.col(ba.xyz)
  rca = matrix.col(aa.xyz)
  rc = matrix.col(da.xyz)
  rcca = rc -rca

  e0 = (rn - rca).normalize()
  e1 = (rcca - (rcca.dot(e0))*e0).normalize()
  e2 = e0.cross(e1)

  pi = math.pi
  alpha = math.radians(av)
  phi = math.radians(dv)

  rh_list = []
  for n in range(0, period):
    rh = rn + bv * (math.sin(alpha)*(math.cos(phi + n*2*pi/period)*e1 +
                                     math.sin(phi + n*2*pi/period)*e2) -
                    math.cos(alpha)*e0)
    rh_list.append(rh)
  return rh_list

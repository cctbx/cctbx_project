from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("cctbx_math_ext")
from cctbx_math_ext import *

def basis_of_mirror_plane_with_normal(u):
  """ Primitive setting assumed """
  assert u != (0,0,0)
  basis = ()
  for t in ((1,0,0), (0,1,0), (0,0,1),
            (1,1,0), (1,0,1), (0,1,1),
            (-1,1,0), (-1,0,1), (0,-1,1),
            (1,1,1),
            (-1,1,1), (1,-1,1), (1,1,-1)):
    if len(basis) == 2: break
    if u[0]*t[0] + u[1]*t[1] + u[2]*t[2] == 0:
      basis += (t,)
  return basis

from __future__ import absolute_import, division, print_function
from libtbx.test_utils import Exception_expected
import scitbx.matrix as mat

def exercise_basis_of_mirror_plane_with_normal():
  from cctbx.math_module import basis_of_mirror_plane_with_normal
  for u in [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,1,1), (0,-1,1), (-1,1,1)]:
    u = mat.col(u)
    v,w = [ mat.col(x) for x in basis_of_mirror_plane_with_normal(u.elems) ]
    assert u.dot(v) == 0
    assert u.dot(w) == 0
    assert abs(v.cross(w)) != 0
  try:
    basis_of_mirror_plane_with_normal((0,0,0))
  except AssertionError:
    pass
  else:
    raise Exception_expected


def run():
  exercise_basis_of_mirror_plane_with_normal()
  print('OK')

if __name__ == '__main__':
  run()

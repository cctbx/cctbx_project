from libtbx.test_utils import Exception_expected
import cctbx.matrix as mat

def exercise_basis_of_mirror_plane_with_normal():
  for u in [(1,0,0), (0,1,0), (0,0,1), (1,1,0), (1,1,1), (0,-1,1), (-1,1,1)]:
    u = mat.col(u)
    v,w = [ mat.col(x) for x in mat.basis_of_mirror_plane_with_normal(u) ]
    assert u.dot(v) == 0
    assert u.dot(w) == 0
    assert v.cross(w) != (0,0,0)
  try:
    mat.basis_of_mirror_plane_with_normal((0,0,0))
  except AssertionError:
    pass
  else:
    raise Exception_expected


def run():
  exercise_basis_of_mirror_plane_with_normal()
  print 'OK'

if __name__ == '__main__':
  run()

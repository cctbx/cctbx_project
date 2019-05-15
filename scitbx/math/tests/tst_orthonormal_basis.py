from __future__ import absolute_import, division, print_function
from scitbx.math import orthonormal_basis
from scitbx import matrix
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from six.moves import range

def input_pairs():
  yield matrix.col((1, 0, 1)), matrix.col((1, 1, 1))
  for i in range(10):
    yield matrix.col(flex.random_double(3)), matrix.col(flex.random_double(3))

def exercise_orthonormal_basis(eps=1e-12):
  for v0, v1 in input_pairs():
    basis = [ matrix.col(e) for e in orthonormal_basis(v0, v1) ]
    for e in basis:
      assert approx_equal(e.length(), 1, eps)
    e0, e1, e2 = basis
    assert approx_equal(e0.cross(e1), e2, eps)

    basis_perm = orthonormal_basis(v0, 2, v1, 1)
    assert approx_equal(basis_perm[2],  e0, eps)
    assert approx_equal(basis_perm[1],  e1, eps)
    assert approx_equal(basis_perm[0], -e2, eps)

    basis_perm = orthonormal_basis(v0, 2, v1, 0)
    assert approx_equal(basis_perm[2],  e0, eps)
    assert approx_equal(basis_perm[0],  e1, eps)
    assert approx_equal(basis_perm[1],  e2, eps)

    left_handed_basis = orthonormal_basis(v0, v1, right_handed=False)
    assert approx_equal(left_handed_basis, (e0, e1, -e2), eps)

def run():
  exercise_orthonormal_basis()
  print('OK')

if __name__ == '__main__':
  run()

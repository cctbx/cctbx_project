from __future__ import absolute_import, division, print_function
from scitbx import linalg
import scitbx.linalg.eigensystem
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import range

def exercise_random_symmetric_matrix(n_trials):
  symm_mat = linalg.random_normal_matrix_generator(4,4)\
           .symmetric_matrix_with_eigenvalues
  for i in range(n_trials):
    lambdas = (1, 2, 3, 4)
    a = symm_mat(flex.double(lambdas)).matrix_packed_u_as_symmetric()
    eigen = linalg.eigensystem.real_symmetric(a)
    assert approx_equal(eigen.values(), sorted(lambdas, reverse=True),
                        eps=1e-12)

def run():
  exercise_random_symmetric_matrix(50)
  print('OK')

if __name__ == '__main__':
  run()

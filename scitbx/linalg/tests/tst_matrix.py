import scitbx.linalg
import scitbx.linalg.eigensystem
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.linalg import matrix_normality_ratio

def exercise_random_normal_matrix():
  for m, n in [ (3,5), (4,5), (5,5), (5,4), (5,3) ]:
    gen = scitbx.linalg.random_normal_matrix_generator(m, n)
    for i in xrange(10):
      assert matrix_normality_ratio(gen.normal_matrix()) < 10

  sigma = flex.double((1, 2, 3))
  for m, n in [ (3,5), (4,5), (5,5), (5,4), (5,3) ]:
    gen = scitbx.linalg.random_normal_matrix_generator(m, n)
    a = gen.matrix_with_singular_values(sigma)
    eig_u = scitbx.linalg.eigensystem.real_symmetric(
      a.matrix_multiply(a.matrix_transpose()))
    eig_v = scitbx.linalg.eigensystem.real_symmetric(
      a.matrix_transpose().matrix_multiply(a))
    assert approx_equal(list(eig_u.values()), [9, 4, 1] + [0,]*(m-3))
    assert approx_equal(list(eig_v.values()),  [9, 4, 1] + [0,]*(n-3))

def exercise_householder():
  a = flex.double([
    ( 4.0, 1.0,-1.0, 2.0),
    (1.0, 4.0, 1.0,-1.0),
    (-1.0, 1.0, 4.0, 1.0),
    (2.0,-1.0, 1.0, 4.0)])
  qr = scitbx.linalg.householder_qr_decomposition(a, may_accumulate_q=True)
  del a
  q = qr.q()

def run():
  exercise_householder()
  exercise_random_normal_matrix()
  print 'OK'

if __name__ == '__main__':
  run()

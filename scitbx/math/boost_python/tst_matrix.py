import scitbx.math
import scitbx.math.svd
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.math import matrix_equality_ratio, matrix_normality_ratio
from libtbx.test_utils import Exception_expected

def exercise_random_normal_matrix():
  for m, n in [ (3,5), (4,5), (5,5), (5,4), (5,3) ]:
    gen = scitbx.math.random_normal_matrix_generator(m, n)
    for i in xrange(10):
      assert matrix_normality_ratio(gen.normal_matrix()) < 10

  sigma = flex.double((1, 2, 3))
  for m, n in [ (3,5), (4,5), (5,5), (5,4), (5,3) ]:
    gen = scitbx.math.random_normal_matrix_generator(m, n)
    a = gen.matrix_with_singular_values(sigma)
    eig_u = scitbx.math.eigensystem.real_symmetric(
      a.matrix_multiply(a.matrix_transpose()))
    eig_v = scitbx.math.eigensystem.real_symmetric(
      a.matrix_transpose().matrix_multiply(a))
    assert approx_equal(list(eig_u.values()), [9, 4, 1] + [0,]*(m-3))
    assert approx_equal(list(eig_v.values()), [9, 4, 1] + [0,]*(n-3))



def exercise_svd():
  a = flex.double(xrange(1,19))
  sigma = [ 45.8945322027251, 1.6407053035305987, 0 ]
  a.resize(flex.grid(6,3))
  svd = scitbx.math.svd.real(
    a.deep_copy(),
    accumulate_u=True,
    accumulate_v=True)
  assert approx_equal(svd.sigma, sigma)
  a1 = svd.reconstruct()
  assert matrix_equality_ratio(a, a1) < 10
  assert matrix_normality_ratio(svd.u) < 10
  assert matrix_normality_ratio(svd.v) < 10
  svd = scitbx.math.svd.real(a.deep_copy(),
                             accumulate_u=False, accumulate_v=False)
  assert approx_equal(svd.sigma, sigma)
  assert not svd.u and not svd.v
  try:
    svd.reconstruct()
    raise Exception_expected
  except AssertionError:
    pass

  a = a.matrix_transpose()
  svd = scitbx.math.svd.real(
    a.deep_copy(),
    accumulate_u=True,
    accumulate_v=True)
  assert approx_equal(svd.sigma, sigma)
  a1 = svd.reconstruct()
  assert matrix_equality_ratio(a, a1) < 10
  assert matrix_normality_ratio(svd.u) < 10
  assert matrix_normality_ratio(svd.v) < 10

  a = flex.double(xrange(1,13))
  sigma = [25.436835633480246818, 1.7226122475210637387, 0]
  a.reshape(flex.grid(3,4))
  svd = scitbx.math.svd.real(
    a.deep_copy(),
    accumulate_u=True,
    accumulate_v=True)
  assert approx_equal(svd.sigma, sigma)
  a1 = svd.reconstruct()
  assert matrix_equality_ratio(a, a1) < 10
  assert matrix_normality_ratio(svd.u) < 10
  assert matrix_normality_ratio(svd.v) < 10

  a = a.matrix_transpose()
  svd = scitbx.math.svd.real(
    a.deep_copy(),
    accumulate_u=True,
    accumulate_v=True)
  assert approx_equal(svd.sigma, sigma)
  a1 = svd.reconstruct()
  assert matrix_equality_ratio(a, a1) < 10
  assert matrix_normality_ratio(svd.u) < 10
  assert matrix_normality_ratio(svd.v) < 10

def run():
  exercise_random_normal_matrix()
  exercise_svd()
  print 'OK'

if __name__ == '__main__':
  run()

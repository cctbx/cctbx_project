import scitbx.math
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.math import matrix_equality_ratio, matrix_normality_ratio
from libtbx.test_utils import Exception_expected

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
  exercise_svd()
  print 'OK'

if __name__ == '__main__':
  run()

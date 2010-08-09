def exercise():
  from scitbx.linalg import lapack_dgesvd_fem
  from scitbx.array_family import flex
  from scitbx import matrix
  from libtbx.test_utils import approx_equal
  #
  for diag in [0, 1]:
    for n in xrange(1, 11):
      a = flex.double(flex.grid(n,n), 0)
      for i in xrange(n):
        a[(i,i)] = diag
      a_inp = a.deep_copy()
      svd = lapack_dgesvd_fem(a=a)
      assert svd.info == 0
      assert approx_equal(svd.s, [diag]*n)
      assert svd.u.all() == (n,n)
      assert svd.vt.all() == (n,n)
  #
  def get_sigma(svd, m, n):
    elems = [0.] * (m*n)
    for i in xrange(min(m,n)):
      elems[i*n+i] = svd.s[i]
    return matrix.rec(elems=elems, n=(m,n))
  #
  mt = flex.mersenne_twister(seed=0)
  for m in xrange(1,11):
    for n in xrange(1,11):
      a = matrix.rec(elems=tuple(mt.random_double(m*n)*4-2), n=(m,n))
      svd = lapack_dgesvd_fem(a=a.transpose().as_flex_double_matrix())
      assert svd.info == 0
      sigma = get_sigma(svd, m, n)
      u = matrix.sqr(svd.u).transpose()
      v = matrix.sqr(svd.vt)
      assert approx_equal(u * sigma * v.transpose(), a)
  #
  a = matrix.rec(elems=[
     0.47,  0.10, -0.21,
    -0.21, -0.03, 0.35], n=(3,2))
  svd = lapack_dgesvd_fem(a=a.transpose().as_flex_double_matrix())
  assert svd.info == 0
  assert approx_equal(svd.s, [0.55981345199567534, 0.35931726783538481])
  assert approx_equal(svd.u, [
    0.81402078804155853, -0.5136261274467826, 0.27121644094748704,
    -0.42424674329757839, -0.20684171439391938, 0.88160717215094342,
    -0.39671760414380258, -0.83270925681813213, -0.38627766719264994])
  assert approx_equal(svd.vt, [
    0.8615633693608673, -0.50765003750177129,
    0.50765003750177129, 0.8615633693608673])

def compare_times():
  import scitbx.linalg.svd
  from scitbx.array_family import flex
  import time
  mt = flex.mersenne_twister(seed=0)
  samples = [100,200]
  print " m   n  real dgesvd"
  for m in samples:
    for n in samples:
      a = mt.random_double(size=m*n)*4-2
      a.reshape(flex.grid(m,n))
      ac = a.deep_copy()
      t0 = time.time()
      svd_real = scitbx.linalg.svd.real(
        ac, accumulate_u=True, accumulate_v=True)
      time_svd_real = time.time() - t0
      at = a.matrix_transpose()
      t0 = time.time()
      dgesvd = scitbx.linalg.lapack_dgesvd_fem(a=at)
      time_dgesvd = time.time() - t0
      print "%3d %3d %4.2f %4.2f" % (m, n, time_svd_real, time_dgesvd)

def run(args):
  assert len(args) == 0
  import scitbx.linalg
  lapack_dgesvd_fem = getattr(scitbx.linalg, "lapack_dgesvd_fem", None)
  if (lapack_dgesvd_fem is None):
    print "Skipping tests: lapack_dgesvd_fem not available."
  else:
    exercise()
    compare_times()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

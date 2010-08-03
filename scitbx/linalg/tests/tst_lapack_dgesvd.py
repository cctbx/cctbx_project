def exercise():
  import scitbx.linalg
  lapack_dgesvd = getattr(scitbx.linalg, "lapack_dgesvd_fem", None)
  if (lapack_dgesvd is None):
    print "Skipping tests: lapack_dgesvd_fem not available."
    return
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
      svd = lapack_dgesvd(a=a)
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
      svd = lapack_dgesvd(a=a.transpose().as_flex_double_matrix())
      assert svd.info == 0
      sigma = get_sigma(svd, m, n)
      u = matrix.sqr(svd.u).transpose()
      v = matrix.sqr(svd.vt)
      assert approx_equal(u * sigma * v.transpose(), a)
  #
  a = matrix.rec(elems=[
     0.47,  0.10, -0.21,
    -0.21, -0.03, 0.35], n=(3,2))
  svd = lapack_dgesvd(a=a.transpose().as_flex_double_matrix())
  assert svd.info == 0
  assert approx_equal(svd.s, [0.55981345199567534, 0.35931726783538481])
  assert approx_equal(svd.u, [
    0.81402078804155853, -0.5136261274467826, 0.27121644094748704,
    -0.42424674329757839, -0.20684171439391938, 0.88160717215094342,
    -0.39671760414380258, -0.83270925681813213, -0.38627766719264994])
  assert approx_equal(svd.vt, [
    0.8615633693608673, -0.50765003750177129,
    0.50765003750177129, 0.8615633693608673])

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

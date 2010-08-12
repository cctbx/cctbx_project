def exercise_impl(svd_impl_name, use_fortran):
  import scitbx.linalg
  svd_impl = getattr(scitbx.linalg, "lapack_%s" % svd_impl_name)
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
      svd = svd_impl(a=a, use_fortran=use_fortran)
      if (svd is None):
        if (not use_fortran):
          print "Skipping tests: lapack_%s not available." % svd_impl_name
        return
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
      svd = svd_impl(
        a=a.transpose().as_flex_double_matrix(), use_fortran=use_fortran)
      assert svd.info == 0
      sigma = get_sigma(svd, m, n)
      u = matrix.sqr(svd.u).transpose()
      v = matrix.sqr(svd.vt)
      assert approx_equal(u * sigma * v.transpose(), a)
  #
  a = matrix.rec(elems=[
     0.47,  0.10, -0.21,
    -0.21, -0.03, 0.35], n=(3,2))
  svd = svd_impl(
    a=a.transpose().as_flex_double_matrix(), use_fortran=use_fortran)
  assert svd.info == 0
  assert approx_equal(svd.s, [0.55981345199567534, 0.35931726783538481])
  assert approx_equal(svd.u, [
    0.81402078804155853, -0.5136261274467826, 0.27121644094748704,
    -0.42424674329757839, -0.20684171439391938, 0.88160717215094342,
    -0.39671760414380258, -0.83270925681813213, -0.38627766719264994])
  assert approx_equal(svd.vt, [
    0.8615633693608673, -0.50765003750177129,
    0.50765003750177129, 0.8615633693608673])

def exercise():
  for svd_impl_name in ["dgesvd", "dgesdd"]:
    for use_fortran in [False, True]:
      exercise_impl(svd_impl_name=svd_impl_name, use_fortran=use_fortran)

def compare_times(
      comprehensive=False,
      svd_impl_name="dgesvd",
      use_fortran=False):
  import scitbx.linalg.svd
  from scitbx.array_family import flex
  import time
  from libtbx.utils import progress_displayed_as_fraction
  mt = flex.mersenne_twister(seed=0)
  samples = []
  if not comprehensive:
    dims = (100, 200)
    for m in dims:
      for n in dims:
        samples.append((m, n))
  else:
    if comprehensive == "comprehensive-timing-1":
      dims = range(100, 600, 100)
      for m in dims:
        for n in dims:
          samples.append((m, n))
    elif comprehensive == "comprehensive-timing-2":
      for k in (1, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                2, 3, 4, 5, 6, 7, 8, 9, 10,
                20, 30, 40, 50, 60, 70, 80, 90, 100):
        for n in (10, 20, 30, 40, 50, 60, 70, 80, 90, 100):
          samples.append((int(k*n), n))
    else:
      raise RuntimeError(comprehensive)
  if not comprehensive:
    header = " m   n  real %s %s" % (
      svd_impl_name, ["fem", "for"][int(use_fortran)])
  else:
    handwritten_wrt_lapack = []
    progress = progress_displayed_as_fraction(len(samples))
  lapack_svd_impl = getattr(scitbx.linalg, "lapack_%s" % svd_impl_name)
  try:
    for m, n in samples:
      if comprehensive: progress.advance()
      a = mt.random_double(size=m*n)*4-2
      a.reshape(flex.grid(m,n))
      ac = a.deep_copy()
      t0 = time.time()
      svd_real = scitbx.linalg.svd.real(
        ac, accumulate_u=True, accumulate_v=True)
      time_svd_real = time.time() - t0
      at = a.matrix_transpose()
      t0 = time.time()
      svd_lapack = lapack_svd_impl(a=at, use_fortran=use_fortran)
      if (svd_lapack is None):
        return
      time_dgesvd = time.time() - t0
      if not comprehensive:
        if (header is not None):
          print header
          header = None
        print "%3d %3d %4.2f %4.2f" % (m, n, time_svd_real, time_dgesvd)
      else:
        handwritten_wrt_lapack.append((m, n, time_svd_real/time_dgesvd))
  finally:
    if comprehensive:
      print "handwrittenwrtlapack%s%s={" % (svd_impl_name,
                                            ['', 'fortran'][use_fortran])
      print ",".join([ "{%3d, %3d, %4.2f}" % (m, n, t)
                       for (m, n, t) in handwritten_wrt_lapack ])
      print "}"

def run(args):
  assert len(args) in (0, 1)
  exercise()
  if (len(args) == 0):
    for svd_impl_name in ["dgesvd", "dgesdd"]:
      for use_fortran in [False, True]:
        compare_times(svd_impl_name=svd_impl_name, use_fortran=use_fortran)
  else:
    compare_times(args[0][2:])
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

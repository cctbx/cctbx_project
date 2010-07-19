def exercise():
  import scitbx.linalg
  lapack_dsyev = getattr(scitbx.linalg, "lapack_dsyev", None)
  if (lapack_dsyev is None):
    print "Skipping tests: lapack_dsyev not available."
    return
  from scitbx.array_family import flex
  from scitbx.math.tests.tst_math import matrix_mul
  from libtbx.test_utils import approx_equal
  import random
  random.seed(0)
  #
  for diag in [0, 1]:
    for n in xrange(1, 11):
      for uplo in ["U", "L"]:
        a = flex.double(flex.grid(n,n), 0)
        for i in xrange(n):
          a[(i,i)] = diag
        w = flex.double(n, -1e100)
        a_inp = a.deep_copy()
        info = lapack_dsyev(jobz="V", uplo=uplo, a=a, w=w)
        assert info == 0
        assert approx_equal(w, [diag]*n)
        if (diag != 0):
          assert approx_equal(a, a_inp)
  #
  for i_trial in xrange(10):
    for n in xrange(1, 11):
      for uplo in ["U", "L"]:
        a = flex.double(flex.grid(n,n))
        for i in xrange(n):
          for j in xrange(i,n):
            a[i*n+j] = random.random() - 0.5
            if (i != j):
              a[j*n+i] = a[i*n+j]
        w = flex.double(n, -1e100)
        a_inp = a.deep_copy()
        info = lapack_dsyev(jobz="V", uplo=uplo, a=a, w=w)
        assert info == 0
        for i in xrange(1,n):
          assert w[i-1] <= w[i]
        for i in xrange(n):
          l = w[i]
          x = a[i*n:i*n+n]
          ax = matrix_mul(a_inp, n, n, x, n, 1)
          lx = [e*l for e in x]
          assert approx_equal(ax, lx)
  #
  a = flex.double([
     0.47,  0.10, -0.21,
     0.10,  0.01, -0.03,
    -0.21, -0.03, 0.35])
  a.reshape(flex.grid(3,3))
  w = flex.double(3, -1e100)
  info = lapack_dsyev(jobz="V", uplo=uplo, a=a, w=w)
  assert info == 0
  assert approx_equal(w, [-0.0114574, 0.1978572, 0.6436002])
  assert approx_equal(a, [
    -0.2236115,  0.9734398, -0.0491212,
    -0.5621211, -0.1699700, -0.8094010,
    -0.7962523, -0.1533793,  0.5851983])

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

import scitbx.linalg
from scitbx.linalg import eigensystem, time_eigensystem_real_symmetric
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from scitbx.math.tests.tst_math import matrix_mul
import random
import time

def exercise_eigensystem():
  s = eigensystem.real_symmetric(
    m=flex.double(flex.grid(0,0)),
    relative_epsilon=1e-6,
    absolute_epsilon=1e-6)
  s.values().all() == (0,)
  assert s.vectors().all() == (0,0)
  u = s.generalized_inverse_as_packed_u()
  assert u.all() == (0,)
  m = u.matrix_packed_u_as_symmetric()
  assert m.all() == (0,0)
  #
  for n in xrange(1,10):
    m = flex.double(flex.grid(n,n))
    s = eigensystem.real_symmetric(m)
    assert approx_equal(tuple(s.values()), [0]*n)
    v = s.vectors()
    for i in xrange(n):
      for j in xrange(n):
        x = 0
        if (i == j): x = 1
        assert approx_equal(v[(i,j)], x)
    v = []
    for i in xrange(n):
      j = (i*13+17) % n
      v.append(j)
      m[i*(n+1)] = j
    s = eigensystem.real_symmetric(m)
    if (n == 3):
      ss = eigensystem.real_symmetric((m[0],m[4],m[8],m[1],m[2],m[5]))
      assert approx_equal(s.values(), ss.values())
      assert approx_equal(s.vectors(), ss.vectors())
    v.sort()
    v.reverse()
    assert approx_equal(s.values(), v)
    if (n > 1):
      assert approx_equal(flex.min(s.vectors()), 0)
    assert approx_equal(flex.max(s.vectors()), 1)
    assert approx_equal(flex.sum(s.vectors()), n)
    for t in xrange(10):
      for i in xrange(n):
        for j in xrange(i,n):
          m[i*n+j] = random.random() - 0.5
          if (i != j):
            m[j*n+i] = m[i*n+j]
      s = eigensystem.real_symmetric(m)
      if (n == 3):
        ss = eigensystem.real_symmetric((m[0],m[4],m[8],m[1],m[2],m[5]))
        assert approx_equal(s.values(), ss.values())
        assert approx_equal(s.vectors(), ss.vectors())
      v = list(s.values())
      v.sort()
      v.reverse()
      assert list(s.values()) == v
      for i in xrange(n):
        l = s.values()[i]
        x = s.vectors()[i*n:i*n+n]
        mx = matrix_mul(m, n, n, x, n, 1)
        lx = [e*l for e in x]
        assert approx_equal(mx, lx)
  #
  m = (1.4573362052597449, 1.7361052947659894, 2.8065584999742659,
       -0.5387293498219814, -0.018204949672480729, 0.44956507395617257)
  n_repetitions = 10000
  t0 = time.time()
  v = time_eigensystem_real_symmetric(m, n_repetitions)
  assert v == (0,0,0)
  print "time_eigensystem_real_symmetric: %.3f micro seconds" % (
    (time.time() - t0)/n_repetitions*1.e6)
  from scitbx.linalg import time_lapack_dsyev
  for use_fortran in [False, True]:
    if (not use_fortran):
      if (not scitbx.linalg.fem_is_available()): continue
      impl_id = "fem"
    else:
      if (not scitbx.linalg.for_is_available()): continue
      impl_id = "for"
    v = time_lapack_dsyev(m, 2, use_fortran)
      # to trigger one-time initialization of SAVE variables
    t0 = time.time()
    v = time_lapack_dsyev(m, n_repetitions, use_fortran)
    assert v == (0,0,0)
    print "time_lapack_dsyev %s: %.3f micro seconds" % (
      impl_id, (time.time() - t0)/n_repetitions*1.e6)
  #
  s = eigensystem.real_symmetric(m=m)
  assert s.min_abs_pivot() > 0
  assert s.min_abs_pivot() < 1.e-10
  assert approx_equal(s.generalized_inverse_as_packed_u(), [
    0.77839538602575065, 0.25063185439711611, -0.03509803174624003,
    0.68162798233326816, -0.10755998636596431, 0.37330996497431423])
  s = eigensystem.real_symmetric(m=m, absolute_epsilon=10)
  assert s.min_abs_pivot() == 10
  assert approx_equal(s.vectors(), [0, 0, 1, 0, 1, 0, 1, 0, 0])
  assert approx_equal(s.values(),
    [2.8065584999742659, 1.7361052947659894, 1.4573362052597449])
  s = eigensystem.real_symmetric(m=m, relative_epsilon=0)
  assert s.min_abs_pivot() == 0
  assert approx_equal(s.values(), [3,2,1])

def compare_times(max_n_power=8):
  from scitbx.linalg import lapack_dsyev
  mt = flex.mersenne_twister(seed=0)
  show_tab_header = True
  tfmt = "%5.2f"
  for n_power in xrange(5,max_n_power+1):
    n = 2**n_power
    l = mt.random_double(size=n*(n+1)//2) * 2 - 1
    a = l.matrix_packed_l_as_symmetric()
    aes = a.deep_copy()
    ala = [a.deep_copy(), a.deep_copy()]
    wla = [flex.double(n, -1e100), flex.double(n, -1e100)]
    t0 = time.time()
    es = eigensystem.real_symmetric(aes)
    tab = [n, tfmt % (time.time() - t0)]
    for i,use_fortran in enumerate([False, True]):
      t0 = time.time()
      info = lapack_dsyev(
        jobz="V", uplo="U", a=ala[i], w=wla[i], use_fortran=use_fortran)
      if (info == 99):
        tab.append(" --- ")
      else:
        assert info == 0
        tab.append(tfmt % (time.time() - t0))
        assert approx_equal(list(reversed(es.values())), wla[i])
    if (show_tab_header):
      print "      time [s]           eigenvalues"
      print " n    es   fem   for     min    max"
      show_tab_header = False
    tab.extend([es.values()[-1], es.values()[0]])
    print "%3d %s %s %s [%6.2f %6.2f]" % tuple(tab)

def run():
  exercise_eigensystem()
  compare_times()
  print 'OK'

if __name__ == '__main__':
  run()

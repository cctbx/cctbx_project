from scitbx.math import eigensystem
from scitbx.array_family import flex
from scitbx.test_utils import approx_equal
import random

def matrix_mul(a, ar, ac, b, br, bc):
  assert br == ac
  result = []
  for i in xrange(ar):
    for k in xrange(bc):
      s = 0
      for j in xrange(ac):
        s += a[i * ac + j] * b[j * bc + k]
      result.append(s)
  return result

def run():
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
  print "OK"

if (__name__ == "__main__"):
  run()

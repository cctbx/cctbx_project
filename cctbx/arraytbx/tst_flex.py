import math
from cctbx_boost.arraytbx import flex

def exercise_flex_grid():
  g = flex.grid()
  assert g.nd() == 0
  assert g.size1d() == 0
  assert g.origin() == ()
  assert g.grid() == ()
  assert g.last() == ()
  assert g.last(1) == ()
  assert g.last(0) == ()
  g = flex.grid((2,3,5))
  assert g.nd() == 3
  assert g.size1d() == 30
  assert g.origin() == (0,0,0)
  assert g.grid() == (2,3,5)
  assert g.last() == (2,3,5)
  assert g.last(1) == (2,3,5)
  assert g.last(0) == (1,2,4)
  assert g((0,0,0)) == 0
  assert g((1,2,4)) == 29
  g = flex.grid((1,2,3), (4,6,8))
  assert g.nd() == 3
  assert g.size1d() == 60
  assert g.origin() == (1,2,3)
  assert g.grid() == (3,4,5)
  assert g.last() == (4,6,8)
  assert g.last(1) == (4,6,8)
  assert g.last(0) == (3,5,7)
  assert g((1,2,3)) == 0
  assert g((3,5,7)) == 59
  assert not g.is_valid_index((0,0,0))
  g = flex.grid((1,2,3), (4,6,8), 0)
  assert g.nd() == 3
  assert g.size1d() == 120
  assert g.origin() == (1,2,3)
  assert g.grid() == (4,5,6)
  assert g.last() == (5,7,9)
  assert g.last(1) == (5,7,9)
  assert g.last(0) == (4,6,8)
  assert g((1,2,3)) == 0
  assert g((4,6,8)) == 119
  assert not g.is_valid_index((0,0,0))
  assert not g.is_valid_index((5,0,0))

def exercise_flex_constructors():
  f = flex.double()
  assert f.size() == 0
  assert f.capacity() == 0
  assert f.accessor().nd() == 1
  assert tuple(f) == ()
  f = flex.double(flex.grid((2,3,5)))
  assert f.size() == 30
  assert f.capacity() == 30
  assert list(f) == [0] * 30
  f = flex.double(flex.grid((2,3,5)), 42)
  assert f.size() == 30
  assert list(f) == [42] * 30
  f = flex.double(2)
  assert f.size() == 2
  assert tuple(f) == (0,0)
  f = flex.double(1, 42)
  assert f.size() == 1
  assert tuple(f) == (42,)
  f = flex.double((1,2,3,4,5))
  assert f.size() == 5
  assert tuple(f) == (1,2,3,4,5)

def exercise_reinterpret():
  fc = flex.complex_double(flex.grid((2,3,5)))
  assert fc.accessor().grid() == (2,3,5)
  fcd = flex.reinterpret_complex_as_real(fc)
  assert fc.id() == fcd.id()
  assert fc.size() * 2 == fcd.size()
  assert fcd.accessor().grid() == (2,3,10)
  fcdc = flex.reinterpret_real_as_complex(fcd)
  assert fc.id() == fcdc.id()
  assert fc.size() == fcdc.size()
  assert fcdc.accessor().grid() == (2,3,5)

def exercise_misc():
  f = flex.double((1,2,3))
  assert f.front() == 1
  assert f.back() == 3
  f.fill(42)
  assert tuple(f) == (42,42,42)
  f = flex.double((1,2,3))
  fc = f.deep_copy()
  assert fc.id() != f.id()
  assert tuple(fc) == (1,2,3)
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert tuple(f1) == (1,2,3)
  f = flex.double(flex.grid((2,3)))
  assert f.accessor().nd() == 2
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert f1.accessor().nd() == 1
  assert f1.size() == 6
  i = flex.double_items()
  f = flex.double(flex.grid((2,3)), 12)
  i = flex.double_items(f)
  assert list(f.indices()) == range(6)
  assert list(f.items()) == zip(range(6), [12] * 6)

def exercise_push_back_etc():
  a = flex.double(3)
  assert a.size() == 3
  assert tuple(a) == (0, 0, 0)
  a = flex.double()
  assert a.size() == 0
  a.assign(3, 1)
  assert tuple(a) == (1, 1, 1)
  a.push_back(2)
  assert tuple(a) == (1, 1, 1, 2)
  a.insert(0, 3)
  assert tuple(a) == (3, 1, 1, 1, 2)
  a.insert(2, 3, 4)
  assert tuple(a) == (3, 1, 4, 4, 4, 1, 1, 2)
  a.pop_back()
  assert tuple(a) == (3, 1, 4, 4, 4, 1, 1)
  a.erase(2)
  assert tuple(a) == (3, 1, 4, 4, 1, 1)
  a.erase(3, 5)
  assert tuple(a) == (3, 1, 4, 1)
  a.resize(6)
  assert tuple(a) == (3, 1, 4, 1, 0, 0)
  a.resize(8, -1)
  assert tuple(a) == (3, 1, 4, 1, 0, 0, -1, -1)
  a.clear()
  assert a.size() == 0
  a = flex.double((0, 1, 2, 3))
  b = flex.double((4, 5, 6))
  a.append(b)
  assert tuple(a) == (0, 1, 2, 3, 4, 5, 6)
  assert tuple(a.indices()) == tuple(xrange(a.size()))
  assert list(a.items()) == zip(xrange(a.size()), a)
  g = flex.grid((2,3,5))
  f = flex.double(g, 11)
  g = flex.grid((2,3,7))
  f.resize(g)
  assert list(f) == [11] * 30 + [0] * 12
  g = flex.grid((2,3,9))
  f.resize(g, 22)
  assert list(f) == [11] * 30 + [0] * 12 + [22] * 12
  g = flex.grid((1,2,3))
  f.resize(g)
  assert list(f) == [11] * 6

def exercise_setitem():
  a = flex.double(2)
  a[0] = 11
  a[1] = 12
  assert tuple(a) == (11, 12)
  g = flex.grid((2,3))
  a = flex.double(g)
  for i in xrange(2):
    for j in xrange(3):
      a[(i,j)] = i * 3 + j
  assert list(a) == range(6)

def exercise_operators():
  a = flex.bool((0, 1, 0, 1))
  b = flex.bool((0, 1, 1, 0))
  assert tuple(~a) == (1, 0, 1, 0)
  assert tuple(a & b) == (0, 1, 0, 0)
  assert tuple(a | b) == (0, 1, 1, 1)
  a = flex.int((4, 9))
  b = flex.int((2, 3))
  assert tuple(-a) == (-4, -9)
  assert tuple(a + b) == (6, 12)
  assert tuple(a - b) == (2, 6)
  assert tuple(a * b) == (8, 27)
  assert tuple(a / b) == (2, 3)
  assert tuple(a % b) == (0, 0)
  assert tuple(a.add(3)) == (7, 12)
  assert tuple(a.sub(4)) == (0, 5)
  assert tuple(a.mul(5)) == (20, 45)
  assert tuple(a.div(2)) == (2, 4)
  assert tuple(a.mod(2)) == (0, 1)
  assert flex.sum(a) == 13
  a = flex.int((4, 9))
  b = flex.int((2, 12))
  assert tuple(a == b) == (0, 0)
  assert tuple(a != b) == (1, 1)
  assert tuple(a < b) == (0, 1)
  assert tuple(a > b) == (1, 0)
  assert tuple(a <= b) == (0, 1)
  assert tuple(a >= b) == (1, 0)
  assert tuple(a == 9) == (0, 1)
  assert tuple(a.as_double()) == (4, 9)
  assert cmp(a, b) > 0
  assert cmp(b, a) < 0
  assert cmp(a, a) == 0
  b = a.deep_copy()
  assert cmp(a, b) == 0
  assert a.cmp(4) > 0
  assert a.cmp(5) < 0
  a = flex.int((1, 1))
  assert a.cmp(1) == 0
  a = flex.miller_Index(((1,2,3), (2,3,4)))
  assert cmp(a, a) == 0
  b = a.deep_copy()
  assert cmp(a, b) == 0

def exercise_bool_inplace_operators():
  a = flex.bool((0, 1, 0, 1))
  b = flex.bool((0, 1, 1, 0))
  a &= b
  assert tuple(a) == (0, 1, 0, 0)
  a |= flex.bool((1, 0, 1, 0))
  assert tuple(a) == (1, 1, 1, 0)
  assert a.count(0) == 1
  assert a.count(1) == 3
  a &= 1 # XXX memory leak with Linux gcc 3.0.4, but not Tru64 cxx
  assert tuple(a) == (1, 1, 1, 0)
  a &= 0 # XXX memory leak
  assert tuple(a) == (0, 0, 0, 0)
  a |= 1 # XXX memory leak
  assert tuple(a) == (1, 1, 1, 1)
  a |= 0 # XXX memory leak
  assert tuple(a) == (1, 1, 1, 1)

def exercise_arith_inplace_operators():
  a = flex.int((4, 9))
  a += 3
  assert tuple(a) == (7, 12)
  a -= 3
  assert tuple(a) == (4, 9)
  a *= 3
  assert tuple(a) == (12, 27)
  a /= 3
  assert tuple(a) == (4, 9)
  a %= 3
  assert tuple(a) == (1, 0)

def exercise_functions():
  a = flex.int((-1, 0, 1))
  assert tuple(flex.abs(a)) == (1, 0, 1)
  a = flex.double((1, 0, 3, 2))
  b = flex.double((4, 5, 6))
  assert flex.min_index(a) == 1
  assert flex.min_index(b) == 0
  assert flex.max_index(a) == 2
  assert flex.max_index(b) == 2
  assert flex.min(a) == 0
  assert flex.min(b) == 4
  assert flex.max(a) == 3
  assert flex.max(b) == 6
  assert flex.sum(a) == 6
  assert flex.sum(b) == 15
  assert flex.product(a) == 0
  assert flex.product(b) == 120
  assert flex.mean(a) == 6. / 4
  assert flex.mean(b) == 15. / 3
  assert flex.mean_sq(a) == (1+3.*3.+2.*2.) / 4
  assert flex.mean_sq(b) == (4.*4.+5.*5.+6.*6.) / 3
  a = flex.double((-2, 0, 3))
  assert tuple(flex.pow(a, 2)) == (4, 0, 9)
  a = flex.double((2, 0, 3))
  assert list(flex.sqrt(a)) == [math.sqrt(x) for x in a]
  b = flex.double((1, 1, 1))
  assert (flex.mean(a) - flex.mean_weighted(a, b)) < 1.e-6
  assert (flex.mean_sq(a) - flex.mean_sq_weighted(a, b)) < 1.e-6

def exercise_complex_functions():
  c = 1+2j
  x = flex.complex_double((c,))
  a = flex.abs(x)
  assert abs(a[0] - abs(c)) < 1.e-6
  p = flex.arg(x)
  y = flex.polar(a, p)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  p = flex.arg(x, 0)
  y = flex.polar(a, p, 0)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  p = flex.arg(x, 1)
  y = flex.polar(a, p, 1)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  y = flex.polar(a, p, 0)
  d = y[0]
  assert abs(d.real - c.real) > 1.e-6
  assert abs(d.imag - c.imag) > 1.e-6

def exercise_select_shuffle():
  a = flex.double((1,2,3,4,5))
  b = flex.bool((0,1,0,1,1))
  assert tuple(a.select(b)) == (2,4,5)
  b = flex.size_t((3,1,0,4,2))
  assert tuple(a.shuffle(b)) == (3,2,5,1,4)

def exercise_exceptions():
  f = flex.double(flex.grid((2,3)))
  try: f.assign(1, 0)
  except: pass
  else: raise AssertionError, "No exception."
  try: f.push_back(0)
  except: pass
  else: raise AssertionError, "No exception."
  try: f[(2,0)]
  except: pass
  else: raise AssertionError, "No exception."

def run(iterations):
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_flex_grid()
    exercise_flex_constructors()
    exercise_reinterpret()
    exercise_misc()
    exercise_push_back_etc()
    exercise_setitem()
    exercise_operators()
    exercise_bool_inplace_operators()
    exercise_arith_inplace_operators()
    exercise_functions()
    exercise_complex_functions()
    exercise_select_shuffle()
    exercise_exceptions()
    i += 1

if (__name__ == "__main__"):
  import sys
  from cctbx.development import debug_utils
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "large",
  ))
  n = 1
  if (len(sys.argv) > 1 + Flags.n):
    n = int(Flags.regular_args[0])
  if (Flags.large):
    pickle_large_arrays(n)
  else:
    run(n)
  print "OK"

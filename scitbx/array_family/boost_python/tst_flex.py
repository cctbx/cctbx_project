import math
from scitbx.array_family import flex
from scitbx.test_utils import approx_equal

def exercise_flex_grid():
  g = flex.grid()
  assert g.nd() == 0
  assert g.size_1d() == 0
  assert g.origin() == ()
  assert g.grid() == ()
  assert g.last() == ()
  assert g.last(1) == ()
  assert g.last(0) == ()
  assert g.is_0_based()
  assert not g.is_padded()
  g = flex.grid((2,3,5))
  assert g.nd() == 3
  assert g.size_1d() == 30
  assert g.origin() == (0,0,0)
  assert g.grid() == (2,3,5)
  assert g.last() == (2,3,5)
  assert g.last(1) == (2,3,5)
  assert g.last(0) == (1,2,4)
  assert g((0,0,0)) == 0
  assert g((1,2,4)) == 29
  assert g.is_0_based()
  assert not g.is_padded()
  g = flex.grid((1,2,3), (4,6,8))
  assert g.nd() == 3
  assert g.size_1d() == 60
  assert g.origin() == (1,2,3)
  assert g.grid() == (3,4,5)
  assert g.last() == (4,6,8)
  assert g.last(1) == (4,6,8)
  assert g.last(0) == (3,5,7)
  assert g((1,2,3)) == 0
  assert g((3,5,7)) == 59
  assert not g.is_valid_index((0,0,0))
  assert not g.is_0_based()
  assert not g.is_padded()
  g = flex.grid((1,2,3), (4,6,8), 0)
  assert g.nd() == 3
  assert g.size_1d() == 120
  assert g.origin() == (1,2,3)
  assert g.grid() == (4,5,6)
  assert g.last() == (5,7,9)
  assert g.last(1) == (5,7,9)
  assert g.last(0) == (4,6,8)
  assert g((1,2,3)) == 0
  assert g((4,6,8)) == 119
  assert not g.is_valid_index((0,0,0))
  assert not g.is_valid_index((5,0,0))
  assert not g.is_0_based()
  assert not g.is_padded()
  g.set_layout((3,-9,5))
  assert not g.is_0_based()
  assert g.is_padded()
  import pickle
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert g.origin() == l.origin()
  assert g.grid() == l.grid()
  assert g.layout() == l.layout()
  assert g == l
  assert not g != l
  l = flex.grid((1,2,4), (4,6,8), 0).set_layout((3,-9,5))
  assert not g == l
  assert g != l
  l = flex.grid((1,2,3), (4,7,8), 0).set_layout((3,-9,5))
  assert not g == l
  assert g != l
  l = flex.grid((1,2,3), (4,6,8), 0).set_layout((4,-9,5))
  assert not g == l
  assert g != l
  g = flex.grid((1,2,3))
  assert g.shift_origin() == g
  g = flex.grid((1,2,3), (4,6,8))
  s = g.shift_origin()
  assert s.origin() == (0,0,0)
  assert s.grid() == g.grid()
  assert s.layout() == ()
  g = flex.grid((1,2,3), (4,6,8)).set_layout((3,5,7))
  s = g.shift_origin()
  assert s.origin() == (0,0,0)
  assert s.grid() == g.grid()
  assert s.layout() == (2,3,4)

def exercise_flex_constructors():
  f = flex.double()
  assert f.size() == 0
  assert f.capacity() == 0
  assert f.accessor().nd() == 1
  assert tuple(f.accessor().origin()) == (0,)
  assert tuple(f.accessor().grid()) == (0,)
  assert tuple(f.accessor().last()) == (0,)
  assert tuple(f.accessor().last(1)) == (0,)
  assert tuple(f.accessor().last(0)) == (-1,)
  assert tuple(f.accessor().layout()) == ()
  assert f.accessor().is_0_based()
  assert not f.accessor().is_padded()
  assert f.nd() == 1
  assert tuple(f.origin()) == (0,)
  assert tuple(f.grid()) == (0,)
  assert tuple(f.last()) == (0,)
  assert tuple(f.last(1)) == (0,)
  assert tuple(f.last(0)) == (-1,)
  assert tuple(f.layout()) == ()
  assert f.is_0_based()
  assert not f.is_padded()
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
  f = flex.double([2,1,3,5])
  assert f.size() == 4
  assert tuple(f) == (2,1,3,5)
  f = flex.double(xrange(10,13))
  assert f.size() == 3
  assert tuple(f) == (10,11,12)

def exercise_misc():
  f = flex.double((1,2,3))
  assert f[0] == 1
  assert f[2] == 3
  assert f[-1] == 3
  assert f[-3] == 1
  assert f.front() == 1
  assert f.back() == 3
  f.fill(42)
  assert tuple(f) == (42,42,42)
  f = flex.double((1,2,3,4,5,6))
  fc = f.deep_copy()
  assert fc.id() != f.id()
  assert tuple(fc) == (1,2,3,4,5,6)
  fc = f.shallow_copy()
  assert fc.id() == f.id()
  assert tuple(fc) == (1,2,3,4,5,6)
  fc.resize(flex.grid((2,3)))
  assert f.nd() == 1
  assert fc.nd() == 2
  f = flex.double((1,2,3))
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert tuple(f1) == (1,2,3)
  f = flex.double(flex.grid((2,3)))
  assert f.nd() == 2
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert f1.nd() == 1
  assert f1.size() == 6
  i = flex.double_items()
  f = flex.double(flex.grid((2,3)), 12)
  i = flex.double_items(f)
  assert list(f.indices()) == range(6)
  assert list(f.items()) == zip(range(6), [12] * 6)
  g = flex.grid((1,2,3), (4,6,8)).set_layout((3,5,7))
  f = flex.double(g)
  assert f.shift_origin().accessor() == g.shift_origin()

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

def exercise_select_shuffle():
  a = flex.double((1,2,3,4,5))
  b = flex.bool((0,1,0,1,1))
  assert tuple(a.select(b)) == (2,4,5)
  b = flex.size_t((3,1,0,4,2))
  assert tuple(a.shuffle(b)) == (3,2,5,1,4)

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
  assert tuple(a + 3) == (7, 12)
  assert tuple(a - 4) == (0, 5)
  assert tuple(a * 5) == (20, 45)
  assert tuple(a / 2) == (2, 4)
  assert tuple(a % 2) == (0, 1)
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
  assert a.all_eq(a)
  assert not a.all_eq(b)
  assert not a.all_eq(4)
  assert a.all_ne(b)
  assert not a.all_ne(a)
  assert not a.all_ne(4)
  assert a.all_ne(5)
  assert not a.all_lt(b)
  assert a.all_lt(10)
  assert not a.all_gt(b)
  assert a.all_gt(3)
  assert not a.all_le(b)
  assert a.all_le(9)
  assert not a.all_ge(b)
  assert a.all_ge(2)
  assert flex.order(a, b) == cmp(tuple(a), tuple(b))
  assert flex.order(b, a) == cmp(tuple(b), tuple(a))
  assert flex.order(a, a) == 0
  b = a.deep_copy()
  assert flex.order(a, b) == 0

def exercise_bool_inplace_operators():
  a = flex.bool((0, 1, 0, 1))
  b = flex.bool((0, 1, 1, 0))
  a &= b
  assert tuple(a) == (0, 1, 0, 0)
  a |= flex.bool((1, 0, 1, 0))
  assert tuple(a) == (1, 1, 1, 0)
  assert a.count(0) == 1
  assert a.count(1) == 3
  a &= 1
  assert tuple(a) == (1, 1, 1, 0)
  a &= 0
  assert tuple(a) == (0, 0, 0, 0)
  a |= 1
  assert tuple(a) == (1, 1, 1, 1)
  a |= 0
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
  assert tuple(flex.pow2(a)) == (1, 0, 1)
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

def exercise_linear_regression():
  x = flex.double((1,2,3))
  r = flex.linear_regression(x, x, 1.e-6)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 0)
  assert approx_equal(r.slope(), 1)
  assert approx_equal(r.cc(), 1)
  y = flex.double((-1./2+1,-2./2+1,-3./2+1))
  r = flex.linear_regression(x, y)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 1)
  assert approx_equal(r.slope(), -1./2)
  assert approx_equal(r.cc(), -1)
  y = flex.double((0,0,0))
  r = flex.linear_regression(x, y)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 0)
  assert approx_equal(r.slope(), 0)
  assert approx_equal(r.cc(), 1)
  r = flex.linear_regression(y, y)
  assert not r.is_well_defined()

def exercise_exceptions():
  f = flex.double(flex.grid((2,3)))
  try: f.assign(1, 0)
  except RuntimeError, e:
    assert str(e) == "Array must be 0-based 1-dimensional."
  else:
    raise AssertionError, "No exception or wrong exception."
  try: f.push_back(0)
  except RuntimeError, e:
    assert str(e) == "Array must be 0-based 1-dimensional."
  else:
    raise AssertionError, "No exception or wrong exception."
  try: f[(2,0)]
  except IndexError, e:
    assert str(e) == "Index out of range."
  else:
    raise AssertionError, "No exception or wrong exception."

def exercise_pickle_single_buffered():
  import pickle
  a = flex.bool((1,0,1))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,0,1)
  a = flex.double(())
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 0
  a = flex.double((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.int((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.float((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.complex_double((1+2j, 2+3j, 4+5j))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1+2j, 2+3j, 4+5j)
  a = flex.double(flex.grid((-1,2,-3), (7,5,3)).set_layout((3,-5,7)), 13)
  a[(1,2,-1)] = -8
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 8 * 3 * 6
  assert tuple(a) == tuple(b)
  assert a.origin() == b.origin()
  assert a.grid() == b.grid()
  assert a.layout() == b.layout()
  assert a.accessor() == b.accessor()

def exercise_pickle_double_buffered():
  import pickle
  a = flex.std_string()
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 0
  assert tuple(b) == ()
  a = flex.std_string(("abc", "bcd", "cde"))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (("abc", "bcd", "cde"))

def pickle_large_arrays(max_exp):
  import os, time
  import pickle, cPickle
  try:
    for array_type in (
        flex.bool,
        flex.size_t,
        flex.int,
        flex.long,
        flex.float,
        flex.double,
        flex.complex_double,
        flex.std_string):
      for e in xrange(max_exp+1):
        n = 2**e
        if (array_type == flex.bool):
          val = 1
        elif (array_type == flex.size_t):
          val = 2147483647
        elif (array_type == flex.int):
          val = -2147483647
        elif (array_type == flex.long):
          val = -9223372036854775808
          if (type(val) == type(1L)):
            val = -2147483647
        elif (array_type == flex.complex_double):
          val = complex(-1.234567809123456e+20, -1.234567809123456e+20)
        elif (array_type in (flex.float, flex.double)):
          val = -1.234567809123456e+20
        elif (array_type == flex.std_string):
          val = "x" * 10
        else:
          raise AssertionError, "Unexpected array type."
        a = array_type(n, val)
        for pickler in (0, pickle, cPickle):
          if (pickler == 0):
            pickler_name = "g/setstate"
            t0 = time.time()
            s = a.__getstate__()
            td = time.time() - t0
            t0 = time.time()
            b = array_type()
            b.__setstate__(s)
            tl = time.time() - t0
          else:
            pickler_name = pickler.__name__
            f = open("pickle.tmp", "wb")
            t0 = time.time()
            pickler.dump(a, f, 1)
            td = time.time() - t0
            f.close()
            f = open("pickle.tmp", "rb")
            t0 = time.time()
            b = pickle.load(f)
            tl = time.time() - t0
            f.close()
          print array_type.__name__, n, pickler_name, "%.2f %.2f" % (td, tl)
          sys.stdout.flush()
  finally:
    try: os.unlink("pickle.tmp")
    except: pass

def run(iterations):
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_flex_grid()
    exercise_flex_constructors()
    exercise_misc()
    exercise_push_back_etc()
    exercise_setitem()
    exercise_select_shuffle()
    exercise_operators()
    exercise_bool_inplace_operators()
    exercise_arith_inplace_operators()
    exercise_functions()
    exercise_complex_functions()
    exercise_linear_regression()
    exercise_exceptions()
    exercise_pickle_single_buffered()
    exercise_pickle_double_buffered()
    i += 1

if (__name__ == "__main__"):
  import sys
  from scitbx.python_utils import command_line
  Flags = command_line.parse_options(sys.argv[1:], (
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

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
  import pickle
  g.set_layout((3,-9,5))
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
  assert f.nd() == 1
  assert tuple(f.origin()) == (0,)
  assert tuple(f.grid()) == (0,)
  assert tuple(f.last()) == (0,)
  assert tuple(f.last(1)) == (0,)
  assert tuple(f.last(0)) == (-1,)
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
  assert fc.grid() == (2,3,5)
  fcd = flex.reinterpret_complex_as_real(fc)
  assert fc.id() == fcd.id()
  assert fc.size() * 2 == fcd.size()
  assert fcd.grid() == (2,3,10)
  fcdc = flex.reinterpret_real_as_complex(fcd)
  assert fc.id() == fcdc.id()
  assert fc.size() == fcdc.size()
  assert fcdc.grid() == (2,3,5)

def exercise_misc():
  f = flex.double((1,2,3))
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

def exercise_type_1_picklers():
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
  a = flex.miller_Index(((1,2,3), (-2,3,-4), (3,-4,5)))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == ((1,2,3), (-2,3,-4), (3,-4,5))
  a = flex.hendrickson_lattman(((1,2,3,4), (-2,3,-4,5), (3,-4,5,-6)))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == ((1,2,3,4), (-2,3,-4,5), (3,-4,5,-6))
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

def exercise_type_2_picklers():
  import pickle
  from cctbx_boost import sftbx
  from cctbx_boost.eltbx.caasf_wk1995 import CAASF_WK1995
  sf = CAASF_WK1995("C")
  sites = (
    sftbx.XrayScatterer("C1", sf, 1+2j, (0.1,0.2,0.3), 1, 0.1),
    sftbx.XrayScatterer("C2", sf, 2+3j, (0.2,0.3,0.4), 0.5, (1,2,3,4,5,6)),
    sftbx.XrayScatterer("C3", sf, 3+4j, (0.3,0.4,0.5), 1, (6,5,4,3,2,1)),
    sftbx.XrayScatterer("C4", sf, 4+5j, (0.4,0.5,0.6), 0.5, 0.2),
  )
  for nd in (1,2,3):
    if (nd == 1):
      a = flex.XrayScatterer(sites)
    elif (nd == 2):
      a.resize(flex.grid((1,2), (3,4)))
    else:
      a.resize(flex.grid((5,1,2), (7,2,4)))
    p = pickle.dumps(a)
    b = pickle.loads(p)
    assert b.size() == len(sites)
    assert a.origin() == b.origin()
    assert a.grid() == b.grid()
    for i in xrange(len(sites)):
      assert sites[i].Label() == b[i].Label()
      assert sites[i].CAASF().Label() == b[i].CAASF().Label()
      assert sites[i].fpfdp() == b[i].fpfdp()
      assert sites[i].Coordinates() == b[i].Coordinates()
      assert sites[i].Occ() == b[i].Occ()
      assert sites[i].isAnisotropic() == b[i].isAnisotropic()
      if (b[i].isAnisotropic()):
        assert sites[i].Uaniso() == b[i].Uaniso()
      else:
        assert sites[i].Uiso() == b[i].Uiso()

def pickle_large_arrays(max_exp):
  import time
  import pickle, cPickle
  for array_type in (
      flex.bool,
      flex.int,
      flex.long,
      flex.float,
      flex.double,
      flex.complex_double):
    for e in xrange(max_exp+1):
      n = 2**e
      if (array_type == flex.bool):
        val = 1
      elif (array_type == flex.int):
        val = -2147483647
      elif (array_type == flex.long):
        val = -9223372036854775808
        if (type(val) == type(1L)):
          val = -2147483647
      elif (array_type == flex.complex_double):
        val = complex(-1.234567809123456e+20, -1.234567809123456e+20)
      else:
        val = -1.234567809123456e+20
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
    exercise_type_1_picklers()
    exercise_type_2_picklers()
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

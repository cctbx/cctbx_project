import math
from cctbx_boost.arraytbx import shared

def exercise_basic():
  a = shared.double()
  b = shared.reinterpret_real_as_complex(a)
  assert a.id() == b.id()
  a = shared.reinterpret_complex_as_real(b)
  assert a.id() == b.id()
  b = a.deep_copy()
  assert a.id() != b.id()

def exercise_push_back_etc():
  a = shared.int(3)
  assert a.size() == 3
  assert tuple(a) == (0, 0, 0)
  a = shared.int()
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
  a.clear()
  assert a.size() == 0
  a = shared.double((0, 1, 2, 3))
  b = shared.double((4, 5, 6))
  a.append(b)
  assert tuple(a) == (0, 1, 2, 3, 4, 5, 6)
  assert tuple(a.indices()) == tuple(xrange(a.size()))
  assert list(a.items()) == zip(xrange(a.size()), a)

def exercise_operators():
  a = shared.bool((0, 1, 0, 1))
  b = shared.bool((0, 1, 1, 0))
  assert tuple(~a) == (1, 0, 1, 0)
  assert tuple(a & b) == (0, 1, 0, 0)
  assert tuple(a | b) == (0, 1, 1, 1)
  a = shared.int((4, 9))
  b = shared.int((2, 3))
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
  assert shared.sum(a) == 13
  a = shared.int((4, 9))
  assert tuple(a.as_double()) == (4, 9)
  b = shared.int((2, 12))
  assert tuple(a == b) == (0, 0)
  assert tuple(a != b) == (1, 1)
  assert tuple(a < b) == (0, 1)
  assert tuple(a > b) == (1, 0)
  assert tuple(a <= b) == (0, 1)
  assert tuple(a >= b) == (1, 0)
  assert tuple(a == 9) == (0, 1)

def exercise_bool_inplace_operators():
  a = shared.bool((0, 1, 0, 1))
  b = shared.bool((0, 1, 1, 0))
  a &= b
  assert tuple(a) == (0, 1, 0, 0)
  a |= shared.bool((1, 0, 1, 0))
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
  a = shared.int((4, 9))
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
  a = shared.int((-1, 0, 1))
  assert tuple(shared.abs(a)) == (1, 0, 1)
  a = shared.double((-2, 0, 3))
  assert tuple(shared.pow(a, 2)) == (4, 0, 9)
  a = shared.double((1, 0, 3, 2))
  b = shared.double((4, 5, 6))
  assert shared.min_index(a) == 1
  assert shared.min_index(b) == 0
  assert shared.max_index(a) == 2
  assert shared.max_index(b) == 2
  assert shared.min(a) == 0
  assert shared.min(b) == 4
  assert shared.max(a) == 3
  assert shared.max(b) == 6
  assert shared.sum(a) == 6
  assert shared.sum(b) == 15
  assert shared.product(a) == 0
  assert shared.product(b) == 120
  assert shared.mean(a) == 6. / 4
  assert shared.mean(b) == 15. / 3
  assert shared.rms(a) == math.sqrt((1+3.*3.+2.*2.) / 4)
  assert shared.rms(b) == math.sqrt((4.*4.+5.*5.+6.*6.) / 3)

def exercise_complex_functions():
  c = 1+2j
  x = shared.complex_double((c,))
  a = shared.abs(x)
  assert abs(a[0] - abs(c)) < 1.e-6
  p = shared.arg(x)
  y = shared.polar(a, p)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  p = shared.arg(x, 0)
  y = shared.polar(a, p, 0)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  p = shared.arg(x, 1)
  y = shared.polar(a, p, 1)
  d = y[0]
  assert abs(d.real - c.real) < 1.e-6
  assert abs(d.imag - c.imag) < 1.e-6
  y = shared.polar(a, p, 0)
  d = y[0]
  assert abs(d.real - c.real) > 1.e-6
  assert abs(d.imag - c.imag) > 1.e-6

def exercise_regression_and_statistics():
  x = shared.double((0, 1, 2, 3))
  y = shared.double((1, 3, 5, 7))
  r = shared.linear_regression(x, y)
  assert r.is_well_defined()
  assert abs(r.b() - 1) < 1.e-6
  assert abs(r.m() - 2) < 1.e-6
  assert abs(r.cc() - 1) < 1.e-6
  s = shared.statistics(x)
  assert (s.min() - 0) < 1.e-6
  assert (s.max() - 3) < 1.e-6
  assert (s.mean() - 6./4.) < 1.e-6
  assert (s.mean2() - 14./4.) < 1.e-6
  assert (s.sigma() - math.sqrt(14./4. - 36./16.)) < 1.e-6

def exercise_type_1_picklers():
  import pickle
  a = shared.bool((1,0,1))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,0,1)
  a = shared.double(())
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 0
  a = shared.double((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = shared.int((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = shared.float((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = shared.complex_double((1+2j, 2+3j, 4+5j))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1+2j, 2+3j, 4+5j)
  a = shared.miller_Index(((1,2,3), (-2,3,-4), (3,-4,5)))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == ((1,2,3), (-2,3,-4), (3,-4,5))
  a = shared.hendrickson_lattman(((1,2,3,4), (-2,3,-4,5), (3,-4,5,-6)))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == ((1,2,3,4), (-2,3,-4,5), (3,-4,5,-6))

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
  a = shared.XrayScatterer(sites)
  p = pickle.dumps(a)
  b = pickle.loads(p)
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
      shared.bool,
      shared.int,
      shared.long,
      shared.float,
      shared.double,
      shared.complex_double):
    for e in xrange(max_exp+1):
      n = 2**e
      if (array_type == shared.bool):
        val = 1
      elif (array_type == shared.int):
        val = -2147483647
      elif (array_type == shared.long):
        val = -9223372036854775808
        if (type(val) == type(1L)):
          val = -2147483647
      elif (array_type == shared.complex_double):
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
    exercise_basic()
    exercise_push_back_etc()
    exercise_operators()
    exercise_bool_inplace_operators()
    exercise_arith_inplace_operators()
    exercise_functions()
    exercise_complex_functions()
    exercise_regression_and_statistics()
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

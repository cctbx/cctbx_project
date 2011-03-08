from boost import rational
from libtbx.test_utils import Exception_expected, approx_equal, show_diff
import pickle

def exercise_int():
  ri = rational.int
  r = ri()
  assert r.numerator() == 0
  assert r.denominator() == 1
  assert r.as_tuple() == (0,1)
  assert int(r) == 0
  assert float(r) == 0
  assert rational.int(rational.int(3)).as_tuple() == (3,1)
  assert rational.int(2).as_tuple() == (2,1)
  assert rational.int(2,3).as_tuple() == (2,3)
  assert str(rational.int()) == "0"
  assert str(rational.int(2)) == "2"
  assert str(rational.int(-2,3)) == "-2/3"
  assert (-rational.int(2,3)).as_tuple() == (-2,3)
  assert (rational.int(2,3) + rational.int(3,4)).as_tuple() == (17,12)
  assert (rational.int(2,3) - rational.int(3,4)).as_tuple() == (-1,12)
  assert (rational.int(2,3) * rational.int(3,4)).as_tuple() == (1,2)
  assert (rational.int(2,3) / rational.int(3,4)).as_tuple() == (8,9)
  assert (rational.int(2,3) // rational.int(3,4)) == 0
  assert (rational.int(2,3) % rational.int(1,2)).as_tuple() == (1,6)
  assert (rational.int(2,3) + 4).as_tuple() == (14,3)
  assert (rational.int(2,3) - 4).as_tuple() == (-10,3)
  assert (rational.int(2,3) * 4).as_tuple() == (8,3)
  assert (rational.int(2,3) / 4).as_tuple() == (1,6)
  assert (rational.int(2,3) // 4) == 0
  assert (rational.int(7,3) % 2).as_tuple() == (1,3)
  assert (5 + rational.int(2,3)).as_tuple() == (17,3)
  assert (5 - rational.int(2,3)).as_tuple() == (13,3)
  assert (5 * rational.int(2,3)).as_tuple() == (10,3)
  assert (5 / rational.int(2,3)).as_tuple() == (15,2)
  assert (5 // rational.int(2,3)) == 7
  assert (5 % rational.int(2,3)).as_tuple() == (1,3)
  assert rational.int(2,3) == rational.int(2,3)
  assert not rational.int(2,3) == rational.int(2,5)
  assert rational.int(2,3) != rational.int(2,5)
  assert not rational.int(2,3) != rational.int(2,3)
  assert rational.int(2,3) < rational.int(3,4)
  assert not rational.int(2,3) < rational.int(2,3)
  assert rational.int(2,3) > rational.int(1,2)
  assert not rational.int(2,3) > rational.int(2,3)
  assert rational.int(2,3) <= rational.int(3,4)
  assert not rational.int(2,3) <= rational.int(1,2)
  assert rational.int(2,3) >= rational.int(1,2)
  assert not rational.int(2,3) >= rational.int(3,4)
  assert rational.int(4,2) == 2
  assert not rational.int(4,2) == 3
  assert rational.int(4,2) != 3
  assert not rational.int(4,2) != 2
  assert rational.int(4,2) < 3
  assert not rational.int(4,2) < 2
  assert rational.int(4,2) > 1
  assert not rational.int(4,2) > 2
  assert rational.int(4,2) <= 3
  assert not rational.int(4,2) <= 1
  assert rational.int(4,2) >= 1
  assert not rational.int(4,2) >= 3
  assert 2 == rational.int(4,2)
  assert not 3 == rational.int(4,2)
  assert 3 != rational.int(4,2)
  assert not 2 != rational.int(4,2)
  assert 3 > rational.int(4,2)
  assert not 2 > rational.int(4,2)
  assert 1 < rational.int(4,2)
  assert not 2 < rational.int(4,2)
  assert 3 >= rational.int(4,2)
  assert not 1 >= rational.int(4,2)
  assert 1 <= rational.int(4,2)
  assert not 3 <= rational.int(4,2)
  r = rational.int(4,3)
  r += 1
  assert r.as_tuple() == (7,3)
  assert approx_equal(float(r), 7./3)
  s = rational.int(4,3)
  assert hash(s) == hash(rational.int(4,3))
  assert hash(s) != hash(r)
  for n in xrange(-100,100):
    assert hash(n) == hash(rational.int(n))
    for d in xrange(1,8):
      assert hash(rational.int(n,d)) == hash(rational.int(n,d))
      assert hash(rational.int(n,d)) == hash(rational.int(3*n,3*d))
      assert hash(rational.int(n,d)) == hash(rational.int(-3*n,-3*d))
  try: int(r)
  except RuntimeError, e:
    assert str(e) == "boost.rational: as_int() conversion error:" \
      " denominator is different from one."
  else: raise Exception_expected
  for n in xrange(-5,6):
    for d in xrange(1,10):
      r = rational.int(n, d)
      p = pickle.dumps(r)
      l = pickle.loads(p)
      assert l == r
      assert str(l) == str(r)
  #
  ee = "bad rational: zero denominator"
  lhs = ri(1)
  for rhs in [ri(0), 0]:
    try: lhs / rhs
    except RuntimeError, e: assert not show_diff(str(e), ee)
    else: raise Exception_expected
    try: lhs % rhs
    except RuntimeError, e: assert not show_diff(str(e), ee)
    else: raise Exception_expected
  #
  try:
    import fractions
  except ImportError:
    fractions = None
  def check(nd1, nd2, expected=None):
    r1, r2 = ri(*nd1), ri(*nd2)
    rm = r1 % r2
    assert (r1 // r2) * r2 + rm == r1
    if (fractions is not None):
      ff = fractions.Fraction
      f1, f2 = ff(*nd1), ff(*nd2)
      fm = f1 % f2
      assert (fm.numerator, fm.denominator) == rm.as_tuple()
    if (expected is not None):
      assert rm.as_tuple() == expected
  check((2,3), (1,2), (1,6))
  check((2,3), (-1,2), (-1,3))
  check((-2,3), (1,2), (1,3))
  check((-2,3), (-1,2), (-1,6))
  for ln in xrange(-7,7+1):
    for rn in xrange(-9,9+1):
      if (rn == 0): continue
      check((ln,3), (rn,4))

def exercise_functions():
  assert rational.gcd(8,6) == 2
  assert rational.lcm(8,6) == 24

def exercise_python_code():
  r = rational.int
  assert rational.from_string("1") == 1
  assert rational.from_string("2/4").as_tuple() == (1,2)
  assert rational.vector((2,3,4), 3) == [r(d,3) for d in (2,3,4)]
  assert rational.vector((2,3,4), (3,4,5)) == [
    r(d,n) for d,n in zip((2,3,4), (3,4,5))]
  assert rational.lcm_denominators(array=[]) == 1
  assert rational.lcm_denominators(array=[r(3,4)]) == 4
  assert rational.lcm_denominators(array=[r(3,4), r(5,6)]) == 12

def run():
  exercise_int()
  exercise_functions()
  exercise_python_code()
  print "OK"

if (__name__ == "__main__"):
  run()

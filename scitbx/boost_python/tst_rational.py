from scitbx import rational
from scitbx.test_utils import approx_equal

def exercise_int():
  r = rational.int()
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
  assert (rational.int(2,3) + 4).as_tuple() == (14,3)
  assert (rational.int(2,3) - 4).as_tuple() == (-10,3)
  assert (rational.int(2,3) * 4).as_tuple() == (8,3)
  assert (rational.int(2,3) / 4).as_tuple() == (1,6)
  assert (5 + rational.int(2,3)).as_tuple() == (17,3)
  assert (5 - rational.int(2,3)).as_tuple() == (13,3)
  assert (5 * rational.int(2,3)).as_tuple() == (10,3)
  assert (5 / rational.int(2,3)).as_tuple() == (15,2)
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
  try: int(r)
  except: pass
  else: raise AssertionError, "Exception expected"

def exercise_functions():
  assert rational.gcd(8,6) == 2
  assert rational.lcm(8,6) == 24

def exercise_python_code():
  assert rational.from_string("1") == 1
  assert rational.from_string("2/4").as_tuple() == (1,2)
  assert rational.vector((2,3,4), 3) == [rational.int(d,3) for d in (2,3,4)]
  assert rational.vector((2,3,4), (3,4,5)) == [
    rational.int(d,n) for d,n in zip((2,3,4), (3,4,5))]

def run():
  exercise_int()
  exercise_functions()
  exercise_python_code()
  print "OK"

if (__name__ == "__main__"):
  run()

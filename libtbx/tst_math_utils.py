from libtbx.test_utils import approx_equal

def exercise_integer():
  from libtbx.math_utils import iround, iceil, ifloor, nearest_integer
  assert iround(0) == 0
  assert iround(1.4) == 1
  assert iround(-1.4) == -1
  assert iround(1.6) == 2
  assert iround(-1.6) == -2
  assert iceil(0) == 0
  assert iceil(1.1) == 2
  assert iceil(-1.1) == -1
  assert iceil(1.9) == 2
  assert iceil(-1.9) == -1
  assert ifloor(0) == 0
  assert ifloor(1.1) == 1
  assert ifloor(-1.1) == -2
  assert ifloor(1.9) == 1
  assert ifloor(-1.9) == -2
  for i in xrange(-3,3+1):
    assert nearest_integer(i+0.3) == i
    assert nearest_integer(i+0.7) == i+1

def exercise_logical():
  from libtbx.math_utils import does_imply, are_equivalent
  #
  assert does_imply(True, True)
  assert not does_imply(True, False)
  assert does_imply(False, True)
  assert does_imply(False, False)
  #
  assert are_equivalent(True, True)
  assert not are_equivalent(True, False)
  assert not are_equivalent(False, True)
  assert are_equivalent(False, False)

def exercise_nested_loop():
  from libtbx.math_utils import nested_loop as nl
  assert [list(i) for i in nl([])] == []
  assert [list(i) for i in nl([1])] == [[0]]
  assert [list(i) for i in nl([1], open_range=False)] == [[0], [1]]
  assert [list(i) for i in nl([3])] == [[0], [1], [2]]
  assert [list(i) for i in nl(begin=[-2], end=[3])] == [
    [-2], [-1], [0], [1], [2]]
  assert [list(i) for i in nl(begin=[-1], end=[1], open_range=False)] == [
    [-1], [0], [1]]
  assert [list(i) for i in nl(begin=[-2,4], end=[3,6])] == [
    [-2, 4], [-2, 5], [-1, 4], [-1, 5], [0, 4], [0, 5], [1, 4], [1, 5],
    [2, 4], [2, 5]]
  assert [list(i) for i in nl(begin=[-2,4], end=[3,6], open_range=False)] == [
    [-2, 4], [-2, 5], [-2, 6], [-1, 4], [-1, 5], [-1, 6], [0, 4], [0, 5],
    [0, 6], [1, 4], [1, 5], [1, 6], [2, 4], [2, 5], [2, 6], [3, 4], [3, 5],
    [3, 6]]
  assert [list(i) for i in nl(begin=[-1,0,-1], end=[1,2,1])] == [
    [-1, 0, -1], [-1, 0, 0], [-1, 1, -1], [-1, 1, 0], [0, 0, -1], [0, 0, 0],
    [0, 1, -1], [0, 1, 0]]

def exercise_next_permutation():
  from libtbx.math_utils import next_permutation
  seq = []
  assert next_permutation(seq) is False
  seq = [0]
  assert next_permutation(seq) is False
  seq = [0,1]
  assert next_permutation(seq)
  assert seq == [1, 0]
  assert not next_permutation(seq)
  assert seq == [0, 1]
  seq = [0,1,2]
  result = []
  while True:
    result.append(tuple(seq))
    if (not next_permutation(seq)):
      break
  assert result == [
    (0, 1, 2),
    (0, 2, 1),
    (1, 0, 2),
    (1, 2, 0),
    (2, 0, 1),
    (2, 1, 0)]
  assert seq == [0,1,2]
  expected_n = 1
  for m in xrange(1,7):
    expected_n *= m
    seq = range(m)
    n = 0
    while True:
      n += 1
      if (not next_permutation(seq)):
        break
    assert seq == range(m)
    assert n == expected_n

def exercise_random_permutation_in_place():
  from libtbx.math_utils import random_permutation_in_place
  import random
  random.seed(0)
  l = range(8)
  for i_trial in xrange(10):
    random_permutation_in_place(list=l)
    if (l != range(8)):
      break
  else:
    raise AssertionError
  assert sorted(l) == range(8)

def exercise_prime_factors_of():
  from libtbx.math_utils import prime_factors_of
  assert prime_factors_of(n=1) == []
  prime_set = set()
  for n in xrange(2, 100):
    primes = prime_factors_of(n)
    pp = 1
    for p in primes:
      pp *= p
    assert pp == n
    prime_set.update(primes)
    if (n == 30):
      assert prime_set == set([2,3,5,7,11,13,17,19,23,29])
  for n in prime_set:
    assert prime_factors_of(n) == [n]
  assert len(prime_set) == 25

def exercise_normalize_angle():
  from libtbx.math_utils import normalize_angle as n
  import math
  for deg,period in [(False, 2*math.pi), (True, 360.)]:
    assert approx_equal(n(0, deg=deg), 0, eps=1.e-12)
    assert approx_equal(n(1.e-8, deg=deg), 1.e-8, eps=1.e-12)
    assert approx_equal(n(-1.e-8, deg=deg), period-1.e-8, eps=1.e-12)
    assert approx_equal(n(1, deg=deg), 1, eps=1.e-12)
    assert approx_equal(n(-1, deg=deg), period-1, eps=1.e-12)
  assert approx_equal(n(1.e+8), 1.9426951384)
  assert approx_equal(n(-1.e+8), 4.34049016878)
  assert approx_equal(n(1.e+8, deg=True), 280)
  assert approx_equal(n(-1.e+8, deg=True), 80)

def exercise():
  exercise_integer()
  exercise_logical()
  exercise_nested_loop()
  exercise_next_permutation()
  exercise_random_permutation_in_place()
  exercise_prime_factors_of()
  exercise_normalize_angle()
  print "OK"

if (__name__ == "__main__"):
  exercise()

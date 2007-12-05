import math

def iround(x):
  if (x < 0): return int(x-0.5)
  return int(x+.5)

def iceil(x):
  return iround(math.ceil(x))

def ifloor(x):
  return iround(math.floor(x))

def does_imply(p,q):
  """ does p => q in the sense of logical implication? """
  return not p or q

def are_equivalent(p,q):
  """ does p <=> q in the sense of logical equivalence? """
  return does_imply(p,q) and does_imply(q,p)

def next_permutation(seq):
  """Emulation of C++ std::next_permutation:
  Treats all permutations of seq as a set of "dictionary" sorted
  sequences. Permutes the current sequence into the next one of this set.
  Returns true if there are more sequences to generate. If the sequence
  is the largest of the set, the smallest is generated and false returned.
"""
  if (len(seq) <= 1): return False
  i = len(seq) - 1
  while True:
    ii = i
    i -= 1
    if (seq[i] < seq[ii]):
      j = len(seq)
      while True:
        j -= 1
        if (seq[i] < seq[j]):
          break
      seq[i], seq[j] = seq[j], seq[i]
      tail = seq[ii:]
      del seq[ii:]
      tail.reverse()
      seq.extend(tail)
      return True
    if (i == 0):
      seq.reverse()
      return False

def exercise_next_permutation():
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

if (__name__ == "__main__"):
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

  assert does_imply(True, True)
  assert not does_imply(True, False)
  assert does_imply(False, True)
  assert does_imply(False, False)

  assert are_equivalent(True, True)
  assert not are_equivalent(True, False)
  assert not are_equivalent(False, True)
  assert are_equivalent(False, False)

  exercise_next_permutation()

  print "OK"

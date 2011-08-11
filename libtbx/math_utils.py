import math

def iround(x):
  if (x < 0): return int(x-0.5)
  return int(x+.5)

def iceil(x):
  return iround(math.ceil(x))

def ifloor(x):
  return iround(math.floor(x))

def nearest_integer(x):
  return ifloor(x+0.5)

def does_imply(p,q):
  """ does p => q in the sense of logical implication? """
  return not p or q

def are_equivalent(p,q):
  """ does p <=> q in the sense of logical equivalence? """
  return does_imply(p,q) and does_imply(q,p)

class nested_loop(object):

  def __init__(O, end, begin=None, open_range=True):
    if (begin is None):
      begin = [0] * len(end)
    else:
      assert len(begin) == len(end)
    if (not open_range):
      end = list(end)
      for i in xrange(len(end)):
        end[i] += 1
    for i in xrange(len(end)):
      assert end[i] >= begin[i]
    O.begin = begin
    O.end = end
    O.current = list(begin)
    for i in xrange(len(end)):
      if (end[i] > begin[i]):
        O.current[-1] -= 1
        break

  def __iter__(O):
    return O

  def next(O):
    b = O.begin
    e = O.end
    c = O.current
    i = len(c)
    while (i > 0):
      i -= 1
      c[i] += 1
      if (c[i] < e[i]): return c
      c[i] = b[i]
    raise StopIteration

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

def random_permutation_in_place(list):
  import random
  n = len(list)
  for i in xrange(n):
    j = random.randrange(n)
    list[i], list[j] = list[j], list[i]

def prime_factors_of(n):
  "http://butunclebob.com/ArticleS.UncleBob.ThePrimeFactorsKata"
  result = []
  candidate = 2
  while (n > 1):
    while (n % candidate == 0):
      result.append(candidate)
      n //= candidate
    candidate += 1
  return result

def normalize_angle(phi, deg=False, zero_centered=False):
  if (deg): period = 360
  else:     period = 2 * math.pi
  phi = math.fmod(phi, period)
  if (phi < 0): phi += period
  if (zero_centered and phi > period/2):
    phi -= period
  return phi

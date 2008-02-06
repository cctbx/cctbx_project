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

def normalize_angle(phi, deg=False):
  if (deg): period = 360
  else:     period = 2 * math.pi
  phi = math.fmod(phi, period)
  if (phi < 0.): phi += period
  return phi

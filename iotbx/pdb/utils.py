from __future__ import division
from six.moves import range
import string
from itertools import product

def all_chain_ids():
  """
  Returns all possible 2 character chain IDs for PDB format files.
  In general, returns single character chains first.
  Also tries to use all combinations of uppercase/numbers before going to lower case.
  """
  chars = string.ascii_uppercase+string.digits
  lowerchars = string.ascii_lowercase
  def permutations(iterable, r=None): # XXX This may go to libtbx or scitbx
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    for indices in product(list(range(n)), repeat=r):
      if(len(indices) == r):
        yield tuple(pool[i] for i in indices)
  both_char_upper = permutations(iterable = chars, r = 2)
  second_char_lower = product(chars, lowerchars)
  first_char_lower = product(lowerchars, chars)
  both_char_lower = permutations(iterable = lowerchars, r = 2)
  #result = [c for c in chars]+[c for c in lowerchars]+\
  result = [" "+c for c in chars]+[" "+c for c in lowerchars]+\
         ["".join(p) for p in both_char_upper]+\
         ["".join(p) for p in second_char_lower]+\
         ["".join(p) for p in first_char_lower]+\
         ["".join(p) for p in both_char_lower]
  return result

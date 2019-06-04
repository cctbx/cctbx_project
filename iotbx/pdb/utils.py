from __future__ import absolute_import, division, print_function
import string
from itertools import product
from six.moves import range

def __permutations(iterable, r=None): # XXX This may go to libtbx or scitbx
  pool = tuple(iterable)
  n = len(pool)
  r = n if r is None else r
  for indices in product(range(n), repeat=r):
    if(len(indices) == r):
      yield tuple(pool[i] for i in indices)

def all_chain_ids():
  """
  Returns all possible 2 character chain IDs for PDB format files.
  In general, returns single character chains first.
  Also tries to use all combinations of uppercase/numbers before going to lower case.
  """
  chars = string.ascii_uppercase+string.digits
  lowerchars = string.ascii_lowercase
  both_char_upper = __permutations(iterable = chars, r = 2)
  second_char_lower = product(chars, lowerchars)
  first_char_lower = product(lowerchars, chars)
  both_char_lower = __permutations(iterable = lowerchars, r = 2)
  #result = [c for c in chars]+[c for c in lowerchars]+\
  result = [" "+c for c in chars]+[" "+c for c in lowerchars]+\
         ["".join(p) for p in both_char_upper]+\
         ["".join(p) for p in second_char_lower]+\
         ["".join(p) for p in first_char_lower]+\
         ["".join(p) for p in both_char_lower]
  return result

def all_label_asym_ids(maximum_length=4):
  chars = string.ascii_uppercase
  rc = ["".join(c) for c in chars]
  for r in range(2, maximum_length+1):
    char_upper = __permutations(iterable = chars, r = r)
    rc += ["".join(p) for p in char_upper]
  return rc

if __name__ == '__main__':
  import time
  l=0
  p=1
  for r in range(1,5):
    t0=time.time()
    rc = all_label_asym_ids(maximum_length=r)
    p*=26
    l+=p
    print('%7d %7d %7d %5s %0.3fs' % (l,p,len(rc), rc[-1], time.time()-t0))
    assert len(rc)==l

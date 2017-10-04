from __future__ import division
import string
from itertools import product

def all_chain_ids():
  chars = string.ascii_letters+string.digits
  def permutations(iterable, r=None): # XXX This may go to libtbx or scitbx
    pool = tuple(iterable)
    n = len(pool)
    r = n if r is None else r
    for indices in product(range(n), repeat=r):
      if(len(set(indices)) == r):
        yield tuple(pool[i] for i in indices)
  result = permutations(iterable = chars, r = 2)
  return [c for c in chars]+["".join(p) for p in result]

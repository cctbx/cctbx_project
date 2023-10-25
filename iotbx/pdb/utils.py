from __future__ import absolute_import, division, print_function
import string
from itertools import product
from six.moves import range

class generate_n_char_string:
  """ Iterator to generate strings of length n_chars, using upper-case,
    lower-case and numbers as desired.
    Allows specialty sets of characters as well

  parameters:
    n_chars:  length of string to produce
    include_upper:  include upper-case letters
    include_lower:  include lower-case letters
    include_numbers:  include numbers
    include_special_chars: include special characters:
       []_,.;:"&<>()/\{}'`~!@#$%*|+-
    end_with_tilde:  return n_chars - 1 plus the character "~"
    reverse_order:  reverse the order so numbers 9-0, lower case z-a,
                  upper case Z-A

  returns:
    n_char string, new string on every next()
    None if no more strings to return

  Tested in iotbx/regression/tst_generate_n_char_string.py

  """
  def __init__(self, n_chars = 1,
      include_upper = True,
      include_lower = True,
      include_numbers = True,
      include_special_chars = False,
      end_with_tilde = False,
      reverse_order = False):
    self._end_with_tilde = end_with_tilde
    if self._end_with_tilde:
      self._n_chars = n_chars - 1
    else: # usual
      self._n_chars = n_chars

    all_chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    all_chars_lc = all_chars.lower()
    all_numbers = '0123456789'
    special_characters = """[]_,.;:"&<>()\/\{}'`~!@#$%*|+-"""
    self._tilde = """~"""

    self._all_everything = ""
    if include_upper:
       self._all_everything += all_chars
    if include_lower:
       self._all_everything += all_chars_lc
    if include_numbers:
       self._all_everything += all_numbers
    if include_special_chars:
       self._all_everything += special_characters

    if reverse_order:
      # Use them in reverse order
      x = []
      for i in self._all_everything:
        x.append(i)
      x.reverse()
      self._all_everything = "".join(x)

    self._n = len(self._all_everything)
    self._indices = []
    for k in range(self._n_chars):
      self._indices.append(0)
  def next(self):
    # Write out current text based on current indices
    value = ""
    for k in range(self._n_chars):
      value += self._all_everything[self._indices[k]]

    # Update indices
    for kk in range(self._n_chars):
      k = self._n_chars - kk - 1 # from last index to first
      self._indices[k] += 1
      if self._indices[k] < self._n: #  current index is in range
        break
      elif k == 0: # current index is out of range and is first index
        return None # no more available
      else: # current index is out of range but is not first index
        self._indices[k] = 0
    if self._end_with_tilde:
      return value + self._tilde
    else: # usual
      return value


def __permutations(iterable, r=None): # XXX This may go to libtbx or scitbx
  pool = tuple(iterable)
  n = len(pool)
  r = n if r is None else r
  for indices in product(range(n), repeat=r):
    if(len(indices) == r):
      yield tuple(pool[i] for i in indices)

def all_chain_ids():
  """
  Test is in
  iotbx/regression/tst_all_chain_ids.py
  There should not be leading whitespace for one letter chain ids.

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
  # result = [" "+c for c in chars]+[" "+c for c in lowerchars]+\
  result = [c for c in chars]+[c for c in lowerchars]+\
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

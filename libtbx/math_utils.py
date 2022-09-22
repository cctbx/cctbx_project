from __future__ import absolute_import, division, print_function
from builtins import object
import math
from scitbx import matrix
from six.moves import range

def cmp(x, y):
  """
  cmp(x, y) -> integer

  Return negative if x<y, zero if x==y, positive if x>y.
  """
  return (x > y) - (x < y)

def round2(x, d=0):
  '''
  Python 3 defaults to rounding to the nearest even number (round half to even),
  so this function keeps the Python 2 behavior, which is rounding half away
  from zero.

  References:
  https://docs.python.org/3.7/library/functions.html#round
  https://en.wikipedia.org/wiki/Rounding#Round_half_to_even

  https://docs.python.org/2.7/library/functions.html#round
  https://en.wikipedia.org/wiki/Rounding#Round_half_away_from_zero

  Function from:
  http://python3porting.com/differences.html#rounding-behavior
  '''
  p = 10 ** d
  if x > 0:
    return float(math.floor((x * p) + 0.5))/p
  else:
    return float(math.ceil((x * p) - 0.5))/p

def roundoff(val, precision=3, as_string=False):
  '''
  round off all floats in a list (or tuple) of list (or tuples)
  recursively using round2() defined above as in:
  >>> math_utils.roundoff( [12.3454, 7.4843, ["foo", (35.3581, -0.3856, [4.2769, 3.2147] )] ])
  [12.345, 7.484, ['foo', (35.358, -0.386, [4.277, 3.215])]]
  If value is less than 10**-precision or greater than 10**precision then return with scientific notation
  '''
  if isinstance(val, float):
    if math.isnan(val):
      return float("nan")
    if abs(val) < float("1e-%d" %precision) or abs(val) > float("9e%d" %precision):
      fstr = "%" + "%d" %precision
      fstr += ".%de" %precision
      val2str = fstr %val
      if as_string:
        return val2str
      return float(val2str)
    return round2(val, precision)
  if isinstance(val, list):
    for i,v in enumerate(val):
      val[i] = roundoff(v, precision)
  if isinstance(val, tuple):
    val = list(val)
    for i,v in enumerate(val):
      val[i] = roundoff(v, precision)
    val = tuple(val)
  if isinstance(val,  matrix.sqr):
    val = list(val)
    for i,v in enumerate(val):
      val[i] = roundoff(v, precision)
    val =  matrix.sqr(val)
  if isinstance(val,  matrix.rec):
    valn = val.n
    val = list(val)
    for i,v in enumerate(val):
      val[i] = roundoff(v, precision)
    val =  matrix.rec(elems=val, n=valn)
  return val

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
      for i in range(len(end)):
        end[i] += 1
    for i in range(len(end)):
      assert end[i] >= begin[i]
    O.begin = begin
    O.end = end
    O.current = list(begin)
    for i in range(len(end)):
      if (end[i] > begin[i]):
        O.current[-1] -= 1
        break

  def __iter__(O):
    return O

  def __next__(O):
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
  for i in range(n):
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

def median(values, sort = None):
  """ Return median value which is same as percentile_based_spread with
    cutoff of 0.5
  """
  return percentile_based_spread(values, pbs_fraction=0.5, sort = sort)

def percentile_based_spread(values, pbs_fraction=0.608, sort = None):
  """
  See Pozharski (2010) Acta. Cryst. D66, 970-978.  The default value of the
  pbs_fraction parameter is for 3D geometries, and should be adjusted as
  circumstances dictate.
  """
  if sort is None:
    if hasattr(values, 'size'):
      n = values.size()
    else:
      n = len(values)
    if n < 100000:  # sort is slow if more than 100000
      sort = True
    else:
      sort = False

  if sort:
    return percentile_based_spread_with_sort(values, pbs_fraction=pbs_fraction)
  else:
    return percentile_based_spread_with_selection(
        values, pbs_fraction=pbs_fraction)

def percentile_based_spread_with_sort(values, pbs_fraction=0.608):
  """
  See Pozharski (2010) Acta. Cryst. D66, 970-978.  The default value of the
  pbs_fraction parameter is for 3D geometries, and should be adjusted as
  circumstances dictate.
  """
  values = sorted(values)
  n = len(values)
  if (n == 0): return 0
  elif (n == 1) : return values[0]
  i_high = min(iceil(n * pbs_fraction),n-1)
  i_low = ifloor(n * pbs_fraction)
  if (i_high == i_low):
    return values[i_high]
  x_high = values[i_high]
  x_low = values[i_low]
  frac_high = i_high / n
  frac_low = i_low / n
  assert (frac_high > frac_low)
  frac_delta = (pbs_fraction - frac_low) / (frac_high - frac_low)
  x_frac = x_low + (frac_delta * (x_high - x_low))
  return x_frac

def percentile_based_spread_with_selection(values, pbs_fraction=0.608,
   tolerance = 0.0001):
  """
  See Pozharski (2010) Acta. Cryst. D66, 970-978.  The default value of the
  pbs_fraction parameter is for 3D geometries, and should be adjusted as
  circumstances dictate.
  This version uses selection to get the pbs within tolerance of true value
   if possible
  """
  from scitbx.array_family import flex
  values = flex.double(values)
  mmm = values.min_max_mean()
  low = mmm.min
  high = mmm.max
  max_tries = values.size()  # absolute limit
  last_value = low
  too_low = True
  for i in range(max_tries):
    if high - low < tolerance:
       break
    if too_low:
      working =  0.5 * (last_value + high)
    else:
      working =  0.5 * (last_value + low)
    frac = (values < working).count(True) / values.size()
    too_low = (frac < pbs_fraction)
    last_value = working
    if too_low:
      low = working
    else:
      high = working
  return last_value

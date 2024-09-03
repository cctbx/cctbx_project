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

def percentile_based_spread_with_selection(values, pbs_fraction=0.608, tolerance = 0.0001):
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

def mahalanobis_using_sklearn(x=None, data=None, cov=None, verbose=False):
  if cov is None:
    cov = covariance_using_sklearn(data, verbose=verbose)
  mahal_cov = cov.mahalanobis(x)
  if verbose: print(mahal_cov)
  return mahal_cov

def covariance_using_sklearn(data=None, choice='Empirical', return_numpy=True, verbose=False):
  from sklearn.covariance import EmpiricalCovariance, MinCovDet
  import numpy as np
  # fit a MCD robust estimator to data
  robust_cov = MinCovDet().fit(data)
  # fit a MLE estimator to data
  emp_cov = EmpiricalCovariance().fit(data)
  if verbose:
    print(
        "Estimated covariance matrix:\nMCD (Robust):\n{}\nMLE:\n{}".format(
            robust_cov.covariance_, emp_cov.covariance_
        )
    )
  # choose one
  if choice=='Empirical':
    cov = emp_cov
  elif choice=='Robust':
    cov = robust_cov
  if return_numpy:
    return np.asarray(cov.covariance_)
  return cov

def covariance_using_numpy(data):
  import numpy as np
  return np.cov(np.transpose(data))

def numpy_mean(data, index=0):
  import numpy as np
  return np.mean(data, index)

def mahalanobis_using_numpy(x, data=None, mu=None, cov=None, inv_cov=None, verbose=False):
  """
  Compute the Mahalanobis Distance between each row of x and the data
    x    : vector or matrix of data with, say, p columns.
    data : ndarray of the distribution from which Mahalanobis distance of each
           observation of x is to be computed.
    mu   : mean of data. If None, computer from data.
    cov  : covariance matrix (p x p) of the distribution. If None, computed from data.
    inv_cov : inverse of covariance matrix. If None, computed from cov
  """
  import numpy as np
  def _mahalanobis(x, data=None, mu=None, cov=None, inv_cov=None):
    if verbose:
      print('x',x)
      print('data',data)
    if mu is None:
      mu = np.mean(data, 0)
    else:
      mu = np.asarray(mu)
    if verbose: print('mu',mu)
    if inv_cov is None:
      if cov is None:
        cov = covariance_using_numpy(data)
      else:
        cov=np.asarray(cov)
      inv_cov = np.linalg.inv(cov)
    if verbose:
      print('Covariance')
      print(cov)
    if verbose:
      print('Inverse Covariance')
      print(inv_cov)
    x_minus_mn = x - mu
    if verbose:
      print('x_minus_mn')
      print(x_minus_mn)
    D_square = np.dot(np.dot(x_minus_mn, inv_cov), np.transpose(x_minus_mn))
    if verbose:
      print('D_square')
      print(D_square)
      print(D_square.diagonal())
    return D_square.diagonal()

  if verbose:
    print(x)
    print(data)
    print(mu)
    print(cov)
    print(inv_cov)
  if data is None:
    assert mu and (cov is not None or inv_cov is not None)
  rc = _mahalanobis(x, data=data, mu=mu, cov=cov, inv_cov=inv_cov)
  return rc

def mahalanobis_p_values_outlier_indices(x=None, data=None, mu=None, cov=None, inv_cov=None, verbose=False):
  rc = mahalanobis_p_values(x=x, data=data, mu=mu, cov=cov, inv_cov=inv_cov, verbose=verbose)
  return list(filter(lambda x: rc[x] <0.01, range(len(rc))))

def mahalanobis_p_values(x=None, data=None, mu=None, cov=None, inv_cov=None, verbose=False):
  from scipy.stats import chi2
  df=len(x[0])-1
  if verbose:
    print(chi2.ppf((1-0.01), df=df))
    #> 9.21 for df=2
  mahal_cov = mahalanobis(x=x, data=data, mu=mu, cov=cov, inv_cov=inv_cov, verbose=verbose)
  if verbose: print(mahal_cov)
  return(1 - chi2.cdf(mahal_cov, df))

def mahalanobis(x=None, data=None, mu=None, cov=None, inv_cov=None, verbose=False):
  if 1:
    return mahalanobis_using_numpy(x=x, data=data, mu=mu, cov=cov, inv_cov=inv_cov, verbose=verbose)
  else:
    return mahalanobis_using_sklearn(x=x, data=data, cov=cov, verbose=verbose)

from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from six.moves import range

class simple_weighted_correlation(object): # used for data merging
  def __init__(self,w, x, y, derivatives_wrt_y_depth=0):
    import math
    assert derivatives_wrt_y_depth==0 # no derivatives implemented presently

    sum_xx = flex.sum(w * x**2)
    sum_yy = flex.sum(w * y**2)
    sum_xy = flex.sum(w * x * y)
    sum_x = flex.sum(w * x)
    sum_y = flex.sum(w * y)
    sum_w = flex.sum(w)
    assert sum_w != 0
      # Linear fit y to x, i.e. find slope and offset such that
      # y = slope * x + offset, optimal in a least-squares sense.
      # see p. 105 in Bevington & Robinson, Data Reduction and Error Analysis for
      #   the Physical Sciences, 3rd edition.  New York: McGraw Hill (2003)
    DELTA = sum_w * sum_xx - sum_x**2
    assert DELTA != 0
    self.slope = (sum_w * sum_xy - sum_x * sum_y) / DELTA
    self.offset = (sum_xx * sum_y - sum_x * sum_xy) / DELTA
    self.corr = (sum_w * sum_xy - sum_x * sum_y) / (math.sqrt(sum_w * sum_xx - sum_x**2) *
                                                    math.sqrt(sum_w * sum_yy - sum_y**2))

def weighted_correlation(w, x, y, derivatives_wrt_y_depth=0):
  "http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation"
  assert derivatives_wrt_y_depth in [0,1,2]
  sum_w = flex.sum(w)
  assert sum_w != 0
  wxm = flex.sum(w * x) / sum_w
  wym = flex.sum(w * y) / sum_w
  xc = x - wxm
  yc = y - wym
  sum_wxy = flex.sum(w * xc * yc)
  sum_wxx = flex.sum(w * xc * xc)
  sum_wyy = flex.sum(w * yc * yc)
  cc_den_sq = sum_wxx * sum_wyy
  cc_den = cc_den_sq**0.5
  cc = sum_wxy / cc_den
  if (derivatives_wrt_y_depth == 0):
    return cc
  #
  # 1st derivatives w.r.t. y
  d_sum_wxy = w * xc
  d_sum_wyy = 2 * w * yc
  d_cc_den_sq = sum_wxx * d_sum_wyy
  d_cc_den = 1/2 * d_cc_den_sq / cc_den
  d_cc_num = d_sum_wxy * cc_den - sum_wxy * d_cc_den
  d_cc = d_cc_num / cc_den_sq
  if (derivatives_wrt_y_depth == 1):
    return cc, d_cc
  #
  # 2nd derivatives w.r.t. y
  d_yc = 1 - w / sum_w
  d2_sum_wyy = 2 * w * d_yc
  d2_cc_den_sq = sum_wxx * d2_sum_wyy
  d2_cc_den = 1/2 * (d2_cc_den_sq / cc_den - d_cc_den_sq * d_cc_den / cc_den_sq)
  d2_cc_num = -sum_wxy * d2_cc_den
  d2_cc = d2_cc_num / cc_den_sq - d_cc_num * d_cc_den_sq / cc_den_sq**2
  return cc, d_cc, d2_cc

def finite_difference_derivatives(w, x, y, depth, eps=1e-6):
  assert depth in [1,2]
  result = flex.double()
  for i in range(len(y)):
    fs = []
    y_orig = y[i]
    for signed_eps in [eps, -eps]:
      y[i] = y_orig + signed_eps
      if (depth == 1):
        fs.append(weighted_correlation(w, x, y))
      else:
        _, d_cc = weighted_correlation(w, x, y, derivatives_wrt_y_depth=1)
        fs.append(d_cc[i])
    y[i] = y_orig
    result.append((fs[0]-fs[1])/(2*eps))
  return result

def exercise():
  mt = flex.mersenne_twister(seed=0)
  sz = 12
  for i_trial in range(10):
    x = mt.random_double(size=sz)*5-1
    y = mt.random_double(size=sz)*3-1
    for i_w,w in enumerate([flex.double(sz, 1), mt.random_double(size=sz)*7]):
      cc, d_cc, d2_cc = weighted_correlation(
        w, x, y, derivatives_wrt_y_depth=2)
      if (i_w == 0):
        cc_w1 = flex.linear_correlation(x, y).coefficient()
        assert approx_equal(cc, cc_w1)
      simple_cc = simple_weighted_correlation(w, x, y).corr
      assert approx_equal(simple_cc, cc)
      d_cc_fd = finite_difference_derivatives(w, x, y, depth=1)
      assert approx_equal(d_cc, d_cc_fd)
      d2_cc_fd = finite_difference_derivatives(w, x, y, depth=2)
      assert approx_equal(d2_cc, d2_cc_fd)

def run(args):
  assert len(args) == 0
  exercise()
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

from __future__ import division
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal

def weighted_correlation(w, x, y, compute_d_cc_d_y=False):
  "http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient#Calculating_a_weighted_correlation"
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
  if (not compute_d_cc_d_y):
    return cc
  #
  # derivatives w.r.t. y
  d_wym = w / sum_w
  d_yc = d_wym.matrix_outer_product(-1)
  d_yc.matrix_diagonal_add_in_place(1)
  d_sum_wxy = d_yc.matrix_multiply(w * xc)
  d_sum_wyy = d_yc.matrix_multiply(2 * w * yc)
  d_cc_den_sq = sum_wxx * d_sum_wyy
  d_cc_den = 1 / (2 * cc_den) * d_cc_den_sq
  d_cc = (d_sum_wxy * cc_den - sum_wxy * d_cc_den) / cc_den**2
  return cc, d_cc

def d_cc_d_y_fd(w, x, y, eps=1e-6):
  result = flex.double()
  for i in xrange(len(y)):
    fs = []
    y_orig = y[i]
    for signed_eps in [eps, -eps]:
      y[i] = y_orig + signed_eps
      fs.append(weighted_correlation(w, x, y))
    y[i] = y_orig
    result.append((fs[0]-fs[1])/(2*eps))
  return result

def exercise():
  mt = flex.mersenne_twister(seed=0)
  sz = 12
  for i_trial in xrange(10):
    x = mt.random_double(size=sz)*5-1
    y = mt.random_double(size=sz)*3-1
    for i_w,w in enumerate([flex.double(sz, 1), mt.random_double(size=sz)*7]):
      cc, d_cc = weighted_correlation(w, x, y, compute_d_cc_d_y=True)
      if (i_w == 0):
        cc_w1 = flex.linear_correlation(x, y).coefficient()
        assert approx_equal(cc, cc_w1)
      d_cc_fd = d_cc_d_y_fd(w, x, y)
      assert approx_equal(d_cc, d_cc_fd)

def run(args):
  assert len(args) == 0
  exercise()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

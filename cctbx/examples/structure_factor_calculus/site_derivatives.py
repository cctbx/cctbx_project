from __future__ import absolute_import, division, print_function
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import math
from six.moves import cStringIO as StringIO
import sys
from six.moves import range

flex.set_random_seed(0)

class cos_alpha:

  def __init__(self, site, ops, hkl):
    self.site = site
    self.ops = ops
    self.hkl = hkl

  def f(self):
    result = 0
    for op in self.ops:
      op_site = op * self.site
      result += math.cos(self.hkl.dot(op_site))
    return result

  def d_site(self):
    result = flex.double(3, 0)
    for op in self.ops:
      op_site = op * self.site
      hkl_op_site = self.hkl.dot(op_site)
      d_op_site = -math.sin(hkl_op_site) * self.hkl
      gtmx = op.transpose()
      d_site = gtmx * d_op_site.transpose()
      result += flex.double(d_site)
    return result

  def d2_site(self):
    result = flex.double(3*3, 0)
    d_alpha_d_site = self.hkl.outer_product()
    for op in self.ops:
      op_site = op * self.site
      hkl_op_site = self.hkl.dot(op_site)
      d2_op_site = -math.cos(hkl_op_site) * d_alpha_d_site
      result += flex.double(op.transpose() * d2_op_site * op)
    return result

def d_cos_alpha_d_site_finite(site, ops, hkl, eps=1.e-6):
  result = flex.double()
  site_eps = list(site)
  for ip in range(3):
    vs = []
    for signed_eps in [eps, -eps]:
      site_eps[ip] = site[ip] + signed_eps
      ca = cos_alpha(site=matrix.col(site_eps), ops=ops, hkl=hkl)
      vs.append(ca.f())
    site_eps[ip] = site[ip]
    result.append((vs[0]-vs[1])/(2*eps))
  return result

def d2_cos_alpha_d_site_finite(site, ops, hkl, eps=1.e-6):
  result = flex.double()
  site_eps = list(site)
  for ip in range(3):
    vs = []
    for signed_eps in [eps, -eps]:
      site_eps[ip] = site[ip] + signed_eps
      ca = cos_alpha(site=matrix.col(site_eps), ops=ops, hkl=hkl)
      vs.append(ca.d_site())
    site_eps[ip] = site[ip]
    result.extend((vs[0]-vs[1])/(2*eps))
  return result

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  for i_trial in range(100):
    ops = []
    for i in range(3):
      ops.append(matrix.sqr(flex.random_double(size=9, factor=4)-2))
    site = matrix.col(flex.random_double(size=3, factor=4)-2)
    hkl = matrix.row(flex.random_double(size=3, factor=4)-2)
    ca = cos_alpha(site=site, ops=ops, hkl=hkl)
    grads_fin = d_cos_alpha_d_site_finite(site=site, ops=ops, hkl=hkl)
    print("grads_fin:", list(grads_fin), file=out)
    grads_ana = ca.d_site()
    print("grads_ana:", list(grads_ana), file=out)
    assert approx_equal(grads_ana, grads_fin)
    curvs_fin = d2_cos_alpha_d_site_finite(site=site, ops=ops, hkl=hkl)
    print("curvs_fin:", list(curvs_fin), file=out)
    curvs_ana = ca.d2_site()
    print("curvs_ana:", list(curvs_ana), file=out)
    assert approx_equal(curvs_ana, curvs_fin)
    print(file=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

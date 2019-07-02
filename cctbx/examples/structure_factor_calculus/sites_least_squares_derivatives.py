from __future__ import absolute_import, division, print_function
from cctbx.examples.exp_i_alpha_derivatives import least_squares
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
import cmath
from six.moves import cStringIO as StringIO
import sys
from six.moves import range
from six.moves import zip

flex.set_random_seed(0)

class exp_i_hx:

  def __init__(self, sites, ops, hkl):
    self.sites = sites
    self.ops = ops
    self.hkl = hkl

  def f(self):
    result = 0
    for site in self.sites:
      for op in self.ops:
        op_site = op * site
        result += cmath.exp(1j*self.hkl.dot(op_site))
    return result

  def d_sites(self):
    result = []
    for site in self.sites:
      d_site = flex.complex_double(3, 0)
      for op in self.ops:
        op_site = op * site
        hkl_op_site = self.hkl.dot(op_site)
        d_op_site = cmath.exp(1j*self.hkl.dot(op_site)) * 1j * self.hkl
        gtmx = op.transpose()
        d_site += flex.complex_double(gtmx * d_op_site.transpose())
      result.append(d_site)
    return result

  def d2_sites(self):
    hkl_outer = self.hkl.outer_product()
    for js,site in enumerate(self.sites):
      d2_site = flex.complex_double(flex.grid(3,3), 0)
      for op in self.ops:
        op_site = op * site
        hkl_op_site = self.hkl.dot(op_site)
        d2_op_site = cmath.exp(1j*self.hkl.dot(op_site)) * (-1) * hkl_outer
        d2_site += flex.complex_double(op.transpose() * d2_op_site * op)
      yield d2_site

  def d_target_d_sites(self, target):
    da, db = target.da(), target.db()
    return flex.double([[da * d.real + db * d.imag
      for d in d_site]
        for d_site in self.d_sites()])

  def d2_target_d_sites(self, target):
    result = []
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    ds = self.d_sites()
    d2s = iter(self.d2_sites())
    for di0,d2i in zip(ds, d2s):
      d2ij_iter = iter(d2i)
      for di in di0:
        row = []
        for dj0 in ds:
          for dj in dj0:
            sum = daa * di.real * dj.real \
                + dbb * di.imag * dj.imag \
                + dab * (di.real * dj.imag + di.imag * dj.real)
            if (di0 is dj0):
              d2ij = next(d2ij_iter)
              sum += da * d2ij.real + db * d2ij.imag
            row.append(sum)
        result.append(row)
    return flex.double(result)

def d_exp_i_hx_d_sites_finite(sites, ops, hkl, obs, eps=1.e-8):
  result = flex.double()
  sites_eps = list(sites)
  for js,site in enumerate(sites):
    site_eps = list(site)
    for jp in range(3):
      vs = []
      for signed_eps in [eps, -eps]:
        site_eps[jp] = site[jp] + signed_eps
        sites_eps[js] = matrix.col(site_eps)
        sf = exp_i_hx(sites=sites_eps, ops=ops, hkl=hkl)
        tf = least_squares(obs=obs, calc=sf.f())
        vs.append(tf.f())
      result.append((vs[0]-vs[1])/(2*eps))
    sites_eps[js] = sites[js]
  return result

def d2_exp_i_hx_d_sites_finite(sites, ops, hkl, obs, eps=1.e-8):
  result = flex.double()
  sites_eps = list(sites)
  for js,site in enumerate(sites):
    site_eps = list(site)
    for jp in range(3):
      vs = []
      for signed_eps in [eps, -eps]:
        site_eps[jp] = site[jp] + signed_eps
        sites_eps[js] = matrix.col(site_eps)
        sf = exp_i_hx(sites=sites_eps, ops=ops, hkl=hkl)
        tf = least_squares(obs=obs, calc=sf.f())
        vs.append(sf.d_target_d_sites(target=tf))
      result.extend(((vs[0]-vs[1])/(2*eps)).as_1d())
    sites_eps[js] = sites[js]
  return result

def compare_derivatives(ana, fin):
  s = max(1, flex.max(flex.abs(ana)))
  assert approx_equal(ana/s, fin/s)

def exercise(args):
  verbose =  "--verbose" in args
  if (not verbose):
    out = StringIO()
  else:
    out = sys.stdout
  for i_trial in range(10):
    for n_sites in range(2,5+1):
      ops = []
      for i in range(3):
        ops.append(matrix.sqr(flex.random_double(size=9, factor=2)-1))
      sites = []
      for i in range(n_sites):
        sites.append(matrix.col(flex.random_double(size=3, factor=4)-2))
      hkl = matrix.row(flex.random_double(size=3, factor=4)-2)
      sf = exp_i_hx(sites=sites, ops=ops, hkl=hkl)
      for obs_factor in [1, 1.1]:
        obs = abs(sf.f()) * obs_factor
        grads_fin = d_exp_i_hx_d_sites_finite(
          sites=sites, ops=ops, obs=obs, hkl=hkl)
        print("grads_fin:", list(grads_fin), file=out)
        tf = least_squares(obs=obs, calc=sf.f())
        grads_ana = sf.d_target_d_sites(target=tf)
        print("grads_ana:", list(grads_ana), file=out)
        compare_derivatives(grads_ana, grads_fin)
        curvs_fin = d2_exp_i_hx_d_sites_finite(
          sites=sites, ops=ops, obs=obs, hkl=hkl)
        print("curvs_fin:", list(curvs_fin), file=out)
        curvs_ana = sf.d2_target_d_sites(target=tf)
        print("curvs_ana:", list(curvs_ana), file=out)
        compare_derivatives(curvs_ana, curvs_fin)
        print(file=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

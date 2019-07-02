from __future__ import absolute_import, division, print_function
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import tensor_rank_2_gradient_transform_matrix
from libtbx.test_utils import approx_equal
import math
from six.moves import cStringIO as StringIO
import sys
from six.moves import range

flex.set_random_seed(0)

mtps = -2 * math.pi**2

class debye_waller:

  def __init__(self, u, ops, hkl):
    self.u = u
    self.ops = ops
    self.hkl = hkl

  def f(self):
    result = 0
    for op in self.ops:
      op_u = (op*matrix.sym(sym_mat3=self.u)*op.transpose()).as_sym_mat3()
      huh = (matrix.row(self.hkl) \
          * matrix.sym(sym_mat3=op_u)).dot(matrix.col(self.hkl))
      result += math.exp(mtps * huh)
    return result

  def d_u(self):
    result = flex.double(6, 0)
    h,k,l = self.hkl
    d_exp_huh_d_u = matrix.col([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    for op in self.ops:
      op_u = (op*matrix.sym(sym_mat3=self.u)*op.transpose()).as_sym_mat3()
      huh = (matrix.row(self.hkl) \
          * matrix.sym(sym_mat3=op_u)).dot(matrix.col(self.hkl))
      d_op_u = math.exp(mtps * huh) * mtps * d_exp_huh_d_u
      gtmx = tensor_rank_2_gradient_transform_matrix(op)
      d_u = gtmx.matrix_multiply(flex.double(d_op_u))
      result += d_u
    return result

  def d2_u(self):
    result = flex.double(flex.grid(6,6), 0)
    h,k,l = self.hkl
    d_exp_huh_d_u = flex.double([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    d2_exp_huh_d_uu = d_exp_huh_d_u.matrix_outer_product(d_exp_huh_d_u)
    for op in self.ops:
      op_u = (op*matrix.sym(sym_mat3=self.u)*op.transpose()).as_sym_mat3()
      huh = (matrix.row(self.hkl) \
          * matrix.sym(sym_mat3=op_u)).dot(matrix.col(self.hkl))
      d2_op_u = math.exp(mtps * huh) * mtps**2 * d2_exp_huh_d_uu
      gtmx = tensor_rank_2_gradient_transform_matrix(op)
      d2_u = gtmx.matrix_multiply(d2_op_u).matrix_multiply(
        gtmx.matrix_transpose())
      result += d2_u
    return result

def d_debye_waller_d_u_finite(u, ops, hkl, eps=1.e-8):
  result = flex.double()
  u_eps = list(u)
  for ip in range(6):
    vs = []
    for signed_eps in [eps, -eps]:
      u_eps[ip] = u[ip] + signed_eps
      a = debye_waller(u=matrix.col(u_eps), ops=ops, hkl=hkl)
      vs.append(a.f())
    u_eps[ip] = u[ip]
    result.append((vs[0]-vs[1])/(2*eps))
  return result

def d2_debye_waller_d_u_finite(u, ops, hkl, eps=1.e-8):
  result = flex.double()
  u_eps = list(u)
  for ip in range(6):
    vs = []
    for signed_eps in [eps, -eps]:
      u_eps[ip] = u[ip] + signed_eps
      a = debye_waller(u=matrix.col(u_eps), ops=ops, hkl=hkl)
      vs.append(a.d_u())
    u_eps[ip] = u[ip]
    result.extend((vs[0]-vs[1])/(2*eps))
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
  for i_trial in range(100):
    ops = []
    for i in range(3):
      ops.append(matrix.sqr(flex.random_double(size=9, factor=4)-2))
    u = matrix.col((flex.random_double(size=6, factor=2)-1)*1.e-3)
    hkl = matrix.row(flex.random_double(size=3, factor=4)-2)
    dw = debye_waller(u=u, ops=ops, hkl=hkl)
    grads_fin = d_debye_waller_d_u_finite(u=u, ops=ops, hkl=hkl)
    print("grads_fin:", list(grads_fin), file=out)
    grads_ana = dw.d_u()
    print("grads_ana:", list(grads_ana), file=out)
    compare_derivatives(grads_ana, grads_fin)
    curvs_fin = d2_debye_waller_d_u_finite(u=u, ops=ops, hkl=hkl)
    print("curvs_fin:", list(curvs_fin), file=out)
    curvs_ana = dw.d2_u()
    print("curvs_ana:", list(curvs_ana), file=out)
    compare_derivatives(curvs_ana, curvs_fin)
    print(file=out)
  print("OK")

if (__name__ == "__main__"):
  exercise(sys.argv[1:])

from cctbx.examples import g_exp_i_alpha_derivatives
from scitbx import matrix
import math

class parameters:

  def __init__(self, xyz, u, w, fp, fdp):
    self.xyz = tuple(xyz)
    self.u = u
    self.w = w
    self.fp = fp
    self.fdp = fdp

  def as_list(self):
    return list(self.xyz) + [self.u, self.w, self.fp, self.fdp]

  def as_g_alpha(self, hkl, d_star_sq):
    return g_exp_i_alpha_derivatives.parameters(
      g = self.w * math.exp(-2 * math.pi**2 * self.u * d_star_sq),
      ffp = 1 + self.fp,
      fdp = self.fdp,
      alpha = 2 * math.pi * matrix.col(self.xyz).dot(matrix.col(hkl)))

class gradients(parameters): pass

class curvatures:

  def __init__(self, uu, uw):
    self.uu = uu
    self.uw = uw

def pack_parameters(params):
  result = []
  for p in params:
    result.extend(p.as_list())
  return result

def pack_gradients(grads):
  return pack_parameters(grads)

class structure_factor:

  def __init__(self, hkl, d_star_sq, params):
    self.hkl = hkl
    self.d_star_sq = d_star_sq
    self.params = params

  def as_exp_i_sum(self):
    return g_exp_i_alpha_derivatives.g_exp_i_alpha_sum(
      params=[p.as_g_alpha(hkl=self.hkl, d_star_sq=self.d_star_sq)
        for p in self.params])

  def f(self):
    return self.as_exp_i_sum().f()

  def d_g_alpha_d_params(self):
    """Mathematica:
         alpha = 2 Pi {h,k,l}.{x,y,z}
         g = w Exp[-2 Pi^2 u dss]
         D[alpha,x]; D[alpha,y]; D[alpha,z]; D[g,u]; D[g,w]"
    """
    result = []
    c = -2 * math.pi**2 * self.d_star_sq
    d_xyz = 2 * math.pi * matrix.col(self.hkl)
    for p in self.params:
      e = math.exp(c * p.u)
      result.append(gradients(xyz=d_xyz, u=p.w*c*e, w=e, fp=1, fdp=1))
    return result

  def d2_g_alpha_d_params(self):
    """Mathematica:

         alpha = 2 Pi {h,k,l}.{x,y,z}
         g = w Exp[-2 Pi^2 u dss]
         D[alpha,x,x]; D[alpha,x,y]; D[alpha,x,z]; D[g,x,u]; D[g,x,w]"
         D[alpha,y,x]; D[alpha,y,y]; D[alpha,y,z]; D[g,y,u]; D[g,y,w]"
         D[alpha,z,x]; D[alpha,z,y]; D[alpha,z,z]; D[g,z,u]; D[g,z,w]"
         D[alpha,u,x]; D[alpha,u,y]; D[alpha,u,z]; D[g,u,u]; D[g,u,w]"
         D[alpha,w,x]; D[alpha,w,y]; D[alpha,w,z]; D[g,w,u]; D[g,w,w]"

    This curvature matrix is symmetric.
    All D[alpha, x|y|z, x|y|z|u|w] are 0.

       D[g,u,u] = (4 dss^2 Pi^4) w Exp[-2 Pi^2 u dss]
       D[g,u,w] = (-2 dss Pi^2)    Exp[-2 Pi^2 u dss]
       D[g,w,w] = 0
    """
    result = []
    c = -2 * math.pi**2 * self.d_star_sq
    for p in self.params:
      e = math.exp(c * p.u)
      result.append(curvatures(uu=c**2*p.w*e, uw=c*e))
    return result

  def d_target_d_params(self, target):
    result = []
    dts = self.as_exp_i_sum().d_target_d_params(target=target)
    ds = self.d_g_alpha_d_params()
    for dt,d in zip(dts, ds):
      result.append(gradients(
        xyz = dt.alpha * matrix.col(d.xyz),
        u = dt.g * d.u,
        w = dt.g * d.w,
        fp = dt.ffp,
        fdp = dt.fdp))
    return result

  def d2_target_d_params(self, target):
    """Combined application of chain rule and product rule.
       d_target_d_.. matrix:

         aa ag a' a"
         ga gg g' g"
         'a 'g '' '"
         "a "g "' ""

       Block in resulting matrix:

         xx xy xz xu xw x' x"
         yx yy yz yu yw y' y"
         zx zy zz zu zw z' z"
         ux uy uz uu uw '' '"
         wx wy wz wu ww "' ""
         'x 'y 'z 'u 'w '' '"
         "x "y "z "u "w "' ""
    """
    result = []
    exp_i_sum = self.as_exp_i_sum()
    dts = exp_i_sum.d_target_d_params(target=target)
    d2ti = iter(exp_i_sum.d2_target_d_params(target=target))
    ds = self.d_g_alpha_d_params()
    d2s = self.d2_g_alpha_d_params()
    idt2 = 0
    for dt,di,d2 in zip(dts, ds, d2s):
      # dx. dy. dz.
      d2ti0 = d2ti.next()
      for dxi in di.xyz:
        row = []; ra = row.append
        d2tij = iter(d2ti0)
        for dj in ds:
          d2t = d2tij.next()
          for dxj in dj.xyz:
            ra(d2t * dxi * dxj)
          d2t = d2tij.next()
          ra(d2t * dxi * dj.u)
          ra(d2t * dxi * dj.w)
          ra(d2tij.next() * dxi)
          ra(d2tij.next() * dxi)
        result.append(row)
      # d2u.
      row = []; ra = row.append
      d2ti0 = d2ti.next()
      d2tij = iter(d2ti0)
      for dj in ds:
        d2t = d2tij.next()
        for dxj in dj.xyz:
          ra(d2t * dxj * di.u)
        d2t = d2tij.next()
        ra(d2t * di.u * dj.u)
        if (di is dj): row[-1] += dt.g * d2.uu
        ra(d2t * di.u * dj.w)
        if (di is dj): row[-1] += dt.g * d2.uw
        ra(d2tij.next() * di.u)
        ra(d2tij.next() * di.u)
      result.append(row)
      # d2w.
      row = []; ra = row.append
      d2tij = iter(d2ti0)
      for dj in ds:
        d2t = d2tij.next()
        for dxj in dj.xyz:
          ra(d2t * dxj * di.w)
        d2t = d2tij.next()
        ra(d2t * di.w * dj.u)
        if (di is dj): row[-1] += dt.g * d2.uw
        ra(d2t * di.w * dj.w)
        ra(d2tij.next() * di.w)
        ra(d2tij.next() * di.w)
      result.append(row)
      # d2'. and d2"
      for ip in [0,1]:
        row = []; ra = row.append
        d2tij = iter(d2ti.next())
        for dj in ds:
          d2t = d2tij.next()
          for dxj in dj.xyz:
            ra(d2t * dxj)
          d2t = d2tij.next()
          ra(d2t * dj.u)
          ra(d2t * dj.w)
          ra(d2tij.next())
          ra(d2tij.next())
        result.append(row)
    return result

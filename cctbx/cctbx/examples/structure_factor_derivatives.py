from cctbx.examples import g_exp_i_alpha_derivatives
from scitbx import matrix
import cmath
import math

class parameters:

  def __init__(self, xyz, u, w):
    self.xyz = tuple(xyz)
    self.u = u
    self.w = w

  def as_list(self):
    return list(self.xyz) + [self.u, self.w]

  def as_g_alpha(self, hkl, d_star_sq):
    return g_exp_i_alpha_derivatives.parameters(
      alpha = 2 * math.pi * matrix.col(self.xyz).dot(matrix.col(hkl)),
      g = self.w * math.exp(-2 * math.pi**2 * self.u * d_star_sq))

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
      result.append(gradients(xyz=d_xyz, u=p.w*c*e, w=e))
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
        w = dt.g * d.w))
    return result

  def d2_target_d_params(self, target):
    """Combined application of chain rule and product rule:
         (d2t_aa da_d. + d2t_ga dg_d.) d. + dt d2
       Combinaion of
         gg ga
         ag aa
       and
         xx xy xz xu xw
         yx yy yz yu yw
         zx zy zz zu zw
         ux uy uz uu uw
         wx wy wz wu ww
    """
    result = []
    exp_i_sum = self.as_exp_i_sum()
    dts = exp_i_sum.d_target_d_params(target=target)
    d2ts = exp_i_sum.d2_target_d_params(target=target)
    ds = self.d_g_alpha_d_params()
    d2s = self.d2_g_alpha_d_params()
    idt2 = 0
    for dt,di,d2 in zip(dts, ds, d2s):
      # dx. dy. dz.
      idt2 += 1
      for dxi in di.xyz:
        row = []
        jdt2 = 0
        for dj in ds:
          jdt2 += 1
          d2t = d2ts[idt2][jdt2]
          for dxj in dj.xyz:
            row.append(d2t * dxi * dxj)
          jdt2 -= 1
          d2t = d2ts[idt2][jdt2]
          row.append(d2t * dxi * dj.u)
          row.append(d2t * dxi * dj.w)
          jdt2 += 2
        result.append(row)
      # d2u.
      idt2 -= 1
      row = []
      jdt2 = 0
      for dj in ds:
        jdt2 += 1
        d2t = d2ts[idt2][jdt2]
        for dxj in dj.xyz:
          row.append(d2t * dxj * di.u)
        jdt2 -= 1
        d2t = d2ts[idt2][jdt2]
        row.append(d2t * di.u * dj.u)
        if (di is dj): row[-1] += dt.g * d2.uu
        row.append(d2t * di.u * dj.w)
        if (di is dj): row[-1] += dt.g * d2.uw
        jdt2 += 2
      result.append(row)
      # d2w.
      row = []
      jdt2 = 0
      for dj in ds:
        jdt2 += 1
        d2t = d2ts[idt2][jdt2]
        for dxj in dj.xyz:
          row.append(d2t * dxj * di.w)
        jdt2 -= 1
        d2t = d2ts[idt2][jdt2]
        row.append(d2t * di.w * dj.u)
        if (di is dj): row[-1] += dt.g * d2.uw
        row.append(d2t * di.w * dj.w)
        jdt2 += 2
      result.append(row)
      idt2 += 2
    return result

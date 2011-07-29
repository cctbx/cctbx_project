from cctbx import xray
from cctbx.examples import g_exp_i_alpha_derivatives
from scitbx import matrix
from scitbx.array_family import flex
import math

def scatterer_as_list(self):
  return list(self.site) + [self.u_iso, self.occupancy, self.fp, self.fdp]

def scatterer_from_list(l):
  return xray.scatterer(
    site=l[:3],
    u=l[3],
    occupancy=l[4],
    scattering_type="const",
    fp=l[5],
    fdp=l[6])

def scatterer_as_g_alpha(scatterer, hkl, d_star_sq):
  return g_exp_i_alpha_derivatives.parameters(
    g = scatterer.occupancy
        * math.exp(-2 * math.pi**2 * scatterer.u_iso * d_star_sq),
    ffp = 1 + scatterer.fp,
    fdp = scatterer.fdp,
    alpha = 2 * math.pi * matrix.col(scatterer.site).dot(matrix.col(hkl)))

class gradients:

  def __init__(self, site, u_iso, occupancy, fp, fdp):
    self.site = site
    self.u_iso = u_iso
    self.occupancy = occupancy
    self.fp = fp
    self.fdp = fdp

class curvatures:

  def __init__(self, uu, uw):
    self.uu = uu
    self.uw = uw

def pack_gradients(grads):
  result = []
  for g in grads:
    result.extend(scatterer_as_list(g))
  return result

class structure_factor:

  def __init__(self, xray_structure, hkl):
    self.unit_cell = xray_structure.unit_cell()
    self.scatterers = xray_structure.scatterers()
    self.hkl = hkl
    self.d_star_sq = self.unit_cell.d_star_sq(hkl)

  def as_exp_i_sum(self):
    params = []
    for scatterer in self.scatterers:
      params.append(scatterer_as_g_alpha(
        scatterer=scatterer, hkl=self.hkl, d_star_sq=self.d_star_sq))
    return g_exp_i_alpha_derivatives.g_exp_i_alpha_sum(params=params)

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
    for scatterer in self.scatterers:
      e = math.exp(c * scatterer.u_iso)
      result.append(gradients(
        site=2*math.pi*matrix.col(self.hkl),
        u_iso=scatterer.occupancy*c*e,
        occupancy=e,
        fp=1,
        fdp=1))
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
    for scatterer in self.scatterers:
      e = math.exp(c * scatterer.u_iso)
      result.append(curvatures(uu=c**2*scatterer.occupancy*e, uw=c*e))
    return result

  def d_target_d_params(self, target):
    result = []
    dts = self.as_exp_i_sum().d_target_d_params(target=target)
    ds = self.d_g_alpha_d_params()
    for dt,d in zip(dts, ds):
      result.append(gradients(
        site = dt.alpha * matrix.col(d.site),
        u_iso = dt.g * d.u_iso,
        occupancy = dt.g * d.occupancy,
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
    for dt,di,d2 in zip(dts, ds, d2s):
      # dx. dy. dz.
      d2ti0 = d2ti.next()
      for dxi in di.site:
        row = []; ra = row.append
        d2tij = iter(d2ti0)
        for dj in ds:
          d2t = d2tij.next()
          for dxj in dj.site:
            ra(d2t * dxi * dxj)
          d2t = d2tij.next()
          ra(d2t * dxi * dj.u_iso)
          ra(d2t * dxi * dj.occupancy)
          ra(d2tij.next() * dxi)
          ra(d2tij.next() * dxi)
        result.append(row)
      # d2u.
      row = []; ra = row.append
      d2ti0 = d2ti.next()
      d2tij = iter(d2ti0)
      for dj in ds:
        d2t = d2tij.next()
        for dxj in dj.site:
          ra(d2t * dxj * di.u_iso)
        d2t = d2tij.next()
        ra(d2t * di.u_iso * dj.u_iso)
        if (di is dj): row[-1] += dt.g * d2.uu
        ra(d2t * di.u_iso * dj.occupancy)
        if (di is dj): row[-1] += dt.g * d2.uw
        ra(d2tij.next() * di.u_iso)
        ra(d2tij.next() * di.u_iso)
      result.append(row)
      # d2w.
      row = []; ra = row.append
      d2tij = iter(d2ti0)
      for dj in ds:
        d2t = d2tij.next()
        for dxj in dj.site:
          ra(d2t * dxj * di.occupancy)
        d2t = d2tij.next()
        ra(d2t * di.occupancy * dj.u_iso)
        if (di is dj): row[-1] += dt.g * d2.uw
        ra(d2t * di.occupancy * dj.occupancy)
        ra(d2tij.next() * di.occupancy)
        ra(d2tij.next() * di.occupancy)
      result.append(row)
      # d2'. and d2"
      for ip in [0,1]:
        row = []; ra = row.append
        d2tij = iter(d2ti.next())
        for dj in ds:
          d2t = d2tij.next()
          for dxj in dj.site:
            ra(d2t * dxj)
          d2t = d2tij.next()
          ra(d2t * dj.u_iso)
          ra(d2t * dj.occupancy)
          ra(d2tij.next())
          ra(d2tij.next())
        result.append(row)
    return result

class structure_factors:

  def __init__(self, xray_structure, miller_set):
    assert xray_structure.is_similar_symmetry(miller_set)
    self.xray_structure = xray_structure
    self.miller_indices = miller_set.indices()
    self.number_of_parameters = xray_structure.scatterers().size()*7

  def fs(self):
    result = flex.complex_double()
    for hkl in self.miller_indices:
      result.append(structure_factor(
        xray_structure=self.xray_structure, hkl=hkl).f())
    return result

  def f(self):
    return flex.sum(self.fs())

  def d_target_d_params(self, f_obs, target_type):
    result = flex.double(self.number_of_parameters, 0)
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      result += flex.double(
        pack_gradients(sf.d_target_d_params(target=target)))
    return result

  def d2_target_d_params(self, f_obs, target_type):
    result = flex.double(self.number_of_parameters**2, 0)
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      result += flex.double(sf.d2_target_d_params(target=target))
    return result

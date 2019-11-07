from __future__ import absolute_import, division, print_function
from cctbx import xray
from scitbx.math import tensor_rank_2_gradient_transform_matrix
from scitbx import matrix
from scitbx.array_family import flex
import cmath
import math
from six.moves import zip

def scatterer_as_list(self):
  if (self.flags.use_u_iso_only()):
    return list(self.site) + [self.u_iso, self.occupancy, self.fp, self.fdp]
  return list(self.site) + list(self.u_star) \
       + [self.occupancy, self.fp, self.fdp]

def scatterer_from_list(l):
  if (len(l) == 7):
    return xray.scatterer(
      site=l[:3],
      u=l[3],
      occupancy=l[4],
      scattering_type="const",
      fp=l[5],
      fdp=l[6])
  return xray.scatterer(
    site=l[:3],
    u=l[3:9],
    occupancy=l[9],
    scattering_type="const",
    fp=l[10],
    fdp=l[11])

class gradients:

  def __init__(self, site, u_iso, u_star, occupancy, fp, fdp):
    self.site = site
    self.u_iso = u_iso
    self.u_star = u_star
    self.flags = xray.scatterer_flags()
    self.flags.set_use_u(
      iso=(u_iso is not None),
      aniso=(u_star is not None))
    self.occupancy = occupancy
    self.fp = fp
    self.fdp = fdp

def pack_gradients(grads):
  result = []
  for g in grads:
    result.extend(scatterer_as_list(g))
  return result

mtps = -2 * math.pi**2

class structure_factor:

  def __init__(self, xray_structure, hkl):
    self.unit_cell = xray_structure.unit_cell()
    self.space_group = xray_structure.space_group()
    self.scatterers = xray_structure.scatterers()
    self.hkl = hkl
    self.d_star_sq = self.unit_cell.d_star_sq(hkl)

  def f(self):
    result = 0
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    for scatterer in self.scatterers:
      assert scatterer.scattering_type == "const"
      w = scatterer.occupancy
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = ffp + 1j * fdp
      for s in self.space_group:
        s_site = s * scatterer.site
        alpha = matrix.col(s_site).dot(tphkl)
        if (scatterer.flags.use_u_aniso()):
          r = s.r().as_rational().as_float()
          s_u_star_s = r*matrix.sym(sym_mat3=scatterer.u_star)*r.transpose()
          huh = (matrix.row(self.hkl) * s_u_star_s).dot(matrix.col(self.hkl))
          dw = math.exp(mtps * huh)
        e = cmath.exp(1j*alpha)
        result += w * dw * ff * e
    return result

  def df_d_params(self):
    result = []
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    h,k,l = self.hkl
    d_exp_huh_d_u_star = matrix.col([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    for scatterer in self.scatterers:
      assert scatterer.scattering_type == "const"
      w = scatterer.occupancy
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = ffp + 1j * fdp
      d_site = matrix.col([0,0,0])
      if (not scatterer.flags.use_u_aniso()):
        d_u_iso = 0
        d_u_star = None
      else:
        d_u_iso = None
        d_u_star = matrix.col([0,0,0,0,0,0])
      d_occ = 0
      d_fp = 0
      d_fdp = 0
      for s in self.space_group:
        r = s.r().as_rational().as_float()
        s_site = s * scatterer.site
        alpha = matrix.col(s_site).dot(tphkl)
        if (scatterer.flags.use_u_aniso()):
          s_u_star_s = r*matrix.sym(sym_mat3=scatterer.u_star)*r.transpose()
          huh = (matrix.row(self.hkl) * s_u_star_s).dot(matrix.col(self.hkl))
          dw = math.exp(mtps * huh)
        e = cmath.exp(1j*alpha)
        site_gtmx = r.transpose()
        d_site += site_gtmx * (
          w * dw * ff * e * 1j * tphkl)
        if (not scatterer.flags.use_u_aniso()):
          d_u_iso += w * dw * ff * e * mtps * self.d_star_sq
        else:
          u_star_gtmx = matrix.sqr(tensor_rank_2_gradient_transform_matrix(r))
          d_u_star += u_star_gtmx * (
            w * dw * ff * e * mtps * d_exp_huh_d_u_star)
        d_occ += dw * ff * e
        d_fp += w * dw * e
        d_fdp += w * dw * e * 1j
      result.append(gradients(
        site=d_site,
        u_iso=d_u_iso,
        u_star=d_u_star,
        occupancy=d_occ,
        fp=d_fp,
        fdp=d_fdp))
    return result

  def d2f_d_params(self):
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    tphkl_outer = tphkl.outer_product()
    h,k,l = self.hkl
    d_exp_huh_d_u_star = matrix.col([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    d2_exp_huh_d_u_star_u_star = d_exp_huh_d_u_star.outer_product()
    for scatterer in self.scatterers:
      assert scatterer.scattering_type == "const"
      w = scatterer.occupancy
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      d2_site_site = flex.complex_double(flex.grid(3,3), 0j)
      if (not scatterer.flags.use_u_aniso()):
        d2_site_u_iso = flex.complex_double(flex.grid(3,1), 0j)
        d2_site_u_star = None
      else:
        d2_site_u_iso = None
        d2_site_u_star = flex.complex_double(flex.grid(3,6), 0j)
      d2_site_occ = flex.complex_double(flex.grid(3,1), 0j)
      d2_site_fp = flex.complex_double(flex.grid(3,1), 0j)
      d2_site_fdp = flex.complex_double(flex.grid(3,1), 0j)
      if (not scatterer.flags.use_u_aniso()):
        d2_u_iso_u_iso = 0j
        d2_u_iso_occ = 0j
        d2_u_iso_fp = 0j
        d2_u_iso_fdp = 0j
      else:
        d2_u_star_u_star = flex.complex_double(flex.grid(6,6), 0j)
        d2_u_star_occ = flex.complex_double(flex.grid(6,1), 0j)
        d2_u_star_fp = flex.complex_double(flex.grid(6,1), 0j)
        d2_u_star_fdp = flex.complex_double(flex.grid(6,1), 0j)
      d2_occ_fp = 0j
      d2_occ_fdp = 0j
      for s in self.space_group:
        r = s.r().as_rational().as_float()
        s_site = s * scatterer.site
        alpha = matrix.col(s_site).dot(tphkl)
        if (scatterer.flags.use_u_aniso()):
          s_u_star_s = r*matrix.sym(sym_mat3=scatterer.u_star)*r.transpose()
          huh = (matrix.row(self.hkl) * s_u_star_s).dot(matrix.col(self.hkl))
          dw = math.exp(mtps * huh)
        e = cmath.exp(1j*alpha)
        site_gtmx = r.transpose()
        d2_site_site += flex.complex_double(
          site_gtmx *
            (w * dw * ff * e * (-1) * tphkl_outer)
               * site_gtmx.transpose())
        if (not scatterer.flags.use_u_aniso()):
          d2_site_u_iso += flex.complex_double(site_gtmx * (
            w * dw * ff * e * 1j * mtps * self.d_star_sq * tphkl))
        else:
          u_star_gtmx = matrix.sqr(tensor_rank_2_gradient_transform_matrix(r))
          d2_site_u_star += flex.complex_double(
              site_gtmx
            * ((w * dw * ff * e * 1j * tphkl).outer_product(
                mtps * d_exp_huh_d_u_star))
            * u_star_gtmx.transpose())
        d2_site_occ += flex.complex_double(site_gtmx * (
          dw * ff * e * 1j * tphkl))
        d2_site_fp += flex.complex_double(site_gtmx * (
          w * dw * e * 1j * tphkl))
        d2_site_fdp += flex.complex_double(site_gtmx * (
          w * dw * e * (-1) * tphkl))
        if (not scatterer.flags.use_u_aniso()):
          d2_u_iso_u_iso += w * dw * ff * e * (mtps * self.d_star_sq)**2
          d2_u_iso_occ += dw * ff * e * mtps * self.d_star_sq
          d2_u_iso_fp += w * dw * e * mtps * self.d_star_sq
          d2_u_iso_fdp += 1j * w * dw * e * mtps * self.d_star_sq
        else:
          d2_u_star_u_star += flex.complex_double(
              u_star_gtmx
            * (w * dw * ff * e * mtps**2 * d2_exp_huh_d_u_star_u_star)
            * u_star_gtmx.transpose())
          d2_u_star_occ += flex.complex_double(u_star_gtmx * (
            dw * ff * e * mtps * d_exp_huh_d_u_star))
          d2_u_star_fp += flex.complex_double(u_star_gtmx * (
            w * dw * e * mtps * d_exp_huh_d_u_star))
          d2_u_star_fdp += flex.complex_double(u_star_gtmx * (
            w * dw * 1j * e * mtps * d_exp_huh_d_u_star))
        d2_occ_fp += dw * e
        d2_occ_fdp += dw * e * 1j
      if (not scatterer.flags.use_u_aniso()):
        i_occ, i_fp, i_fdp, np = 4, 5, 6, 7
      else:
        i_occ, i_fp, i_fdp, np = 9, 10, 11, 12
      dp = flex.complex_double(flex.grid(np,np), 0j)
      paste = dp.matrix_paste_block_in_place
      paste(d2_site_site, 0,0)
      if (not scatterer.flags.use_u_aniso()):
        paste(d2_site_u_iso, 0,3)
        paste(d2_site_u_iso.matrix_transpose(), 3,0)
      else:
        paste(d2_site_u_star, 0,3)
        paste(d2_site_u_star.matrix_transpose(), 3,0)
      paste(d2_site_occ, 0,i_occ)
      paste(d2_site_occ.matrix_transpose(), i_occ,0)
      paste(d2_site_fp, 0,i_fp)
      paste(d2_site_fp.matrix_transpose(), i_fp,0)
      paste(d2_site_fdp, 0,i_fdp)
      paste(d2_site_fdp.matrix_transpose(), i_fdp,0)
      if (not scatterer.flags.use_u_aniso()):
        dp[3*7+3] = d2_u_iso_u_iso
        dp[3*7+4] = d2_u_iso_occ
        dp[4*7+3] = d2_u_iso_occ
        dp[3*7+5] = d2_u_iso_fp
        dp[5*7+3] = d2_u_iso_fp
        dp[3*7+6] = d2_u_iso_fdp
        dp[6*7+3] = d2_u_iso_fdp
      else:
        paste(d2_u_star_u_star, 3,3)
        paste(d2_u_star_occ, 3, 9)
        paste(d2_u_star_occ.matrix_transpose(), 9, 3)
        paste(d2_u_star_fp, 3, 10)
        paste(d2_u_star_fp.matrix_transpose(), 10, 3)
        paste(d2_u_star_fdp, 3, 11)
        paste(d2_u_star_fdp.matrix_transpose(), 11, 3)
      dp[i_occ*np+i_fp] = d2_occ_fp
      dp[i_fp*np+i_occ] = d2_occ_fp
      dp[i_occ*np+i_fdp] = d2_occ_fdp
      dp[i_fdp*np+i_occ] = d2_occ_fdp
      yield dp

  def d_target_d_params(self, target):
    da, db = target.da(), target.db()
    return flex.double([[da * d.real + db * d.imag
      for d in scatterer_as_list(d_scatterer)]
        for d_scatterer in self.df_d_params()])

  def d2_target_d_params(self, target):
    result = []
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    ds = self.df_d_params()
    d2s = self.d2f_d_params()
    for di0,d2i in zip(ds, d2s):
      d2ij_iter = iter(d2i)
      for di in scatterer_as_list(di0):
        row = []
        for dj0 in ds:
          for dj in scatterer_as_list(dj0):
            sum = daa * di.real * dj.real \
                + dbb * di.imag * dj.imag \
                + dab * (di.real * dj.imag + di.imag * dj.real)
            if (di0 is dj0):
              d2ij = next(d2ij_iter)
              sum += da * d2ij.real + db * d2ij.imag
            row.append(sum)
        result.append(row)
    return flex.double(result)

class structure_factors:

  def __init__(self, xray_structure, miller_set):
    assert xray_structure.is_similar_symmetry(miller_set)
    self.xray_structure = xray_structure
    self.miller_indices = miller_set.indices()
    np = 0
    for scatterer in xray_structure.scatterers():
      if (not scatterer.flags.use_u_aniso()):
        np += 7
      else:
        np += 12
    self.number_of_parameters = np

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
      result += sf.d_target_d_params(target=target)
    return result

  def d2_target_d_params(self, f_obs, target_type):
    np = self.number_of_parameters
    result = flex.double(flex.grid(np, np), 0)
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      result += sf.d2_target_d_params(target=target)
    return result

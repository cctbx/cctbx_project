from __future__ import generators
from cctbx import xray
from scitbx import matrix
from scitbx.array_family import flex
import cmath
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

class gradients:

  def __init__(self, site, u_iso, occupancy, fp, fdp):
    self.site = site
    self.u_iso = u_iso
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
      assert not scatterer.anisotropic_flag
      w = scatterer.occupancy
      dw = math.exp(mtps * scatterer.u_iso * self.d_star_sq)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      wdwff = w * dw * ff
      for s in self.space_group:
        s_site = s * scatterer.site
        alpha = matrix.col(s_site).dot(tphkl)
        result += wdwff * cmath.exp(1j*alpha)
    return result

  def df_d_params(self):
    result = []
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    for scatterer in self.scatterers:
      assert scatterer.scattering_type == "const"
      assert not scatterer.anisotropic_flag
      w = scatterer.occupancy
      dw = math.exp(mtps * scatterer.u_iso * self.d_star_sq)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      d_site = matrix.col([0,0,0])
      d_u_iso = 0
      d_occupancy = 0
      d_fp = 0
      d_fdp = 0
      for smx in self.space_group:
        site_s = smx * scatterer.site
        alpha = matrix.col(site_s).dot(tphkl)
        e = cmath.exp(1j*alpha)
        site_gtmx = smx.r().as_rational().as_float().transpose()
        d_site += site_gtmx * (w * dw * ff * e * 1j * tphkl)
        d_u_iso += w * dw * ff * e * mtps * self.d_star_sq
        d_occupancy += dw * ff * e
        d_fp += w * dw * e
        d_fdp += w * dw * e * 1j
      result.append(gradients(
        site=d_site,
        u_iso=d_u_iso,
        occupancy=d_occupancy,
        fp=d_fp,
        fdp=d_fdp))
    return result

  def d2f_d_params(self):
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    tphkl_outer = tphkl.outer_product()
    for scatterer in self.scatterers:
      assert scatterer.scattering_type == "const"
      assert not scatterer.anisotropic_flag
      w = scatterer.occupancy
      dw = math.exp(mtps * scatterer.u_iso * self.d_star_sq)
      ffp = 1 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      d2_site_site = flex.complex_double(flex.grid(3,3), 0j)
      d2_site_u_iso = flex.complex_double(flex.grid(1,3), 0j)
      d2_site_occupancy = flex.complex_double(flex.grid(1,3), 0j)
      d2_site_fp = flex.complex_double(flex.grid(1,3), 0j)
      d2_site_fdp = flex.complex_double(flex.grid(1,3), 0j)
      d2_u_iso_u_iso = 0j
      d2_u_iso_occupancy = 0j
      d2_u_iso_fp = 0j
      d2_u_iso_fdp = 0j
      d2_occupancy_fp = 0j
      d2_occupancy_fdp = 0j
      for smx in self.space_group:
        site_s = smx * scatterer.site
        alpha = matrix.col(site_s).dot(tphkl)
        e = cmath.exp(1j*alpha)
        site_gtmx = smx.r().as_rational().as_float().transpose()
        d2_site_site += flex.complex_double(
          site_gtmx *
            (w * dw * ff * e * (-1) * tphkl_outer)
               * site_gtmx.transpose())
        d2_site_u_iso += flex.complex_double(site_gtmx * (
          w * dw * ff * e * 1j * mtps * self.d_star_sq * tphkl))
        d2_site_occupancy += flex.complex_double(site_gtmx * (
          dw * ff * e * 1j * tphkl))
        d2_site_fp += flex.complex_double(site_gtmx * (
          w * dw * e * 1j * tphkl))
        d2_site_fdp += flex.complex_double(site_gtmx * (
          w * dw * e * (-1) * tphkl))
        d2_u_iso_u_iso += w * dw * ff * e * (mtps * self.d_star_sq)**2
        d2_u_iso_occupancy += dw * ff * e * mtps * self.d_star_sq
        d2_u_iso_fp += w * dw * e * mtps * self.d_star_sq
        d2_u_iso_fdp += 1j * w * dw * e * mtps * self.d_star_sq
        d2_occupancy_fp += dw * e
        d2_occupancy_fdp += dw * e * 1j
      dp = flex.complex_double(flex.grid(7,7), 0j)
      dp.matrix_paste_block_in_place(d2_site_site, 0,0)
      dp.matrix_paste_block_in_place(d2_site_u_iso.matrix_transpose(), 0,3)
      dp.matrix_paste_block_in_place(d2_site_u_iso, 3,0)
      dp.matrix_paste_block_in_place(d2_site_occupancy.matrix_transpose(), 0,4)
      dp.matrix_paste_block_in_place(d2_site_occupancy, 4,0)
      dp.matrix_paste_block_in_place(d2_site_fp.matrix_transpose(), 0,5)
      dp.matrix_paste_block_in_place(d2_site_fp, 5,0)
      dp.matrix_paste_block_in_place(d2_site_fdp.matrix_transpose(), 0,6)
      dp.matrix_paste_block_in_place(d2_site_fdp, 6,0)
      dp[3*7+3] = d2_u_iso_u_iso
      dp[3*7+4] = d2_u_iso_occupancy
      dp[4*7+3] = d2_u_iso_occupancy
      dp[3*7+5] = d2_u_iso_fp
      dp[5*7+3] = d2_u_iso_fp
      dp[3*7+6] = d2_u_iso_fdp
      dp[6*7+3] = d2_u_iso_fdp
      dp[4*7+5] = d2_occupancy_fp
      dp[5*7+4] = d2_occupancy_fp
      dp[4*7+6] = d2_occupancy_fdp
      dp[6*7+4] = d2_occupancy_fdp
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
              d2ij = d2ij_iter.next()
              sum += da * d2ij.real + db * d2ij.imag
            row.append(sum)
        result.append(row)
    return flex.double(result)

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

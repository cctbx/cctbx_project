from scitbx.math import tensor_rank_2_gradient_transform_matrix
from scitbx import matrix
from scitbx.array_family import flex
import cmath
import math

mtps = -2 * math.pi**2

class structure_factor:

  def __init__(self, xray_structure, hkl):
    self.unit_cell = xray_structure.unit_cell()
    self.space_group = xray_structure.space_group()
    self.scatterers = xray_structure.scatterers()
    self.site_symmetry_table = xray_structure.site_symmetry_table()
    self.scattering_type_registry = xray_structure.scattering_type_registry()
    self.hkl = hkl
    self.d_star_sq = self.unit_cell.d_star_sq(hkl)

  def f(self):
    result = 0
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    for scatterer in self.scatterers:
      w = scatterer.weight()
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      gaussian = self.scattering_type_registry.gaussian_not_optional(
        scattering_type=scatterer.scattering_type)
      f0 = gaussian.at_d_star_sq(self.d_star_sq)
      ffp = f0 + scatterer.fp
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
    tphkl = 2 * math.pi * matrix.col(self.hkl)
    h,k,l = self.hkl
    d_exp_huh_d_u_star = matrix.col([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    for i_scatterer,scatterer in enumerate(self.scatterers):
      site_symmetry_ops = None
      if (self.site_symmetry_table.is_special_position(i_scatterer)):
        site_symmetry_ops = self.site_symmetry_table.get(i_scatterer)
        site_constraints = site_symmetry_ops.site_constraints()
        if (scatterer.flags.use_u_aniso()):
          adp_constraints = site_symmetry_ops.adp_constraints()
      w = scatterer.weight()
      wwo = scatterer.weight_without_occupancy()
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      gaussian = self.scattering_type_registry.gaussian_not_optional(
        scattering_type=scatterer.scattering_type)
      f0 = gaussian.at_d_star_sq(self.d_star_sq)
      ffp = f0 + scatterer.fp
      fdp = scatterer.fdp
      ff = ffp + 1j * fdp
      d_site = matrix.col([0,0,0])
      if (not scatterer.flags.use_u_aniso()):
        d_u_iso = 0
        d_u_star = None
      else:
        d_u_iso = None
        d_u_star = matrix.col([0,0,0,0,0,0])
      d_occ = 0j
      d_fp = 0j
      d_fdp = 0j
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
        d_occ += wwo * dw * ff * e
        d_fp += w * dw * e
        d_fdp += w * dw * e * 1j
      if (site_symmetry_ops is not None):
        gsm = site_constraints.gradient_sum_matrix()
        gsm = matrix.rec(elems=gsm, n=gsm.focus())
        d_site = gsm * d_site
        if (scatterer.flags.use_u_aniso()):
          gsm = adp_constraints.gradient_sum_matrix()
          gsm = matrix.rec(elems=gsm, n=gsm.focus())
          d_u_star = gsm * d_u_star
      result = flex.complex_double(d_site)
      if (not scatterer.flags.use_u_aniso()):
        result.append(d_u_iso)
      else:
        result.extend(flex.complex_double(d_u_star))
      result.extend(flex.complex_double([d_occ, d_fp, d_fdp]))
      yield result

  def d2f_d_params(self):
    tphkl = 2 * math.pi * flex.double(self.hkl)
    tphkl_outer = tphkl.matrix_outer_product(tphkl) \
      .matrix_symmetric_as_packed_u()
    h,k,l = self.hkl
    d_exp_huh_d_u_star = flex.double([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    d2_exp_huh_d_u_star_u_star = d_exp_huh_d_u_star.matrix_outer_product(
      d_exp_huh_d_u_star).matrix_symmetric_as_packed_u()
    for i_scatterer,scatterer in enumerate(self.scatterers):
      site_symmetry_ops = None
      if (self.site_symmetry_table.is_special_position(i_scatterer)):
        site_symmetry_ops = self.site_symmetry_table.get(i_scatterer)
        site_constraints = site_symmetry_ops.site_constraints()
        if (scatterer.flags.use_u_aniso()):
          adp_constraints = site_symmetry_ops.adp_constraints()
      w = scatterer.weight()
      wwo = scatterer.weight_without_occupancy()
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      gaussian = self.scattering_type_registry.gaussian_not_optional(
        scattering_type=scatterer.scattering_type)
      f0 = gaussian.at_d_star_sq(self.d_star_sq)
      ffp = f0 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      d2_site_site = flex.complex_double(3*(3+1)//2, 0j)
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
        d2_u_star_u_star = flex.complex_double(6*(6+1)//2, 0j)
        d2_u_star_occ = flex.complex_double(flex.grid(6,1), 0j)
        d2_u_star_fp = flex.complex_double(flex.grid(6,1), 0j)
        d2_u_star_fdp = flex.complex_double(flex.grid(6,1), 0j)
      d2_occ_fp = 0j
      d2_occ_fdp = 0j
      for s in self.space_group:
        r = s.r().as_rational().as_float()
        s_site = s * scatterer.site
        alpha = tphkl.dot(flex.double(s_site))
        if (scatterer.flags.use_u_aniso()):
          s_u_star_s = r*matrix.sym(sym_mat3=scatterer.u_star)*r.transpose()
          huh = (matrix.row(self.hkl) * s_u_star_s).dot(matrix.col(self.hkl))
          dw = math.exp(mtps * huh)
        e = cmath.exp(1j*alpha)
        site_gtmx = flex.double(r.transpose())
        site_gtmx.reshape(flex.grid(3,3))
        d2_site_site += (w * dw * ff * e * (-1)) * (
          site_gtmx.matrix_multiply_packed_u_multiply_lhs_transpose(
            tphkl_outer))
        if (not scatterer.flags.use_u_aniso()):
          d2_site_u_iso += (w * dw * ff * e * 1j * mtps * self.d_star_sq) \
            * site_gtmx.matrix_multiply(tphkl)
        else:
          u_star_gtmx = tensor_rank_2_gradient_transform_matrix(r)
          d2_site_u_star += (w * dw * ff * e * 1j * mtps) \
            * site_gtmx.matrix_multiply(
                tphkl.matrix_outer_product(d_exp_huh_d_u_star)) \
                  .matrix_multiply(u_star_gtmx.matrix_transpose())
        site_gtmx_tphkl = site_gtmx.matrix_multiply(tphkl)
        d2_site_occ += (wwo * dw * ff * e * 1j) * site_gtmx_tphkl
        d2_site_fp += (w * dw * e * 1j) * site_gtmx_tphkl
        d2_site_fdp += (w * dw * e * (-1)) * site_gtmx_tphkl
        if (not scatterer.flags.use_u_aniso()):
          d2_u_iso_u_iso += w * dw * ff * e * (mtps * self.d_star_sq)**2
          d2_u_iso_occ += wwo * dw * ff * e * mtps * self.d_star_sq
          d2_u_iso_fp += w * dw * e * mtps * self.d_star_sq
          d2_u_iso_fdp += 1j * w * dw * e * mtps * self.d_star_sq
        else:
          d2_u_star_u_star +=(w * dw * ff * e * mtps**2) \
            * u_star_gtmx.matrix_multiply_packed_u_multiply_lhs_transpose(
                d2_exp_huh_d_u_star_u_star)
          u_star_gtmx_d_exp_huh_d_u_star = u_star_gtmx.matrix_multiply(
            d_exp_huh_d_u_star)
          d2_u_star_occ += (wwo * dw * ff * e * mtps) \
            * u_star_gtmx_d_exp_huh_d_u_star
          d2_u_star_fp += (w * dw * e * mtps) \
            * u_star_gtmx_d_exp_huh_d_u_star
          d2_u_star_fdp += (w * dw * 1j * e * mtps) \
            * u_star_gtmx_d_exp_huh_d_u_star
        d2_occ_fp += wwo * dw * e
        d2_occ_fdp += wwo * dw * e * 1j
      if (site_symmetry_ops is None):
        i_u = 3
      else:
        i_u = site_constraints.n_independent_params()
      if (not scatterer.flags.use_u_aniso()):
        i_occ = i_u + 1
      elif (site_symmetry_ops is None):
        i_occ = i_u + 6
      else:
        i_occ = i_u + adp_constraints.n_independent_params()
      i_fp, i_fdp, np = i_occ+1, i_occ+2, i_occ+3
      if (site_symmetry_ops is not None):
        gsm = site_constraints.gradient_sum_matrix()
        d2_site_site = gsm.matrix_multiply_packed_u_multiply_lhs_transpose(
          packed_u=d2_site_site)
        if (not scatterer.flags.use_u_aniso()):
          d2_site_u_iso = gsm.matrix_multiply(d2_site_u_iso)
        else:
          d2_site_u_star = gsm.matrix_multiply(d2_site_u_star)
        d2_site_occ = gsm.matrix_multiply(d2_site_occ)
        d2_site_fp = gsm.matrix_multiply(d2_site_fp)
        d2_site_fdp = gsm.matrix_multiply(d2_site_fdp)
        if (scatterer.flags.use_u_aniso()):
          gsm = adp_constraints.gradient_sum_matrix()
          d2_site_u_star = d2_site_u_star.matrix_multiply(
            gsm.matrix_transpose())
          d2_u_star_u_star = gsm \
            .matrix_multiply_packed_u_multiply_lhs_transpose(
              packed_u=d2_u_star_u_star)
          d2_u_star_occ = gsm.matrix_multiply(d2_u_star_occ)
          d2_u_star_fp = gsm.matrix_multiply(d2_u_star_fp)
          d2_u_star_fdp = gsm.matrix_multiply(d2_u_star_fdp)
      dp = flex.complex_double(flex.grid(np,np), 0j)
      paste = dp.matrix_paste_block_in_place
      paste(d2_site_site.matrix_packed_u_as_symmetric(), 0,0)
      if (not scatterer.flags.use_u_aniso()):
        paste(d2_site_u_iso, 0,i_u)
        paste(d2_site_u_iso.matrix_transpose(), i_u,0)
      else:
        paste(d2_site_u_star, 0,i_u)
        paste(d2_site_u_star.matrix_transpose(), i_u,0)
      paste(d2_site_occ, 0,i_occ)
      paste(d2_site_occ.matrix_transpose(), i_occ,0)
      paste(d2_site_fp, 0,i_fp)
      paste(d2_site_fp.matrix_transpose(), i_fp,0)
      paste(d2_site_fdp, 0,i_fdp)
      paste(d2_site_fdp.matrix_transpose(), i_fdp,0)
      if (not scatterer.flags.use_u_aniso()):
        dp[i_u*np+i_u] = d2_u_iso_u_iso
        dp[i_u*np+i_occ] = d2_u_iso_occ
        dp[i_occ*np+i_u] = d2_u_iso_occ
        dp[i_u*np+i_fp] = d2_u_iso_fp
        dp[i_fp*np+i_u] = d2_u_iso_fp
        dp[i_u*np+i_fdp] = d2_u_iso_fdp
        dp[i_fdp*np+i_u] = d2_u_iso_fdp
      else:
        paste(d2_u_star_u_star.matrix_packed_u_as_symmetric(), i_u, i_u)
        paste(d2_u_star_occ, i_u, i_occ)
        paste(d2_u_star_occ.matrix_transpose(), i_occ, i_u)
        paste(d2_u_star_fp, i_u, i_fp)
        paste(d2_u_star_fp.matrix_transpose(), i_fp, i_u)
        paste(d2_u_star_fdp, i_u, i_fdp)
        paste(d2_u_star_fdp.matrix_transpose(), i_fdp, i_u)
      dp[i_occ*np+i_fp] = d2_occ_fp
      dp[i_fp*np+i_occ] = d2_occ_fp
      dp[i_occ*np+i_fdp] = d2_occ_fdp
      dp[i_fdp*np+i_occ] = d2_occ_fdp
      yield dp

  def d2f_d_params_diag(self):
    tphkl = 2 * math.pi * flex.double(self.hkl)
    tphkl_outer = tphkl.matrix_outer_product(tphkl) \
      .matrix_symmetric_as_packed_u()
    h,k,l = self.hkl
    d_exp_huh_d_u_star = flex.double([h**2, k**2, l**2, 2*h*k, 2*h*l, 2*k*l])
    d2_exp_huh_d_u_star_u_star = d_exp_huh_d_u_star.matrix_outer_product(
      d_exp_huh_d_u_star).matrix_symmetric_as_packed_u()
    for i_scatterer,scatterer in enumerate(self.scatterers):
      site_symmetry_ops = None
      if (self.site_symmetry_table.is_special_position(i_scatterer)):
        site_symmetry_ops = self.site_symmetry_table.get(i_scatterer)
        site_constraints = site_symmetry_ops.site_constraints()
        if (scatterer.flags.use_u_aniso()):
          adp_constraints = site_symmetry_ops.adp_constraints()
      w = scatterer.weight()
      if (not scatterer.flags.use_u_aniso()):
        huh = scatterer.u_iso * self.d_star_sq
        dw = math.exp(mtps * huh)
      gaussian = self.scattering_type_registry.gaussian_not_optional(
        scattering_type=scatterer.scattering_type)
      f0 = gaussian.at_d_star_sq(self.d_star_sq)
      ffp = f0 + scatterer.fp
      fdp = scatterer.fdp
      ff = (ffp + 1j * fdp)
      d2_site_site = flex.complex_double(3*(3+1)//2, 0j)
      if (not scatterer.flags.use_u_aniso()):
        d2_u_iso_u_iso = 0j
      else:
        d2_u_star_u_star = flex.complex_double(6*(6+1)//2, 0j)
      for s in self.space_group:
        r = s.r().as_rational().as_float()
        s_site = s * scatterer.site
        alpha = tphkl.dot(flex.double(s_site))
        if (scatterer.flags.use_u_aniso()):
          s_u_star_s = r*matrix.sym(sym_mat3=scatterer.u_star)*r.transpose()
          huh = (matrix.row(self.hkl) * s_u_star_s).dot(matrix.col(self.hkl))
          dw = math.exp(mtps * huh)
        e = cmath.exp(1j*alpha)
        site_gtmx = flex.double(r.transpose())
        site_gtmx.reshape(flex.grid(3,3))
        d2_site_site += (w * dw * ff * e * (-1)) * (
          site_gtmx.matrix_multiply_packed_u_multiply_lhs_transpose(
            tphkl_outer))
        if (not scatterer.flags.use_u_aniso()):
          d2_u_iso_u_iso += w * dw * ff * e * (mtps * self.d_star_sq)**2
        else:
          u_star_gtmx = tensor_rank_2_gradient_transform_matrix(r)
          d2_u_star_u_star +=(w * dw * ff * e * mtps**2) \
            * u_star_gtmx.matrix_multiply_packed_u_multiply_lhs_transpose(
                d2_exp_huh_d_u_star_u_star)
      if (site_symmetry_ops is None):
        i_u = 3
      else:
        i_u = site_constraints.n_independent_params()
      if (not scatterer.flags.use_u_aniso()):
        i_occ = i_u + 1
      elif (site_symmetry_ops is None):
        i_occ = i_u + 6
      else:
        i_occ = i_u + adp_constraints.n_independent_params()
      np = i_occ+3
      if (site_symmetry_ops is not None):
        gsm = site_constraints.gradient_sum_matrix()
        d2_site_site = gsm.matrix_multiply_packed_u_multiply_lhs_transpose(
          packed_u=d2_site_site)
        if (scatterer.flags.use_u_aniso()):
          gsm = adp_constraints.gradient_sum_matrix()
          d2_u_star_u_star = gsm \
            .matrix_multiply_packed_u_multiply_lhs_transpose(
              packed_u=d2_u_star_u_star)
      #
      dpd = flex.complex_double(flex.grid(np,1), 0j)
      def paste(d, i):
        d.reshape(flex.grid(d.size(),1))
        dpd.matrix_paste_block_in_place(d, i,0)
      paste(d2_site_site.matrix_packed_u_diagonal(), 0)
      if (not scatterer.flags.use_u_aniso()):
        dpd[i_u] = d2_u_iso_u_iso
      else:
        paste(d2_u_star_u_star.matrix_packed_u_diagonal(), i_u)
      yield dpd

  def d_target_d_params(self, target):
    result = flex.double()
    da, db = target.da(), target.db()
    for d_scatterer in self.df_d_params():
      result.extend(flex.double([da * d.real + db * d.imag
        for d in d_scatterer]))
    return result

  def d2_target_d_params(self, target):
    result = []
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    ds = list(self.df_d_params())
    d2s = self.d2f_d_params()
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
              d2ij = d2ij_iter.next()
              sum += da * d2ij.real + db * d2ij.imag
            row.append(sum)
        result.append(row)
    return flex.double(result)

  def d2_target_d_params_diag(self, target):
    result = flex.double()
    da, db = target.da(), target.db()
    daa, dbb, dab = target.daa(), target.dbb(), target.dab()
    ds = self.df_d_params()
    d2sd = self.d2f_d_params_diag()
    for i_scatterer,(di0,d2id) in enumerate(zip(ds, d2sd)):
      for di,d2ij in zip(di0, d2id):
        sum = daa * di.real * di.real \
            + dbb * di.imag * di.imag \
            + dab * 2 * di.real * di.imag \
            + da * d2ij.real + db * d2ij.imag
        result.append(sum)
    return result

class structure_factors:

  def __init__(self, xray_structure, miller_set):
    assert xray_structure.is_similar_symmetry(miller_set)
    self.xray_structure = xray_structure
    self.miller_indices = miller_set.indices()

  def fs(self):
    result = flex.complex_double()
    for hkl in self.miller_indices:
      result.append(structure_factor(
        xray_structure=self.xray_structure, hkl=hkl).f())
    return result

  def f(self):
    return flex.sum(self.fs())

  def d_target_d_params(self, f_obs, target_type):
    result = None
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      contribution = sf.d_target_d_params(target=target)
      if (result is None): result = contribution
      else:                result += contribution
    return result

  def d2_target_d_params(self, f_obs, target_type):
    result = None
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      contribution = sf.d2_target_d_params(target=target)
      if (result is None): result = contribution
      else:                result += contribution
    return result

  def d2_target_d_params_diag(self, f_obs, target_type):
    result = None
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      contribution = sf.d2_target_d_params_diag(target=target)
      if (result is None): result = contribution
      else:                result += contribution
    return result

  def d2_target_d_params_diag_cpp(self, f_obs, target_type):
    da_db = flex.complex_double()
    daa_dbb_dab = flex.vec3_double()
    for hkl,obs in zip(self.miller_indices, f_obs.data()):
      sf = structure_factor(xray_structure=self.xray_structure, hkl=hkl)
      target = target_type(obs=obs, calc=sf.f())
      da_db.append(complex(target.da(), target.db()))
      daa_dbb_dab.append((target.daa(), target.dbb(), target.dab()))
    return self.xray_structure.grads_and_curvs_target_simple(
      miller_indices=f_obs.indices(), da_db=da_db, daa_dbb_dab=daa_dbb_dab)

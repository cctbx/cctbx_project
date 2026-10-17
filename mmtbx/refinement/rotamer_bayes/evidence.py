"""
Evidence from the diffraction data for moving a group of atoms: the change
in log-likelihood (dLL, natural log units) when the atoms of a side chain
move from their current sites to new ones, with the rest of the model, the
bulk-solvent mask, the scales and the error model held fixed.

Any target of mmtbx.f_model.manager can be used: its core is evaluated on
F_model + k dF_calc (k = k_isotropic k_anisotropic, as
d_target_d_f_calc_work). The targets are normalised by the number of working
reflections, so dLL = -n_work dT.

exact(...)      dF_calc by direct summation over the moved atoms, then the
                target: the rigid-move dLL itself.
predicted(...)  second order in dF_calc, without structure factors:

    dT ~ sum_h Re(g_h* dF_h) + 1/2 sum_h c_h |dF_h|^2

  The linear term is the occupancy derivative of the target for each atom
  at its old and new site (one FFT gradient calculation for a batch of
  atoms). The Hessian of the target in F_calc is block diagonal (one 2x2
  block per reflection); c_h is the mean of its radial and tangential
  curvatures, from finite differences of the gradient (target-independent),
  summed in resolution shells. With
    |dF_h|^2 ~ sum_jk q_j q_k f_j f_k T_j T_k sum_R cos(2 pi hR.(x_j - x_k))
  and the shell average of the symmetry sum, n_sym sin(x)/x with
  x = 2 pi s r_jk (cross terms between symmetry copies neglected), the
  quadratic term is a sum over atom pairs.
search_term(...) the linear term from precomputed gradient maps (carbon form
  factor, a few B values, interpolated), for searching chi angles cheaply.
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np
from libtbx import group_args
from cctbx.array_family import flex
from cctbx import maptbx, xray, adptbx

B_BINS = (5., 10., 20., 35., 55., 80., 120., 180., 270.)

class evidence(object):

  def __init__(self, fmodel, n_shells=40, eps=1.e-3):
    self.fmodel = fmodel
    self.xrs = fmodel.xray_structure
    self.functor = fmodel.target_functor()
    self.f_model = fmodel.f_model()
    self.k = (fmodel.k_isotropic() * fmodel.k_anisotropic())
    r0 = self.functor.core(f_calc=self.f_model, compute_gradients=True)
    self.t0 = r0.target_work()
    self.work = ~fmodel.r_free_flags().data()
    self.n_work = self.work.count(True)
    g_model = r0.gradients_work()
    k_work = self.k.select(self.work)
    self.g_calc = g_model * flex.complex_double(k_work,
      flex.double(k_work.size(), 0))
    self.f_obs_work = fmodel.f_obs_work()
    self.n_sym = self.xrs.space_group().order_z()
    self._curvatures(g_model, eps, n_shells)
    self._registry = self.xrs.scattering_type_registry()
    self._maps = {}

  # --- curvature of the target in F_calc, by resolution shell -------------
  def _curvatures(self, g_model, eps, n_shells):
    fm = self.f_model.data()
    amp = flex.abs(fm)
    small = amp < 1.e-6 * max(1.e-30, flex.max(amp))
    amp.set_selected(small, 1.)
    unit = fm / flex.complex_double(amp, flex.double(amp.size(), 0))
    unit.set_selected(small, 1 + 0j)
    step = amp * eps
    c = flex.double(self.n_work, 0)
    for direction in (unit, unit * 1j):
      d = direction * flex.complex_double(step, flex.double(step.size(), 0))
      r = self.functor.core(f_calc=self.f_model.customized_copy(
        data=fm + d), compute_gradients=True)
      dg = r.gradients_work() - g_model
      dw = d.select(self.work)
      num = flex.real(dg) * flex.real(dw) + flex.imag(dg) * flex.imag(dw)
      c += num / flex.pow2(flex.abs(dw)) / 2
    # on the F_calc scale
    c *= flex.pow2(self.k.select(self.work))
    self.c_work = c
    dss = self.f_obs_work.d_star_sq().data()
    order = flex.sort_permutation(dss)
    s2, csum = [], []
    n = dss.size()
    for b in range(n_shells):
      sel = order[b * n // n_shells:(b + 1) * n // n_shells]
      if(sel.size() == 0): continue
      s2.append(flex.mean(dss.select(sel)))
      csum.append(flex.sum(c.select(sel)))
    self.shell_s2 = np.array(s2)
    self.shell_c = np.array(csum)
    self.shell_s = np.sqrt(self.shell_s2)

  # --- form factors ---------------------------------------------------------
  def form_factors(self, scattering_types, b_isos, occupancies):
    """ q_j f_j(s) exp(-B_j s^2/4) at the shell centres (n_atoms, n_shells) """
    out = np.zeros((len(scattering_types), len(self.shell_s2)))
    cache = {}
    for i, (t, b, q) in enumerate(zip(scattering_types, b_isos, occupancies)):
      if(t not in cache):
        g = self._registry.gaussian_not_optional(scattering_type=t)
        cache[t] = np.array([g.at_d_star_sq(x) for x in self.shell_s2])
      out[i] = q * cache[t] * np.exp(-b * self.shell_s2 / 4)
    return out

  def quadratic(self, xyz, ff):
    """ 1/2 sum_h c_h |dF_h|^2 for atoms at xyz (n, 3) with signed form
    factors ff (n, n_shells) """
    d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2))
    x = 2 * math.pi * d[:, :, None] * self.shell_s[None, None, :]
    sinc = np.where(x > 1.e-8, np.sin(x) / np.where(x > 1.e-8, x, 1.), 1.)
    pair = ff[:, None, :] * ff[None, :, :] * sinc
    return 0.5 * self.n_sym * float((pair.sum(axis=(0, 1)) *
      self.shell_c).sum())

  # --- linear term ---------------------------------------------------------
  def occupancy_gradients(self, i_seqs_list, sites_list):
    """ dT/docc (per unit occupancy) for copies of the scatterers i_seqs at
    new sites, all in one FFT gradient calculation. Returns a list of
    numpy arrays matching sites_list. """
    all_i = flex.size_t()
    all_sites = flex.vec3_double()
    for i_seqs, sites in zip(i_seqs_list, sites_list):
      all_i.extend(flex.size_t([int(i) for i in i_seqs]))
      all_sites.extend(flex.vec3_double([tuple(s) for s in sites]))
    if(all_i.size() == 0): return []
    xrs = self.xrs.select(all_i)
    xrs.set_sites_cart(all_sites)
    sc = xrs.scatterers()
    sc.set_occupancies(flex.double(sc.size(), 1.))
    xray.set_scatterer_grad_flags(scatterers=sc, occupancy=True)
    sfg = self.fmodel.structure_factor_gradients_w(
      u_iso_refinable_params=None,
      d_target_d_f_calc=self.g_calc,
      xray_structure=xrs,
      n_parameters=0,
      miller_set=self.f_obs_work,
      algorithm=self.fmodel.sfg_params.algorithm,
      extra_params=self.fmodel.sfg_params.extra)
    g = np.array(sfg.d_target_d_occupancy())
    out, k = [], 0
    for i_seqs in i_seqs_list:
      out.append(g[k:k + len(i_seqs)]); k += len(i_seqs)
    return out

  def predicted(self, i_seqs, xyz_old, xyz_new, lin_old=None, lin_new=None):
    """ Predicted dLL for moving scatterers i_seqs from xyz_old to xyz_new.
    lin_old, lin_new: occupancy gradients if already computed. """
    sc = self.xrs.scatterers()
    occ = np.array([sc[int(i)].occupancy for i in i_seqs])
    if(lin_old is None or lin_new is None):
      lin_old, lin_new = self.occupancy_gradients(
        [i_seqs, i_seqs], [xyz_old, xyz_new])
    lin = float((occ * (lin_new - lin_old)).sum())
    quad = self.quadratic_move(i_seqs, xyz_old, xyz_new)
    return -self.n_work * (lin + quad)

  def atom_params(self, i_seqs):
    sc = self.xrs.scatterers()
    uc = self.xrs.unit_cell()
    types = [sc[int(i)].scattering_type for i in i_seqs]
    b = [sc[int(i)].u_iso_or_equiv(uc) * 8 * math.pi ** 2 for i in i_seqs]
    occ = [sc[int(i)].occupancy for i in i_seqs]
    return types, b, occ

  def quadratic_move(self, i_seqs, xyz_old, xyz_new, ff=None):
    if(ff is None):
      ff = self.form_factors(*self.atom_params(i_seqs))
    xyz = np.concatenate([xyz_old, xyz_new])
    return self.quadratic(xyz, np.concatenate([-ff, ff]))

  # --- exact -----------------------------------------------------------------
  def f_calc_of(self, i_seqs, xyz, occupancies=None):
    xrs = self.xrs.select(flex.size_t([int(i) for i in i_seqs]))
    xrs.set_sites_cart(flex.vec3_double([tuple(s) for s in xyz]))
    if(occupancies is not None):
      xrs.set_occupancies(flex.double([float(o) for o in occupancies]))
    return self.fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure=xrs, algorithm="direct").f_calc().data()

  def exact(self, i_seqs, xyz_old, xyz_new, f_old=None, occ_new=None,
            occ_old=None, mask_exclude=None):
    """ dLL for the rigid move (direct summation, then the target); occ_new:
    occupancies of the atoms at their new sites (default: as now). With
    mask_exclude (the residue's i_seqs) and setup_mask() done, the change
    of the bulk-solvent mask is included (local recomputation). """
    if(f_old is None): f_old = self.f_calc_of(i_seqs, xyz_old)
    d = self.f_calc_of(i_seqs, xyz_new, occupancies=occ_new) - f_old
    if(mask_exclude is not None and getattr(self, "mask_env", None)):
      sc = self.xrs.scatterers()
      oo = occ_old if occ_old is not None else [sc[int(i)].occupancy
        for i in i_seqs]
      on = occ_new if occ_new is not None else oo
      dm = self.delta_f_mask(i_seqs, xyz_old, xyz_new, oo, on, mask_exclude)
      if(dm is not None):
        d = d + dm   # dm carries k_mask (added to F_calc before k)
    d = d * flex.complex_double(self.k, flex.double(self.k.size(), 0))
    r = self.functor.core(f_calc=self.f_model.customized_copy(
      data=self.f_model.data() + d), compute_gradients=False)
    return -self.n_work * (r.target_work() - self.t0)

  # --- information -----------------------------------------------------------
  def llg_total(self):
    """ LLG of the whole model: LL(F_model) - LL(no model) """
    zero = self.f_model.customized_copy(
      data=flex.complex_double(self.f_model.size(), 0))
    r = self.functor.core(f_calc=zero, compute_gradients=False)
    return -self.n_work * (self.t0 - r.target_work())

  def scattering_power(self, i_seqs=None):
    """ sum of (occ f(0))^2 over non-hydrogen scatterers (all, or i_seqs) """
    sc = self.xrs.scatterers()
    if(i_seqs is None): i_seqs = range(sc.size())
    f0 = {}
    total = 0.
    for i in i_seqs:
      s = sc[int(i)]
      t = s.scattering_type
      if(t in ("H", "D")): continue
      if(t not in f0):
        f0[t] = self._registry.gaussian_not_optional(
          scattering_type=t).at_d_star_sq(0)
      total += (s.occupancy * f0[t]) ** 2
    return total

  def exact_removal(self, i_seqs, xyz_old, f_old=None):
    """ dLL for removing the atoms (direct summation, then the target) """
    if(f_old is None): f_old = self.f_calc_of(i_seqs, xyz_old)
    d = -f_old * flex.complex_double(self.k, flex.double(self.k.size(), 0))
    r = self.functor.core(f_calc=self.f_model.customized_copy(
      data=self.f_model.data() + d), compute_gradients=False)
    return -self.n_work * (r.target_work() - self.t0)

  # --- bulk-solvent mask change ----------------------------------------------
  def setup_mask(self, environment, spacing=0.4, cell=1.0):
    """ For the local mask correction: the mask parameters and atom radii of
    the f_model manager, and the environment (side_chain.environment, all
    present heavy atoms) for the atoms around a moving group. """
    from cctbx.masks import vdw_radii_from_xray_structure
    mp = self.fmodel.mask_params
    self.mask_probe = mp.solvent_radius
    self.mask_shrink = mp.shrink_truncation_radius
    self.mask_radii = np.array(vdw_radii_from_xray_structure(
      xray_structure=self.xrs))
    self.mask_env = environment
    self.mask_spacing = spacing
    self.mask_cell = cell
    self.k_mask = self.fmodel.k_masks()[0]
    from cctbx.eltbx import xray_scattering
    self._const = xray_scattering.gaussian(1.0)

  def _excluded(self, pts, xyz, radii):
    """ grid points within r + probe of any of the atoms """
    out = np.zeros(len(pts), dtype=bool)
    for x, r in zip(xyz, radii):
      rr = (r + self.mask_probe) ** 2
      d = pts - x
      out |= (d * d).sum(axis=1) < rr
    return out

  def _shrink(self, excluded, shape, spacing):
    """ solvent: points within the shrink radius of an accessible point """
    from scipy import ndimage
    dist = ndimage.distance_transform_edt(excluded.reshape(shape),
      sampling=spacing).ravel()
    return dist <= self.mask_shrink

  def mask_context(self, i_seqs, xyz_sets, occ_old, exclude_i_seqs,
                   spacing=None):
    """ Local mask machinery for one moving group: a grid covering the
    current sites and all those in xyz_sets (lists of sites of i_seqs), the
    exclusion by the atoms that stay (environment and the residue's other
    heavy atoms), and the solvent with the current atoms (those of occ_old
    > 0). """
    if(spacing is None): spacing = self.mask_spacing
    sc = self.xrs.scatterers()
    keep = np.array([k for k, i in enumerate(i_seqs) if
      sc[int(i)].scattering_type not in ("H", "D")], dtype=int)
    if(len(keep) == 0): return None
    i_mv = np.array([int(i_seqs[k]) for k in keep])
    r_mv = self.mask_radii[i_mv]
    allx = np.concatenate([np.array(x)[keep] for x in xyz_sets])
    reach = self.mask_probe + self.mask_shrink + 2.2
    lo = allx.min(axis=0) - reach; hi = allx.max(axis=0) + reach
    n = np.ceil((hi - lo) / spacing).astype(int) + 1
    pts = lo + np.mgrid[0:n[0], 0:n[1], 0:n[2]].reshape(3, -1).T * spacing
    centre = (lo + hi) / 2
    radius = np.sqrt(((hi - lo) ** 2).sum()) / 2 + 3.
    env_xyz, env_i = self.mask_env.near(centre, radius, exclude_i_seqs,
      return_i_seqs=True)
    mv = set(i_mv.tolist())
    others = [int(i) for i in exclude_i_seqs if int(i) not in mv
      and sc[int(i)].scattering_type not in ("H", "D")
      and sc[int(i)].occupancy > 0]
    fix_xyz = np.concatenate([env_xyz, np.array(
      self.xrs.sites_cart().select(flex.size_t(others))).reshape(-1, 3)])
    fix_r = np.concatenate([self.mask_radii[env_i], self.mask_radii[others]])
    excl_fixed = self._excluded(pts, fix_xyz, fix_r)
    ctx = group_args(pts=pts, shape=tuple(n), lo=lo, hi=hi, keep=keep,
      r_mv=r_mv, excl_fixed=excl_fixed, spacing=spacing,
      inside=np.all((pts > lo + self.mask_shrink + spacing) &
                    (pts < hi - self.mask_shrink - spacing), axis=1))
    return ctx

  def mask_solvent(self, ctx, xyz, occ):
    """ solvent on the context grid with the moving atoms at xyz (those of
    occ > 0) """
    x = np.array(xyz)[ctx.keep]; on = np.array(occ)[ctx.keep] > 0
    excl = ctx.excl_fixed | self._excluded(ctx.pts, x[on], ctx.r_mv[on])
    return self._shrink(excl, ctx.shape, ctx.spacing)

  def mask_cells(self, ctx, solvent_old, solvent_new):
    """ the change of solvent as pseudo-atoms (cells of mask_cell A):
    positions (n, 3) and signed volumes (n,) """
    dm = solvent_new.astype(int) - solvent_old.astype(int)
    dm[~ctx.inside] = 0
    sel = np.nonzero(dm)[0]
    if(len(sel) == 0): return np.zeros((0, 3)), np.zeros(0)
    p = ctx.pts[sel]; d = dm[sel]
    c = np.floor((p - ctx.lo) / self.mask_cell).astype(int)
    key = (c[:, 0] * 100000 + c[:, 1]) * 100000 + c[:, 2]
    key = key * 2 + (d > 0)
    u, inv = np.unique(key, return_inverse=True)
    cnt = np.bincount(inv)
    pos = np.stack([np.bincount(inv, weights=p[:, k]) for k in range(3)],
      axis=1) / cnt[:, None]
    sign = np.where(u % 2 == 1, 1., -1.)
    return pos, sign * cnt * ctx.spacing ** 3

  def _const_structure(self, pos, w):
    xrs = xray.structure(crystal_symmetry=self.xrs.crystal_symmetry())
    b = 8 * math.pi ** 2 * self.mask_cell ** 2 / 12
    uc = self.xrs.unit_cell()
    for x, o in zip(pos, w):
      xrs.add_scatterer(xray.scatterer(label="M", scattering_type="CONST",
        site=uc.fractionalize(tuple(x)), u=adptbx.b_as_u(b), occupancy=o))
    xrs.scattering_type_registry(custom_dict={"CONST": self._const})
    return xrs

  def delta_f_mask(self, i_seqs, xyz_old, xyz_new, occ_old, occ_new,
                   exclude_i_seqs):
    """ k_mask (F_mask(new) - F_mask(old)) for moving the atoms i_seqs, from
    the mask recomputed on a local grid (old and new atoms present only
    where their occupancy is > 0; exclude_i_seqs: the residue's atoms, which
    are not in the environment). Returns None if nothing changes. """
    ctx = self.mask_context(i_seqs, [xyz_old, xyz_new], occ_old,
      exclude_i_seqs)
    if(ctx is None): return None
    pos, w = self.mask_cells(ctx, self.mask_solvent(ctx, xyz_old, occ_old),
      self.mask_solvent(ctx, xyz_new, occ_new))
    if(len(w) == 0): return None
    f = self.fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure=self._const_structure(pos, w),
      algorithm="direct").f_calc().data()
    return f * flex.complex_double(self.k_mask, flex.double(
      self.k_mask.size(), 0))

  def mask_gradients(self, pos_list):
    """ dT/d(volume) of solvent pseudo-atoms (constant scattering, B of the
    cell size, scaled by k_mask) at the positions of each array in
    pos_list, in one FFT gradient calculation """
    sizes = [len(p) for p in pos_list]
    if(sum(sizes) == 0): return [np.zeros(0) for p in pos_list]
    pos = np.concatenate([p for p in pos_list if len(p)])
    xrs = self._const_structure(pos, np.ones(len(pos)))
    xray.set_scatterer_grad_flags(scatterers=xrs.scatterers(), occupancy=True)
    k_mask_work = self.k_mask.select(self.work)
    g = self.g_calc * flex.complex_double(k_mask_work,
      flex.double(k_mask_work.size(), 0))
    sfg = self.fmodel.structure_factor_gradients_w(
      u_iso_refinable_params=None, d_target_d_f_calc=g, xray_structure=xrs,
      n_parameters=0, miller_set=self.f_obs_work,
      algorithm=self.fmodel.sfg_params.algorithm,
      extra_params=self.fmodel.sfg_params.extra)
    out = np.array(sfg.d_target_d_occupancy())
    res, k = [], 0
    for n in sizes:
      res.append(out[k:k + n]); k += n
    return res

  def mask_form_factors(self, w):
    """ form factors of solvent pseudo-atoms at the shell centres (signed
    volume x shell-mean k_mask x cell smearing) """
    if(not hasattr(self, "shell_kmask")):
      # mean k_mask by shell: same binning as the curvatures
      km = self.k_mask.select(self.work)
      dss = self.f_obs_work.d_star_sq().data()
      order = flex.sort_permutation(dss)
      n = dss.size(); nb = len(self.shell_s2)
      self.shell_kmask = np.array([flex.mean(km.select(
        order[b * n // nb:(b + 1) * n // nb])) for b in range(nb)])
    b = 8 * math.pi ** 2 * self.mask_cell ** 2 / 12
    return np.outer(w, self.shell_kmask * np.exp(-b * self.shell_s2 / 4))

  # --- maps for searching ---------------------------------------------------
  def _map(self, b):
    if(b in self._maps): return self._maps[b]
    ms = self.f_obs_work
    dss = ms.d_star_sq().data()
    g = self._registry.gaussian_not_optional(scattering_type="C")
    f = flex.double([g.at_d_star_sq(x) for x in dss]) * flex.exp(-b * dss / 4)
    coeffs = ms.array(data=self.g_calc * flex.complex_double(f,
      flex.double(f.size(), 0)))
    cg = maptbx.crystal_gridding(
      unit_cell=self.xrs.unit_cell(),
      space_group_info=self.xrs.space_group_info(),
      d_min=ms.d_min(), resolution_factor=0.25,
      symmetry_flags=maptbx.use_space_group_symmetry)
    m = coeffs.fft_map(crystal_gridding=cg).real_map_unpadded()
    self._maps[b] = np.array(m.as_numpy_array())
    return self._maps[b]

  def _interp(self, grid, frac):
    n = np.array(grid.shape)
    x = (frac % 1.0) * n
    i0 = np.floor(x).astype(int)
    t = x - i0
    out = np.zeros(len(frac))
    for dx in (0, 1):
      for dy in (0, 1):
        for dz in (0, 1):
          w = (np.where(dx, t[:, 0], 1 - t[:, 0]) *
               np.where(dy, t[:, 1], 1 - t[:, 1]) *
               np.where(dz, t[:, 2], 1 - t[:, 2]))
          idx = (i0 + [dx, dy, dz]) % n
          out += w * grid[idx[:, 0], idx[:, 1], idx[:, 2]]
    return out

  def search_setup(self, i_seqs):
    """ Per-atom (bin, weight) for the search maps, and the map scale fitted
    to the exact occupancy gradients of these atoms at their sites """
    types, b, occ = self.atom_params(i_seqs)
    b_bin = []
    for x in b:
      b_bin.append(min(B_BINS, key=lambda y: abs(math.log(y / max(x, 1.)))))
    w = []
    gc = self._registry.gaussian_not_optional(scattering_type="C")
    for t in types:
      w.append(self._registry.gaussian_not_optional(
        scattering_type=t).at_d_star_sq(0) / gc.at_d_star_sq(0))
    return np.array(b_bin), np.array(w), np.array(occ)

  def search_linear(self, xyz, setup, scale):
    b_bin, w, occ = setup
    frac = xyz.dot(np.array(self.xrs.unit_cell().fractionalization_matrix()
      ).reshape(3, 3).T)
    out = np.zeros(len(xyz))
    for b in set(b_bin.tolist()):
      sel = b_bin == b
      out[sel] = self._interp(self._map(b), frac[sel])
    return scale * w * out

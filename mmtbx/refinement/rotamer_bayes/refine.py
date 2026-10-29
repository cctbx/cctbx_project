"""
Local refinement of one side chain for rotamer_bayes (data gathering; may be
too expensive for routine use): the moving atoms of a residue (positions,
and one B shift for all of them) are refined by LBFGS against

    -LL(data) + E_geo/2

with the rest of the model, the bulk-solvent mask, the scales and the error
model fixed. -LL is the target of the f_model manager evaluated on
F_model + k (F_moving(new) - F_moving(current)) (moving atoms by direct
summation); E_geo is the geometry restraint target (sum w delta^2, so E/2
is -ln prior up to a constant) of the restraints within local_radius of the
residue. Comparing a candidate rotamer refined in this way with the current
conformer refined in the same way gives the evidence with the side chain's
own adjustments included (not those of its neighbours).
"""
from __future__ import absolute_import, division, print_function
import math, time
import numpy as np
import scitbx.lbfgs
from libtbx import group_args
from cctbx import xray
from cctbx.array_family import flex

B_SCALE = 10.   # LBFGS parameter for the B shift is dB / B_SCALE

class local_refinement(object):

  def __init__(self, ev, geometry, sites_cart, i_seqs, local_radius=10.):
    """ ev: evidence.evidence; geometry: the whole model's geometry
    restraints manager (or None); sites_cart: all sites (flex.vec3_double);
    i_seqs: the moving atoms. """
    self.ev = ev
    self.i_seqs = flex.size_t([int(i) for i in i_seqs])
    self.xrs_cur = ev.xrs.select(self.i_seqs)
    self.f_cur = ev.fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure=self.xrs_cur, algorithm="direct").f_calc().data()
    self.k_c = flex.complex_double(ev.k, flex.double(ev.k.size(), 0))
    self.k_work = flex.complex_double(ev.k.select(ev.work),
      flex.double(ev.work.count(True), 0))
    self.u_cur = self.xrs_cur.extract_u_iso_or_u_equiv()
    self.geometry = None
    if(geometry is not None):
      xyz = np.array(sites_cart)
      centre = xyz[np.array(list(self.i_seqs))].mean(axis=0)
      d = np.sqrt(((xyz - centre) ** 2).sum(axis=1))
      local = flex.size_t(np.nonzero(d < local_radius)[0].tolist())
      self.local = local
      self.geometry = geometry.select(iselection=local)
      self.local_sites = sites_cart.select(local)
      pos = dict((int(i), k) for k, i in enumerate(local))
      self.moving_in_local = flex.size_t([pos[int(i)] for i in self.i_seqs])
    self.sfg = xray.structure_factors.gradients(
      miller_set=ev.f_obs_work, cos_sin_table=False)

  def _state(self, xyz, db):
    xrs = self.xrs_cur.deep_copy_scatterers()
    xrs.set_sites_cart(flex.vec3_double([tuple(s) for s in xyz]))
    u = self.u_cur + db / (8 * math.pi ** 2)
    u.set_selected(u < 1. / (8 * math.pi ** 2), 1. / (8 * math.pi ** 2))
    xrs.set_u_iso(values=u)
    return xrs

  def data(self, xyz, db, compute_gradients=True):
    """ -LL (relative to the target's zero) and its gradients wrt the moving
    sites (n, 3) and db """
    ev = self.ev
    xrs = self._state(xyz, db)
    f = ev.fmodel.f_obs().structure_factors_from_scatterers(
      xray_structure=xrs, algorithm="direct").f_calc().data()
    fm = ev.f_model.data() + (f - self.f_cur) * self.k_c
    r = ev.functor.core(f_calc=ev.f_model.customized_copy(data=fm),
      compute_gradients=compute_gradients)
    t = ev.n_work * r.target_work()
    if(not compute_gradients): return t, None, None
    g_calc = r.gradients_work() * self.k_work
    xray.set_scatterer_grad_flags(scatterers=xrs.scatterers(), site=True,
      u_iso=True)
    gr = self.sfg(xray_structure=xrs, u_iso_refinable_params=None,
      miller_set=ev.f_obs_work, d_target_d_f_calc=g_calc,
      n_parameters=xrs.n_parameters(), algorithm="direct")
    p = np.array(gr.packed()).reshape(-1, 4)
    g_xyz = ev.n_work * p[:, :3]   # packed site gradients are Cartesian
    g_db = ev.n_work * p[:, 3].sum() / (8 * math.pi ** 2)
    return t, g_xyz, g_db

  def geo(self, xyz, compute_gradients=True):
    if(self.geometry is None): return 0., np.zeros(xyz.shape)
    s = self.local_sites.deep_copy()
    s.set_selected(self.moving_in_local,
      flex.vec3_double([tuple(x) for x in xyz]))
    e = self.geometry.energies_sites(sites_cart=s,
      compute_gradients=compute_gradients)
    g = np.array(e.gradients.select(self.moving_in_local)) \
      if compute_gradients else None
    return e.target / 2, (g / 2 if g is not None else None)

  def set_well(self, residue_xyz, moving, chi_atoms, centres,
               half_width=30., sigma=5., periods=None):
    """ Keep the chi angles within half_width of centres (degrees) with a
    flat-bottomed harmonic restraint (sigma outside): the refinement then
    stays in one rotamer well. residue_xyz: all atoms of the residue;
    moving: mask of the moving atoms (in the order of i_seqs); chi_atoms:
    residue atom indices of each chi. None: no restraint. """
    if(centres is None):
      self.well = None; return
    self.well = group_args(xyz=np.array(residue_xyz, dtype=float),
      moving=np.array(moving), quads=chi_atoms, centres=list(centres),
      periods=periods if periods is not None else [360.] * len(centres),
      half_width=half_width, weight=1. / sigma ** 2)

  def well_penalty(self, xyz):
    """ restraint target and gradients wrt the moving atoms """
    g = np.zeros(xyz.shape)
    w = getattr(self, "well", None)
    if(w is None): return 0., g
    from mmtbx.refinement.rotamer_bayes.side_chain import dihedral, \
      dihedral_gradients
    full = w.xyz.copy(); full[w.moving] = xyz
    index = -np.ones(len(full), dtype=int)
    index[np.nonzero(w.moving)[0]] = np.arange(w.moving.sum())
    e = 0.
    for quad, c, per in zip(w.quads, w.centres, w.periods):
      d = (dihedral(full[quad]) - c + per / 2) % per - per / 2
      excess = max(abs(d) - w.half_width, 0.) * (1 if d > 0 else -1)
      if(excess == 0.): continue
      e += w.weight * excess ** 2
      dg = 2 * w.weight * excess * dihedral_gradients(full[quad])
      for k, a in enumerate(quad):
        if(index[a] >= 0): g[index[a]] += dg[k]
    return e, g

  def run(self, xyz_start, max_iterations=25, refine_b=True):
    t0 = time.time()
    n = len(xyz_start)
    self._refine_b = refine_b
    self._n = n
    self.x = flex.double(list(np.array(xyz_start).ravel()) + [0.])
    self.n_evaluations = 0
    try:
      scitbx.lbfgs.run(target_evaluator=self,
        termination_params=scitbx.lbfgs.termination_parameters(
          max_iterations=max_iterations),
        exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
          ignore_line_search_failed_step_at_upper_bound=True,
          ignore_line_search_failed_maxfev=True,
          ignore_line_search_failed_xtol=True,
          ignore_search_direction_not_descent=True))
    except RuntimeError:
      pass
    xyz, db = self._unpack(self.x)
    ll, _, _ = self.data(xyz, db, compute_gradients=False)
    geo, _ = self.geo(xyz, compute_gradients=False)
    well, _ = self.well_penalty(xyz)
    return group_args(minus_ll=ll, geo=geo, well=well, xyz=xyz, db=db,
      shift_rms=float(np.sqrt(((xyz - np.array(xyz_start)) ** 2).sum(
        axis=1).mean())), n_evaluations=self.n_evaluations,
      time=time.time() - t0)

  def _unpack(self, x):
    a = np.array(x)
    xyz = a[:3 * self._n].reshape(-1, 3)
    db = a[3 * self._n] * B_SCALE if self._refine_b else 0.
    return xyz, db

  def compute_functional_and_gradients(self):
    self.n_evaluations += 1
    xyz, db = self._unpack(self.x)
    t, g_xyz, g_db = self.data(xyz, db)
    e, g_e = self.geo(xyz)
    ew, g_w = self.well_penalty(xyz)
    e += ew
    g = list((g_xyz + g_e + g_w).ravel())
    g.append(g_db * B_SCALE if self._refine_b else 0.)
    return t + e, flex.double(g)

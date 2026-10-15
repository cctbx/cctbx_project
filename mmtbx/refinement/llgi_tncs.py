from __future__ import absolute_import, division, print_function
""" Translational NCS (tNCS) for the LLGI target.

Copies of a molecule related by tNCS (the same molecule in nearly the same
orientation, after the symmetry operator that brings each copy closest to a
pure translation of the group's reference copy) modulate the intensities.
With the copies modelled explicitly, the model's expected intensity is
modulated by

  M(h) = 1 + A(h),
  A(h) = 2 sum_g sum_{m<n in g} sqrt(f_m(s) f_n(s)) c_mn(s) Re Pbar_mn(h)

  Pbar_mn(h) = (1/Nsym) sum_k sum_j w_j exp(2 pi i (T_k^T h).(x_jm - x_jn))
               / sum_j w_j,    w_j = Z_j^2 exp(-(B_jm + B_jn) s^2/4)

(f_m(s): copy m's fraction of the model's scattering at resolution s,
sum_j Z_j^2 exp(-B_j s^2/2) over its atoms over the same sum for all
atoms, so that copies with different B-factors weigh differently at
different resolutions; c_mn(s) = sum_j w_j / sqrt(sum_j Z_j^2 exp(-B_jm
s^2/2) sum_j Z_j^2 exp(-B_jn s^2/2)) <= 1, the loss of coherence from
B-factor differences between corresponding atoms; sums over j are over
the atoms common to all copies of the group, in the frame of the group's
reference copy; T_k the rotations of the space group's representative
operators, Nsym of them, without the inversion: with a centre of symmetry,
R and -R give the same real part). Pbar is the structure factor
of pseudo-atoms at the difference vectors, with constant scattering Z^2
and B = B_m + B_n, in the symmorphic group with the same rotations, so one
FFT structure-factor calculation per pair gives it (pair_modulation).

The true structure's expected intensity is modulated by T(h) = 1 +
Gamma(s) A_cons(h): A_cons is the pair modulation of the group consensus
(the model's copies superposed, averaged and placed back by each copy's
operator, so that the model's own differences between copies do not
weaken it), and Gamma(s) = gamma exp(-2 pi^2 v^2 s^2) is the coherence
of the true copies relative to the consensus (gamma <= 1; v the rms
difference per direction; one curve shared by all groups). The error
variance is then modulated by u(h) = (T - sigmaA^2 M)/(1 - sigmaA^2), and
the expected intensity by t(h) = sigmaA^2 M(h) + (1 - sigmaA^2) u(h) = T,
sigmaA being the coefficient of the unit-normalised model E-value. With
coherence = tau, u(h) = 1 + tau(s) A(h) on the model's own pair terms
instead. The LLGI functors take u and t per reflection (cctbx/xray/targets/llgi.h, tncs_frame); the data are normalised
by t (renormalise). Parameter fits use the u-frame (t = u), in which the
frame is valid for any sigmaA and which, with the exact likelihood, gives
the same result as any other frame.
"""
import math
import itertools
from libtbx import group_args
from cctbx import crystal, miller, sgtbx, xray
from cctbx.array_family import flex
from cctbx.eltbx import xray_scattering
from cctbx.xray import ext as xray_ext
from scitbx.math import superpose
import iotbx.pdb
import iotbx.phil

def _polymer_atoms(chain):
  """ {(resseq, icode, atom name): (i_seq, resname)} for the first
  conformer of the amino-acid and nucleotide residues of a chain, without
  hydrogens """
  out = {}
  for rg in chain.residue_groups():
    ag = rg.atom_groups()[0]
    cls = iotbx.pdb.common_residue_names_get_class(ag.resname)
    if cls not in ("common_amino_acid", "common_rna_dna"):
      continue
    for a in ag.atoms():
      if a.element.strip().upper() in ("H", "D"):
        continue
      out.setdefault((rg.resseq_as_int(), rg.icode, a.name.strip()),
                     (a.i_seq, ag.resname))
  return out

def apply_op(op, sites_frac):
  """ rt_mx applied to fractional coordinates """
  return op.r().as_double() * sites_frac + op.t().as_double()

def _rotation_angle_after_op(sites_frac_a, sites_frac_b, op, unit_cell):
  """ rotation angle (degrees) of the superposition of copy A onto op(B) """
  b = apply_op(op, sites_frac_b)
  lsq = superpose.least_squares_fit(
    reference_sites=unit_cell.orthogonalize(b),
    other_sites=unit_cell.orthogonalize(sites_frac_a))
  return lsq.r.rotation_angle()

def _best_op(sites_frac_a, sites_frac_b, crystal_symmetry):
  """ symmetry operator bringing copy B closest to a pure translation of
  copy A: (rotation angle, rt_mx) """
  uc = crystal_symmetry.unit_cell()
  best = None
  for op in crystal_symmetry.space_group().all_ops():
    angle = _rotation_angle_after_op(sites_frac_a, sites_frac_b, op, uc)
    if best is None or angle < best[0]:
      best = (angle, op)
  return best

class tncs_groups(object):
  """ Groups of chains related by tNCS in a model.

  Chains are copies of the same molecule if at least min_identity of the
  residue names agree at their shared residue numbers, which cover at
  least min_overlap of the shorter chain: corresponding atoms are matched
  by residue number, insertion code and atom name, so copies must be
  numbered alike. Two copies are tNCS-related if,
  after the symmetry operator that brings one closest to a pure
  translation of the other, they differ by less than angle_max degrees;
  groups are the connected components of that relation. Within a group,
  every copy is placed by the operator that brings it closest to a pure
  translation of the group's reference copy (the first), and the atoms
  common to all copies are used.

  .groups: list of group_args(chain_ids, i_seqs (one flex.size_t per copy,
  aligned: the common atoms), chain_i_seqs (one flex.size_t per copy: all
  its polymer atoms, for the scattering fractions), ops (rt_mx per copy),
  angles (degrees, relative to the reference copy)); .n_atoms: the number
  of atoms in the model, which the i_seqs index.
  """
  def __init__(self, pdb_hierarchy, crystal_symmetry, angle_max=10.,
               min_identity=0.9, min_overlap=0.5, min_atoms=50):
    self.crystal_symmetry = crystal_symmetry
    self.n_atoms = pdb_hierarchy.atoms().size()
    sites_frac = crystal_symmetry.unit_cell().fractionalize(
      pdb_hierarchy.atoms().extract_xyz())
    chains = []
    for ch in pdb_hierarchy.models()[0].chains():
      atoms = _polymer_atoms(ch)
      if len(atoms) >= min_atoms:
        chains.append(group_args(id=ch.id, atoms=atoms))
    n = len(chains)
    def ca_keys(c):
      return set(k for k in c.atoms if k[2] in ("CA", "P"))
    parent = list(range(n))
    def root(i):
      while parent[i] != i: i = parent[i]
      return i
    for i, j in itertools.combinations(range(n), 2):
      a, b = chains[i], chains[j]
      ka, kb = ca_keys(a), ca_keys(b)
      shared = sorted(ka & kb)
      if len(shared) < max(3, min_overlap * min(len(ka), len(kb))):
        continue
      same = [k for k in shared if a.atoms[k][1] == b.atoms[k][1]]
      if len(same) < min_identity * len(shared):
        continue
      angle, op = _best_op(
        sites_frac.select(flex.size_t([a.atoms[k][0] for k in same])),
        sites_frac.select(flex.size_t([b.atoms[k][0] for k in same])),
        crystal_symmetry)
      if angle < angle_max:
        parent[root(i)] = root(j)
    components = {}
    for i in range(n):
      components.setdefault(root(i), []).append(i)
    self.groups = []
    for members in sorted(components.values()):
      if len(members) < 2: continue
      common = set(chains[members[0]].atoms)
      for i in members[1:]:
        common &= set(chains[i].atoms)
      ref = chains[members[0]]
      common = sorted(k for k in common
        if all(chains[i].atoms[k][1] == ref.atoms[k][1] for i in members))
      i_seqs = [flex.size_t([chains[i].atoms[k][0] for k in common])
                for i in members]
      ops, angles = [sgtbx.rt_mx()], [0.]
      ref_sites = sites_frac.select(i_seqs[0])
      for sel in i_seqs[1:]:
        angle, op = _best_op(ref_sites, sites_frac.select(sel),
                             crystal_symmetry)
        ops.append(op); angles.append(angle)
      self.groups.append(group_args(
        chain_ids=[chains[i].id for i in members], i_seqs=i_seqs,
        chain_i_seqs=[flex.size_t(sorted(v[0] for v in chains[i].atoms.values()))
                      for i in members],
        ops=ops, angles=angles))

  def __len__(self):
    return len(self.groups)

  def show(self, log):
    if not self.groups:
      print("tNCS: no groups of chains related by tNCS", file=log)
      return
    for g in self.groups:
      print("tNCS group: chains %s, %d common atoms" % (
        " ".join(g.chain_ids), g.i_seqs[0].size()), file=log)
      for cid, op, angle in zip(g.chain_ids[1:], g.ops[1:], g.angles[1:]):
        print("  %s -> %s by %s: rotation %.2f degrees" % (
          cid, g.chain_ids[0], op.as_xyz(), angle), file=log)

def _atomic_weight_sq(xray_structure):
  """ (f0 at s = 0)^2 * occupancy per scatterer """
  reg = xray_structure.scattering_type_registry()
  z = flex.double([reg.gaussian_not_optional(sc.scattering_type).at_d_star_sq(0)
                   for sc in xray_structure.scatterers()])
  return z * z * xray_structure.scatterers().extract_occupancies()

def symmorphic_group(space_group):
  """ the space group with the same rotations (representative operators,
  no centring) and no translational parts """
  result = sgtbx.space_group()
  for op in space_group.smx():
    result.expand_smx(sgtbx.rt_mx(op.r(), sgtbx.tr_vec((0, 0, 0), op.t().den())))
  return result

def _pseudo_atom_structure(crystal_symmetry, sites_frac, z_sq, b):
  sps = crystal.special_position_settings(crystal_symmetry,
    min_distance_sym_equiv=0)
  scatterers = flex.xray_scatterer()
  for site, occ, bb in zip(sites_frac, z_sq, b):
    scatterers.append(xray.scatterer(label="Q", scattering_type="Q",
      site=site, b=bb, occupancy=occ))
  xrs = xray.structure(special_position_settings=sps, scatterers=scatterers)
  xrs.scattering_type_registry(
    custom_dict={"Q": xray_scattering.gaussian(1.0)})
  return xrs

def _weight_sums(d_star_sq, z_sq, b, n_grid=400):
  """ sum_j z_sq_j exp(-b_j s^2/4), tabulated on a grid of s^2 and
  interpolated linearly """
  import numpy as np
  ss_max = flex.max(d_star_sq) * 1.0001 if d_star_sq.size() else 1.
  grid = np.linspace(0, ss_max, n_grid)
  tab = [flex.sum(z_sq * flex.exp(-b * (g / 4.))) for g in grid]
  return flex.double(np.interp(np.array(d_star_sq), grid, tab))

def pair_term(miller_set, sites_frac_m, sites_frac_n, z_sq, b_sum):
  """ Re Pbar_mn(h) on miller_set's indices, by one FFT structure-factor
  calculation in the symmorphic group """
  sym = crystal.symmetry(unit_cell=miller_set.unit_cell(),
    space_group=symmorphic_group(miller_set.space_group()))
  d = sites_frac_m - sites_frac_n
  d = d - flex.vec3_double([tuple(math.floor(x) for x in v) for v in d])
  xrs = _pseudo_atom_structure(sym, d, z_sq, b_sum)
  ms = miller.set(sym, miller_set.indices(), anomalous_flag=False)
  # quality_factor 1e5: errors in Pbar ~0.001 rms (the default, ~0.02)
  fc = ms.structure_factors_from_scatterers(xray_structure=xrs,
    algorithm="fft", quality_factor=100000).f_calc().data()
  nsym = miller_set.space_group().n_smx()   # the operators in symmorphic_group
  norm = _weight_sums(miller_set.d_star_sq().data(), z_sq, b_sum)
  return flex.real(fc) / (nsym * norm)

def consensus_sites(group, sites_frac, unit_cell):
  """ The group's consensus copy placed on each copy (fractional
  coordinates of the common atoms, per copy): the copies (after their
  symmetry operators) superposed on the reference copy, averaged, and
  placed back by each copy's least-squares operator. """
  placed = [apply_op(op, sites_frac.select(sel))
            for op, sel in zip(group.ops, group.i_seqs)]
  xyz = [unit_cell.orthogonalize(p) for p in placed]
  fits = [superpose.least_squares_fit(reference_sites=xyz[0], other_sites=x)
          for x in xyz]
  mean = flex.vec3_double(xyz[0].size(), (0, 0, 0))
  for f, x in zip(fits, xyz):
    mean += f.r.elems * x + f.t.elems
  mean = mean * (1. / len(xyz))
  return [unit_cell.fractionalize(f.r.inverse().elems * (mean - f.t.elems))
          for f in fits]

def pair_modulation(miller_set, groups, xray_structure, basis="model"):
  """ A(h) for the current coordinates and B-factors: the model's
  expected-intensity modulation is M = 1 + A (basis "model"); with basis
  "consensus", the same for the group consensus (consensus_sites), the
  pair modulation of the true structure's expected intensity. Returns
  group_args(a (flex.double on miller_set's indices), fractions (per
  group, the scattering fraction of each copy at s = 0), f_total (the
  largest fraction of the scattering in tNCS groups at any
  resolution)). """
  result = pair_modulations(miller_set, groups, xray_structure, (basis,))
  return group_args(a=result.a[basis], fractions=result.fractions,
    f_total=result.f_total)

def pair_modulations(miller_set, groups, xray_structure,
      bases=("model", "consensus")):
  """ pair_modulation for several bases at once (the weights are computed
  once): group_args(a ({basis: flex.double}), fractions, f_total) """
  assert all(basis in ("model", "consensus") for basis in bases)
  assert xray_structure.scatterers().size() == groups.n_atoms, (
    "tNCS groups were found for a model with %d atoms; this structure has %d"
    % (groups.n_atoms, xray_structure.scatterers().size()))
  w_sq = _atomic_weight_sq(xray_structure)
  w_total = flex.sum(w_sq)
  sites_frac = xray_structure.sites_frac()
  b_iso = xray_structure.extract_u_iso_or_u_equiv() * (8 * math.pi ** 2)
  dss = miller_set.d_star_sq().data()
  def weights(sel):  # sum_j Z_j^2 exp(-B_j s^2/2) over sel, per reflection
    return _weight_sums(dss, w_sq.select(sel), 2 * b_iso.select(sel))
  w_total_s = _weight_sums(dss, w_sq, 2 * b_iso)
  a = dict((basis, flex.double(dss.size(), 0)) for basis in bases)
  f_in_groups = flex.double(dss.size(), 0)
  fractions = []
  for g in groups.groups:
    fractions.append([flex.sum(w_sq.select(sel)) / w_total
                      for sel in g.chain_i_seqs])
    f_s = [weights(sel) / w_total_s for sel in g.chain_i_seqs]
    for f in f_s: f_in_groups += f
    w_common = [weights(sel) for sel in g.i_seqs]
    placed = {}
    for basis in bases:
      if(basis == "consensus"):
        placed[basis] = consensus_sites(g, sites_frac, xray_structure.unit_cell())
      else:
        placed[basis] = [apply_op(op, sites_frac.select(sel))
                         for op, sel in zip(g.ops, g.i_seqs)]
    z_sq = w_sq.select(g.i_seqs[0])
    for m, n in itertools.combinations(range(len(g.i_seqs)), 2):
      b_sum = b_iso.select(g.i_seqs[m]) + b_iso.select(g.i_seqs[n])
      weight = 2 * flex.sqrt(f_s[m] * f_s[n]) * _weight_sums(dss, z_sq, b_sum) \
        / flex.sqrt(w_common[m] * w_common[n])   # including the coherence
      for basis in bases:
        a[basis] += weight * pair_term(miller_set, placed[basis][m],
          placed[basis][n], z_sq, b_sum)
  return group_args(a=a, fractions=fractions,
    f_total=max(flex.max(f_in_groups), sum(sum(f) for f in fractions)))

# ---------------------------------------------------------------------------
# Data in a tNCS frame, and the fits

def renormalise(i_obs, sig_i_obs, resn, t, centric_flags):
  """ The intensities normalised by RESN^2 t, with moment-matched Dobs and
  Eeff (llgi_hybrid.rice_dobs_eeff) and Feff = Eeff RESN sqrt(t), the
  amplitude on the scale of the data. All arguments flex arrays.
  Returns group_args(e_obs_sq, sig_e_obs_sq, force_exact, dobs, eeff,
  feff). """
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  e_obs_sq, sig_e_obs_sq = llgi_hybrid.e_obs_sq_and_sigma(
    i_obs=i_obs, sig_i_obs=sig_i_obs, resn=resn, teps=t)
  fe = llgi_hybrid.force_exact_flags(e_obs_sq, sig_e_obs_sq, centric_flags)
  dobs, eeff = llgi_hybrid.rice_dobs_eeff(fe.rice)
  return group_args(e_obs_sq=e_obs_sq, sig_e_obs_sq=sig_e_obs_sq,
    force_exact=fe.force_exact, dobs=dobs, eeff=eeff,
    feff=eeff * resn * flex.sqrt(t))

def null_log_z(e_obs_sq, sig_e_obs_sq, centric_flags, c):
  """ M(0, b) of llgi_exact.h for the observations rescaled by c: with
  the LLGI gain in the c-normalised frame, log p(I) = gain + M(0, b) +
  terms that do not depend on c. """
  return xray_ext.llgi_hybrid(e_obs_sq=e_obs_sq / c,
    sig_e_obs_sq=sig_e_obs_sq / c, force_exact=flex.bool(e_obs_sq.size(), False),
    centric_flags=centric_flags, rice_kappa=0.).null_log_z

class linear_in_d_star_sq(object):
  """ p0 + p1 x with x = d*^2/d*^2_max: the form of the tNCS strengths
  rho(s) and tau(s) """
  def __init__(self, d_star_sq_max, p0, p1=0.):
    self.d_star_sq_max = d_star_sq_max
    self.p0, self.p1 = float(p0), float(p1)  # not numpy: flex arithmetic

  def at(self, d_star_sq):
    return self.p0 + self.p1 * (d_star_sq / self.d_star_sq_max)

  def end_values(self):
    return (self.p0, self.p0 + self.p1)

def _fit_strength(f, f_total, x0, linear):
  """ minimise f(p0, p1) over the strength curve p0 + p1 x: p1 = 0 and a
  bounded 1-D search unless linear (Nelder-Mead from x0) """
  if not linear:
    from scipy.optimize import minimize_scalar
    hi = 1. / max(f_total, 1e-6)
    r = minimize_scalar(lambda p0: f(p0, 0.), bounds=(-0.5, hi),
      method="bounded", options=dict(xatol=1e-4))
    return [r.x, 0.]
  from scipy.optimize import minimize
  r = minimize(lambda x: f(x[0], x[1]), x0, method="Nelder-Mead",
    options=dict(xatol=1e-4, fatol=1e-9, maxiter=400, maxfev=400,
                 initial_simplex=[x0, [x0[0] + 0.2, x0[1]],
                                  [x0[0], x0[1] + 0.2]]))
  return list(r.x)

def _strength_ok(curve, f_total):
  """ the modulation 1 + strength*A is positive everywhere, and the
  strength is between -0.5 and the positivity limit 1/f_total """
  lo, hi = curve.end_values()
  if min(lo, hi) < -0.5 or max(lo, hi) > 1. / max(f_total, 1e-6):
    return False
  return True

def fit_wilson_modulation(e_obs_sq, sig_e_obs_sq, centric_flags, a,
      d_star_sq, f_total, selection=None, max_reflections=20000, seed=0,
      linear=False):
  """ Data-only fit of the expected-intensity modulation t = 1 + rho(s) A
  (rho constant, or linear in d*^2): the Wilson likelihood, with
  measurement error, of the intensities normalised by RESN^2. Returns the
  rho curve (linear_in_d_star_sq). """
  sel = _subset(e_obs_sq.size(), selection, max_reflections, seed)
  eo, sg = e_obs_sq.select(sel), sig_e_obs_sq.select(sel)
  cen, aa, ss = centric_flags.select(sel), a.select(sel), d_star_sq.select(sel)
  ss_max = flex.max(d_star_sq)
  def f(p0, p1):
    curve = linear_in_d_star_sq(ss_max, p0, p1)
    if not _strength_ok(curve, f_total): return 1e10
    t = 1 + curve.at(ss) * aa
    if flex.min(t) <= 0.01: return 1e10
    return -flex.mean(null_log_z(eo, sg, cen, t))
  p = _fit_strength(f, f_total, [0.5 / max(f_total, 0.5), 0.], linear)
  return linear_in_d_star_sq(ss_max, p[0], p[1])

def _subset(n, selection, max_reflections, seed):
  if selection is None:
    selection = flex.bool(n, True)
  isel = selection.iselection()
  if isel.size() > max_reflections:
    import random
    rnd = random.Random(seed)
    isel = flex.size_t(sorted(rnd.sample(list(isel), max_reflections)))
  return isel

def fit_tau(i_obs, sig_i_obs, resn, centric_flags, e_model, sigmaa, a,
      d_star_sq, f_total, selection, tau_start=None, linear=False):
  """ tau(s) (constant, or linear in d*^2) of the error-variance modulation u = 1 +
  tau A, by the exact E-scale likelihood over the selected reflections,
  with sigmaA (E-scale, per reflection) fixed. The objective is the
  absolute log-likelihood: the gain in the u-frame (llgi_e functor with
  tncs_u = tncs_t = u, valid for any sigmaA) plus the u-dependent null
  term. Returns the tau curve (linear_in_d_star_sq). """
  isel = selection.iselection()
  io, si, rn = i_obs.select(isel), sig_i_obs.select(isel), resn.select(isel)
  cen, em = centric_flags.select(isel), e_model.select(isel)
  sa, aa, ss = sigmaa.select(isel), a.select(isel), d_star_sq.select(isel)
  ss_max = flex.max(d_star_sq)
  n = isel.size()
  def f(p0, p1):
    curve = linear_in_d_star_sq(ss_max, p0, p1)
    if not _strength_ok(curve, f_total): return 1e10
    u = 1 + curve.at(ss) * aa
    if flex.min(u) <= 0.01: return 1e10
    d = renormalise(io, si, rn, u, cen)
    hybrid = xray_ext.llgi_hybrid(e_obs_sq=d.e_obs_sq,
      sig_e_obs_sq=d.sig_e_obs_sq, force_exact=d.force_exact,
      centric_flags=cen, rice_kappa=0.)
    r = xray_ext.llgi_e_sigmaa_target_and_gradients(e_eff=d.eeff,
      selection=flex.bool(n, True), e_model=em, dobs=d.dobs, sigmaa=sa,
      centric_flags=cen, hybrid=hybrid, tncs_u=u, tncs_t=u)
    return r.target() - flex.mean(hybrid.null_log_z)
  x0 = [0.5, 0.] if tau_start is None else list(tau_start)
  p = _fit_strength(f, f_total, x0, linear)
  return linear_in_d_star_sq(ss_max, p[0], p[1])

class coherence_curve(object):
  """ Gamma(s) = gamma exp(-2 pi^2 v^2 s^2): the coherence of the true
  tNCS-related copies relative to the model's consensus, gamma its
  low-resolution limit and v (A) the rms difference per direction """
  def __init__(self, gamma, v):
    self.gamma, self.v = float(gamma), abs(float(v))

  def at(self, d_star_sq):
    return self.gamma * flex.exp(
      -2 * math.pi ** 2 * self.v ** 2 * d_star_sq)

  def summary(self):
    return "coherence gamma = %.3f, rms difference %.2f A" % (
      self.gamma, self.v)

# u is kept at least U_MIN where the model's own modulation exceeds the
# true one by more than the coherence allows (T < sigmaA^2 M); the same
# floor applies to t0 = 1 + rho A and to u = 1 + tau A
U_MIN = 0.02

# fewer free reflections than this: the coherence (or tau) is fitted to all
MIN_FIT_REFLECTIONS = 200

def error_variance_modulation(sigmaa, a_model, a_cons, gamma):
  """ u = (T - sigmaA^2 M)/(1 - sigmaA^2), with T = 1 + Gamma A_cons and
  M = 1 + A_model (Gamma per reflection), at least U_MIN. Returns
  (u, number of reflections raised to U_MIN). """
  sa2 = flex.pow2(sigmaa)
  sa2.set_selected(sa2 > 0.999, 0.999)
  u = (1 + gamma * a_cons - sa2 * (1 + a_model)) / (1 - sa2)
  low = u < U_MIN
  u.set_selected(low, U_MIN)
  return u, low.count(True)

def fit_coherence(e_obs_sq, sig_e_obs_sq, centric_flags, e_model, sigmaa,
      a_model, a_cons, d_star_sq, selection, start=None):
  """ Gamma(s) (coherence_curve) by the exact E-scale likelihood over the
  selected reflections, with sigmaA fixed, as fit_tau but with u from
  error_variance_modulation; 0 <= gamma <= 1, 0 <= v <= 3 A. """
  from scipy.optimize import minimize
  isel = selection.iselection()
  io, si = e_obs_sq.select(isel), sig_e_obs_sq.select(isel)
  cen, em = centric_flags.select(isel), e_model.select(isel)
  sa, am = sigmaa.select(isel), a_model.select(isel)
  ac, ss = a_cons.select(isel), d_star_sq.select(isel)
  n = isel.size()
  one = flex.double(n, 1.)
  def f(x):
    if not (0 <= x[0] <= 1 and 0 <= x[1] <= 3): return 1e10
    u = error_variance_modulation(sa, am, ac,
      coherence_curve(x[0], x[1]).at(ss))[0]
    d = renormalise(io, si, one, u, cen)
    hybrid = xray_ext.llgi_hybrid(e_obs_sq=d.e_obs_sq,
      sig_e_obs_sq=d.sig_e_obs_sq, force_exact=d.force_exact,
      centric_flags=cen, rice_kappa=0.)
    r = xray_ext.llgi_e_sigmaa_target_and_gradients(e_eff=d.eeff,
      selection=flex.bool(n, True), e_model=em, dobs=d.dobs, sigmaa=sa,
      centric_flags=cen, hybrid=hybrid, tncs_u=u, tncs_t=u)
    return r.target() - flex.mean(hybrid.null_log_z)
  x0 = [0.95, 0.2] if start is None else [min(start.gamma, 0.99),
    max(start.v, 0.05)]
  r = minimize(f, x0, method="Nelder-Mead", options=dict(xatol=1e-4,
    fatol=1e-9, maxiter=400, maxfev=400, initial_simplex=[x0,
    [x0[0] - 0.05, x0[1]], [x0[0], x0[1] + 0.1]]))
  return coherence_curve(r.x[0], r.x[1])

def expected_intensity_modulation(sigmaa, a, u):
  """ t = sigmaA^2 M + (1 - sigmaA^2) u, with M = 1 + A: the tNCS
  modulation of the expected intensity. t exceeds (1 - sigmaA^2) u
  wherever sigmaA > 0 and M > 0, so the frame of the LLGI functors is valid
  there, and in that frame the model term's coefficient is the correlation
  sigmaA sqrt(M/t). """
  return sigmaa * sigmaa * (1 + a) + (1 - sigmaa * sigmaa) * u

# ---------------------------------------------------------------------------
# Functor calls with or without tNCS

def e_sigmaa_functor(e_eff, selection, e_model, dobs, sigmaa, centric_flags,
      hybrid=None, tncs_u=None):
  """ cctbx.xray.ext.llgi_e_sigmaa_target_and_gradients, in the u-frame
  (tncs_t = tncs_u) when tncs_u is given """
  args = dict(e_eff=e_eff, selection=selection, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags, hybrid=hybrid)
  if(tncs_u is not None):
    args.update(tncs_u=tncs_u, tncs_t=tncs_u)
  return xray_ext.llgi_e_sigmaa_target_and_gradients(**args)

def sigmaa_scatfrac_functor(tncs_u=None, **args):
  """ cctbx.xray.ext.llgi_sigmaa_scatfrac_target_and_gradients, in the
  u-frame (tncs_t = tncs_u) when tncs_u is given """
  if(tncs_u is not None):
    args.setdefault("hybrid", None)
    args.update(tncs_u=tncs_u, tncs_t=tncs_u)
  return xray_ext.llgi_sigmaa_scatfrac_target_and_gradients(**args)

def frame(a, u, t):
  """ llgi.h tncs_frame for arrays: (a', s), both 0 where the frame is
  invalid (no contribution) """
  r = xray_ext.llgi_tncs_frame(a=a, u=u, t=t)
  return r.a_p, r.s

def arrays(llgi_data):
  """ (u, t) as flex.double from llgi_data, or (None, None) without tNCS """
  if(llgi_data is None or getattr(llgi_data, "tncs_u", None) is None):
    return None, None
  return llgi_data.tncs_u.data(), llgi_data.tncs_t.data()

# ---------------------------------------------------------------------------
# Setting up and updating the tNCS treatment of fmodel's LLGI data

llgi_tncs_params = iotbx.phil.parse("""\
  enabled = Auto
    .type = bool
    .short_caption = Model translational NCS in the LLGI target
    .help = "Account for translational NCS (copies of a molecule in nearly " \\
            "the same orientation) in the LLGI target. Auto: if the model " \\
            "has groups of chains related by tNCS."
  angle_max = 10.
    .type = float
    .short_caption = Largest rotation between tNCS-related copies (degrees)
  ncs_restraints = False
    .type = bool
    .short_caption = NCS restraints for tNCS-related copies
    .help = "With tNCS groups in the model, switch on NCS restraints " \\
            "(pdb_interpretation.ncs_search.enabled) in phenix.refine. " \\
            "They pull the model's copies together, which costs accuracy " \\
            "where the true copies differ."
  coherence = *physical tau
    .type = choice
    .short_caption = Model of the coherence between tNCS-related copies
    .help = "physical: the true copies' coherence relative to the model's " \\
            "consensus, gamma exp(-2 pi^2 v^2 s^2), on the consensus pair " \\
            "terms, with the error-variance modulation derived from it. " \\
            "tau: the error-variance modulation 1 + tau A on the model's " \\
            "own pair terms (constant tau, or linear in d*^2)."
  tau_linear_in_resolution = False
    .type = bool
    .short_caption = Let tau vary linearly with d*^2 (coherence = tau only)
  n_cycles = 2
    .type = int
    .short_caption = sigmaA/coherence fit cycles per update
  fit_reflections = all *free
    .type = choice
    .short_caption = Reflections for the coherence (or tau) fit
    .help = "free: the R-free set, as for sigmaA; all: every reflection " \\
            "(more precise, but refinement fits the working set, which can " \\
            "inflate the coherence)."
""")

def renormalise_e(e_obs_sq_base, sig_e_obs_sq_base, resn, t, centric_flags):
  """ renormalise for observations already on the E scale without tNCS
  (I/RESN^2): Feff = Eeff RESN sqrt(t) """
  one = flex.double(e_obs_sq_base.size(), 1.)
  d = renormalise(e_obs_sq_base, sig_e_obs_sq_base, one, t, centric_flags)
  d.feff = d.feff * resn
  return d

def _hybrid(d, centric_flags, rice_kappa):
  return xray_ext.llgi_hybrid(e_obs_sq=d.e_obs_sq,
    sig_e_obs_sq=d.sig_e_obs_sq, force_exact=d.force_exact,
    centric_flags=centric_flags, rice_kappa=rice_kappa)

def _framed_llgi_data(llgi_data, f_obs, d, u, t, **changes):
  """ llgi_data with the data d (renormalise_e) in the t-frame and the
  tNCS arrays; the cached hybrid is dropped """
  return replace_llgi_data(llgi_data,
    teps=f_obs.array(data=flex.double(f_obs.size(), 1.)),
    dobs=f_obs.array(data=d.dobs), feff=f_obs.array(data=d.feff),
    e_obs_sq=f_obs.array(data=d.e_obs_sq),
    sig_e_obs_sq=f_obs.array(data=d.sig_e_obs_sq),
    force_exact=f_obs.array(data=d.force_exact),
    tncs_u=f_obs.array(data=u),
    tncs_t=f_obs.array(data=t), _hybrid_object=None, **changes)

def replace_llgi_data(llgi_data, **changes):
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  return llgi_hybrid.replace_llgi_data(llgi_data, **changes)

def setup(llgi_data, f_obs, pdb_hierarchy, xray_structure, params=None,
      log=None):
  """ The tNCS treatment for llgi_data (on f_obs's index set), if the model
  has groups of chains related by tNCS: the observations are kept on the E
  scale without tNCS (nacelle's TEPS is not used), and renormalised by
  t0 = 1 + rho A, rho fitted to the data alone, until sigmaA and the
  coherence (or tau) are first fitted (fit_and_renormalise). Returns the
  new llgi_data (Feff in it is the new f_obs), or None if there are no tNCS
  groups or the treatment is switched off. """
  if(params is None):
    params = llgi_tncs_params.extract()
  if(params.enabled is False):
    return None
  assert getattr(llgi_data, "tncs_state", None) is None  # t would apply twice
  if(getattr(llgi_data, "e_obs_sq", None) is None):
    from libtbx.utils import Sorry
    raise Sorry("The LLGI treatment of translational NCS needs the "
      "intensities (I, SIGI); the LLGI data have none.")
  groups = tncs_groups(pdb_hierarchy, f_obs.crystal_symmetry(),
    angle_max=params.angle_max)
  if(log is not None): groups.show(log)
  if(len(groups) == 0):
    return None
  centric = f_obs.centric_flags().data()
  dss = f_obs.d_star_sq().data()
  teps = llgi_data.teps.data()
  e_base = llgi_data.e_obs_sq.data() * teps
  sig_base = llgi_data.sig_e_obs_sq.data() * teps
  mod = pair_modulation(f_obs, groups, xray_structure)
  rho = fit_wilson_modulation(e_base, sig_base, centric, mod.a, dss,
    mod.f_total)
  t0 = 1 + rho.at(dss) * mod.a
  t0.set_selected(t0 < U_MIN, U_MIN)   # rho was checked on a subset only
  d = renormalise_e(e_base, sig_base, llgi_data.resn.data(), t0, centric)
  if(log is not None):
    print("tNCS: %.3f of the model's scattering in tNCS groups; data-only "
      "modulation rho = %.3f%s" % (mod.f_total, rho.p0,
      "" if rho.p1 == 0 else " + %.3f d*^2/d*^2_max" % rho.p1), file=log)
  return _framed_llgi_data(llgi_data, f_obs, d, u=t0, t=t0,
    e_obs_sq_base=f_obs.array(data=e_base),
    sig_e_obs_sq_base=f_obs.array(data=sig_base),
    tncs_state=group_args(groups=groups, tau=rho,
      coherence=coherence_curve(min(max(rho.p0, 0.), 0.99), 0.2),
      f_total=mod.f_total, params=params))

def _rice_kappa(llgi_data):
  return llgi_data.hybrid_params.rice_kappa \
    if getattr(llgi_data, "hybrid_params", None) is not None else 0.1

def _fit_scatfrac(fmodel, d, u, sigmaa, scatfrac_params, rice_kappa,
      warnings):
  """ ScatFrac (F scale, working set) in the u-frame, for the data d
  (renormalise_e by u) """
  import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
  if(scatfrac_params is None):
    scatfrac_params = llgi_scatfrac.llgi_scatfrac_params.extract()
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  result = llgi_scatfrac.estimate_llgi_scatfrac(
    f_eff=d.feff, working_selection=~fmodel.r_free_flags().data(),
    f_calc=fmodel.f_model().data(), dobs=d.dobs, sigmaa=sigmaa,
    teps=flex.double(f_obs.size(), 1.), resn=fmodel.llgi_data().resn.data(),
    centric_flags=centric, d_star_sq=f_obs.d_star_sq().data(),
    scale_factor=fmodel.scale_ml_wrapper(), params=scatfrac_params,
    hybrid=_hybrid(d, centric, rice_kappa), tncs_u=u)
  if(result.lbfgs_error is not None):
    warnings.append("ScatFrac fit: L-BFGS stopped early: %s" %
      result.lbfgs_error)
  return result

def refit_sigmaa_scatfrac(fmodel, scatfrac_params=None, e_params=None):
  """ sigmaA and ScatFrac for fmodel's current model in the current tNCS
  frame (u and t, and so the data, are kept; A(h) is not recomputed): for
  refits on a copy of fmodel whose model is changed temporarily, e.g. with
  atoms removed or B-factors shaken (mmtbx.refinement.weights.
  adp_gradients). Returns a group_args as fit_and_renormalise. """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  llgi_data = fmodel.llgi_data()
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  resn = llgi_data.resn.data()
  rice_kappa = _rice_kappa(llgi_data)
  u = llgi_data.tncs_u.data()
  d = renormalise_e(llgi_data.e_obs_sq_base.data(),
    llgi_data.sig_e_obs_sq_base.data(), resn, u, centric)
  warnings = []
  sigmaa_result = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel(fmodel,
    dobs=d.dobs, feff=d.eeff * resn, resn=resn, params=e_params,
    hybrid=_hybrid(d, centric, rice_kappa), tncs_u=u)
  if(sigmaa_result.lbfgs_error is not None):
    warnings.append("sigmaA fit: L-BFGS stopped early: %s" %
      sigmaa_result.lbfgs_error)
  sigmaa = sigmaa_result.sigmaa
  scatfrac_result = _fit_scatfrac(fmodel, d, u, sigmaa, scatfrac_params,
    rice_kappa, warnings)
  fmodel.set_llgi_data(replace_llgi_data(llgi_data,
    sigmaa=f_obs.array(data=sigmaa), sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=f_obs.array(data=scatfrac_result.scatfrac),
    _sigmaa_f_model=fmodel.f_model().data().deep_copy()))
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac_result.scatfrac, target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac, warnings=warnings, tau=None,
    coherence=None, tncs_summary=None)

def _current_groups(state, crystal_symmetry, xray_structure, pdb_hierarchy,
      warnings):
  """ state.groups, or, if the model has gained or lost atoms since they were
  found (e.g. waters, ions), the groups found again in pdb_hierarchy, which
  must match xray_structure """
  from libtbx.utils import Sorry
  n_atoms = xray_structure.scatterers().size()
  if(state.groups.n_atoms == n_atoms):
    return state.groups
  if(pdb_hierarchy is None or pdb_hierarchy.atoms().size() != n_atoms):
    raise Sorry("tNCS: the model now has %d atoms, not the %d for which the "
      "tNCS groups were found, and no model with the current atoms was "
      "given to find them again." % (n_atoms, state.groups.n_atoms))
  groups = tncs_groups(pdb_hierarchy, crystal_symmetry,
    angle_max=state.params.angle_max)
  if(len(groups) == 0):
    raise Sorry("tNCS: no groups of chains related by tNCS in the model with "
      "%d atoms (they were found in the model with %d atoms)." % (n_atoms,
      state.groups.n_atoms))
  warnings.append("tNCS: groups found again for the model's %d atoms "
    "(previously %d)" % (n_atoms, state.groups.n_atoms))
  return groups

def fit_and_renormalise(fmodel, scatfrac_params=None, e_params=None,
      pdb_hierarchy=None):
  """ Per macrocycle, for fmodel's current model and scales: A(h) (and
  A_cons(h) for coherence = physical); sigmaA (E scale, R-free set) and
  the coherence or tau (params.fit_reflections), alternately, in the
  u-frame; ScatFrac (F scale, working set) in the u-frame; then the data
  renormalised by t = sigmaA^2 (1 + A) + (1 - sigmaA^2) u, with Feff the
  new f_obs. Returns a group_args like llgi_e_sigmaa.
  estimate_sigmaa_e_then_scatfrac_f's, plus .tau or .coherence and
  .tncs_summary. pdb_hierarchy: the model, matching fmodel's
  xray_structure, in which the tNCS groups are found again if the number of
  atoms has changed. """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  llgi_data = fmodel.llgi_data()
  state = llgi_data.tncs_state
  params = state.params
  warnings = []
  groups = _current_groups(state, fmodel.f_obs().crystal_symmetry(),
    fmodel.xray_structure, pdb_hierarchy, warnings)
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  dss = f_obs.d_star_sq().data()
  n = f_obs.size()
  resn = llgi_data.resn.data()
  e_base = llgi_data.e_obs_sq_base.data()
  sig_base = llgi_data.sig_e_obs_sq_base.data()
  rice_kappa = _rice_kappa(llgi_data)
  physical = params.coherence == "physical"
  pm = pair_modulations(f_obs, groups, fmodel.xray_structure,
    ("model", "consensus") if physical else ("model",))
  a, a_cons = pm.a["model"], pm.a.get("consensus")
  tau, coherence = state.tau, state.coherence
  f_total = state.f_total if groups is state.groups else pm.f_total
  fit_selection = flex.bool(n, True)
  if(params.fit_reflections == "free"):
    free = fmodel.r_free_flags().data()
    if(free.count(True) >= MIN_FIT_REFLECTIONS):
      fit_selection = free
    else:
      warnings.append("tNCS: only %d free reflections; the coherence is "
        "fitted to all reflections" % free.count(True))
  n_low = 0
  u = llgi_data.tncs_u.data()     # the current frame, for the first sigmaA
  for cycle in range(max(1, params.n_cycles)):
    d = renormalise_e(e_base, sig_base, resn, u, centric)
    hybrid_u = _hybrid(d, centric, rice_kappa)
    sigmaa_result = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel(fmodel,
      dobs=d.dobs, feff=d.eeff * resn, resn=resn, params=e_params,
      hybrid=hybrid_u, tncs_u=u)
    sigmaa = sigmaa_result.sigmaa
    e_model = sigmaa_result.e_model
    if(physical):
      coherence = fit_coherence(e_base, sig_base, centric, e_model, sigmaa,
        a, a_cons, dss, fit_selection, start=coherence)
      u, n_low = error_variance_modulation(sigmaa, a, a_cons,
        coherence.at(dss))
    else:
      tau = fit_tau(e_base, sig_base, flex.double(n, 1.), centric, e_model,
        sigmaa, a, dss, f_total, fit_selection,
        tau_start=[tau.p0, tau.p1], linear=params.tau_linear_in_resolution)
      u = 1 + tau.at(dss) * a   # tau was checked on fit_selection only
      n_low = (u < U_MIN).count(True)
      u.set_selected(u < U_MIN, U_MIN)
  if(sigmaa_result.lbfgs_error is not None):
    warnings.append("sigmaA fit: L-BFGS stopped early: %s" %
      sigmaa_result.lbfgs_error)
  if(n_low > 0):
    warnings.append("tNCS: the model's modulation exceeds the true one for "
      "%d reflections (error-variance modulation held at %g)" % (n_low,
      U_MIN))
  d = renormalise_e(e_base, sig_base, resn, u, centric)
  scatfrac_result = _fit_scatfrac(fmodel, d, u, sigmaa, scatfrac_params,
    rice_kappa, warnings)
  scatfrac = scatfrac_result.scatfrac
  t = expected_intensity_modulation(sigmaa, a, u)
  d = renormalise_e(e_base, sig_base, resn, t, centric)
  new = _framed_llgi_data(llgi_data, f_obs, d, u=u, t=t,
    sigmaa=f_obs.array(data=sigmaa),
    sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=f_obs.array(data=scatfrac),
    tncs_state=group_args(groups=groups, tau=tau,
      coherence=coherence, f_total=f_total, params=params),
    _sigmaa_f_model=fmodel.f_model().data().deep_copy())
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=f_obs.customized_copy(data=d.feff))
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac, target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac,
    warnings=warnings, tau=None if physical else tau,
    coherence=coherence if physical else None,
    tncs_summary=coherence.summary() if physical else "tau = %.3f%s" % (
      tau.p0, "" if tau.p1 == 0 else " + %.3f d*^2/d*^2_max" % tau.p1))

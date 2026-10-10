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
reference copy; T_k the rotations of the space group). Pbar is the
structure factor
of pseudo-atoms at the difference vectors, with constant scattering Z^2
and B = B_m + B_n, in the symmorphic group with the same rotations, so one
FFT structure-factor calculation per pair gives it (pair_modulation).

The error variance is modulated by u(h) = 1 + tau(s) A(h) (one tau shared
by all groups), and the expected intensity by t(h) = m M(h) + (1 - a^2) u(h),
where a is the target's effective sigmaA and m the fraction of the expected
intensity explained by the model. The LLGI functors take u and t per
reflection (cctbx/xray/targets/llgi.h, tncs_frame); the data are normalised
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
  least min_overlap of the shorter chain. Two copies are tNCS-related if,
  after the symmetry operator that brings one closest to a pure
  translation of the other, they differ by less than angle_max degrees;
  groups are the connected components of that relation. Within a group,
  every copy is placed by the operator that brings it closest to a pure
  translation of the group's reference copy (the first), and the atoms
  common to all copies are used.

  .groups: list of group_args(chain_ids, i_seqs (one flex.size_t per copy,
  aligned: the common atoms), chain_i_seqs (one flex.size_t per copy: all
  its polymer atoms, for the scattering fractions), ops (rt_mx per copy),
  angles (degrees, relative to the reference copy)).
  """
  def __init__(self, pdb_hierarchy, crystal_symmetry, angle_max=10.,
               min_identity=0.9, min_overlap=0.5, min_atoms=50):
    self.crystal_symmetry = crystal_symmetry
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
  nsym = miller_set.space_group().order_p()
  norm = _weight_sums(miller_set.d_star_sq().data(), z_sq, b_sum)
  return flex.real(fc) / (nsym * norm)

def pair_modulation(miller_set, groups, xray_structure):
  """ A(h) for the current coordinates and B-factors: the model's
  expected-intensity modulation is M = 1 + A. Returns group_args(a (flex.
  double on miller_set's indices), fractions (per group, the scattering
  fraction of each copy at s = 0), f_total (the largest fraction of the
  scattering in tNCS groups at any resolution)). """
  w_sq = _atomic_weight_sq(xray_structure)
  w_total = flex.sum(w_sq)
  sites_frac = xray_structure.sites_frac()
  b_iso = xray_structure.extract_u_iso_or_u_equiv() * (8 * math.pi ** 2)
  dss = miller_set.d_star_sq().data()
  def weights(sel):  # sum_j Z_j^2 exp(-B_j s^2/2) over sel, per reflection
    return _weight_sums(dss, w_sq.select(sel), 2 * b_iso.select(sel))
  w_total_s = _weight_sums(dss, w_sq, 2 * b_iso)
  a = flex.double(dss.size(), 0)
  f_in_groups = flex.double(dss.size(), 0)
  fractions = []
  for g in groups.groups:
    fractions.append([flex.sum(w_sq.select(sel)) / w_total
                      for sel in g.chain_i_seqs])
    f_s = [weights(sel) / w_total_s for sel in g.chain_i_seqs]
    for f in f_s: f_in_groups += f
    w_common = [weights(sel) for sel in g.i_seqs]
    placed = [apply_op(op, sites_frac.select(sel))
              for op, sel in zip(g.ops, g.i_seqs)]
    z_sq = w_sq.select(g.i_seqs[0])
    for m, n in itertools.combinations(range(len(g.i_seqs)), 2):
      b_sum = b_iso.select(g.i_seqs[m]) + b_iso.select(g.i_seqs[n])
      coherence = _weight_sums(dss, z_sq, b_sum) / flex.sqrt(
        w_common[m] * w_common[n])
      a += 2 * flex.sqrt(f_s[m] * f_s[n]) * coherence * pair_term(
        miller_set, placed[m], placed[n], z_sq, b_sum)
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

def expected_intensity_modulation(model_fraction, a_eff, a, u):
  """ t = m M + (1 - a_eff^2) u, with M = 1 + A: the tNCS modulation of
  the expected intensity, m the fraction of the expected intensity
  explained by the model (E-scale sigmaA^2) and a_eff the effective sigmaA
  of the target. t exceeds (1 - a_eff^2) u wherever m > 0, so the frame
  of the LLGI functors is valid. """
  return model_fraction * (1 + a) + (1 - a_eff * a_eff) * u

# ---------------------------------------------------------------------------
# Functor calls with or without tNCS

def e_sigmaa_functor(e_eff, selection, e_model, dobs, sigmaa, centric_flags,
      hybrid=None, tncs_u=None, tncs_t=None):
  """ cctbx.xray.ext.llgi_e_sigmaa_target_and_gradients, with the tNCS
  arrays when tncs_u is given (tncs_t defaults to tncs_u: the u-frame) """
  from cctbx.xray import ext as xray_ext
  args = dict(e_eff=e_eff, selection=selection, e_model=e_model, dobs=dobs,
    sigmaa=sigmaa, centric_flags=centric_flags, hybrid=hybrid)
  if(tncs_u is not None):
    args.update(tncs_u=tncs_u, tncs_t=tncs_u if tncs_t is None else tncs_t)
  return xray_ext.llgi_e_sigmaa_target_and_gradients(**args)

def sigmaa_scatfrac_functor(tncs_u=None, tncs_t=None, **args):
  """ cctbx.xray.ext.llgi_sigmaa_scatfrac_target_and_gradients, with the
  tNCS arrays when tncs_u is given (tncs_t defaults to tncs_u) """
  from cctbx.xray import ext as xray_ext
  if(tncs_u is not None):
    args.setdefault("hybrid", None)
    args.update(tncs_u=tncs_u, tncs_t=tncs_u if tncs_t is None else tncs_t)
  return xray_ext.llgi_sigmaa_scatfrac_target_and_gradients(**args)

def frame(a, u, t):
  """ llgi.h tncs_frame for arrays: (a', s), with a' = 0 where the frame
  is invalid (no contribution) """
  a_p = flex.double(a.size(), 0)
  s = flex.double(a.size(), 0)
  a_p_sq = 1 - (1 - a * a) * u / t
  ok = (a > 0) & (a < 1) & (u > 0) & (t > 0) & (a_p_sq > 0)
  a_p.set_selected(ok, flex.sqrt(a_p_sq.select(ok)))
  s.set_selected(ok, a.select(ok) / (flex.sqrt(t.select(ok)) * a_p.select(ok)))
  return a_p, s

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
  tau_linear_in_resolution = False
    .type = bool
    .short_caption = Let tau vary linearly with d*^2
    .help = "tau, the correlation of the model errors in tNCS-related " \\
            "copies, is constant by default."
  n_cycles = 2
    .type = int
    .short_caption = sigmaA/tau fit cycles per update
  tau_reflections = *all free
    .type = choice
    .short_caption = Reflections for the tau fit
    .help = "all: every reflection (precise, but refinement fits the " \
            "working set, which can inflate tau); free: the R-free set, as " \
            "for sigmaA."
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

def _framed_llgi_data(llgi_data, f_obs, d, u, t, a, **changes):
  """ llgi_data with the data d (renormalise_e) in the t-frame and the
  tNCS arrays; the cached hybrid is dropped """
  return replace_llgi_data(llgi_data,
    teps=f_obs.array(data=flex.double(f_obs.size(), 1.)),
    dobs=f_obs.array(data=d.dobs), feff=f_obs.array(data=d.feff),
    e_obs_sq=f_obs.array(data=d.e_obs_sq),
    sig_e_obs_sq=f_obs.array(data=d.sig_e_obs_sq),
    force_exact=f_obs.array(data=d.force_exact),
    tncs_a=f_obs.array(data=a), tncs_u=f_obs.array(data=u),
    tncs_t=f_obs.array(data=t), _hybrid_object=None, **changes)

def replace_llgi_data(llgi_data, **changes):
  import mmtbx.refinement.llgi_hybrid as llgi_hybrid
  return llgi_hybrid.replace_llgi_data(llgi_data, **changes)

def setup(llgi_data, f_obs, pdb_hierarchy, xray_structure, params=None,
      log=None):
  """ The tNCS treatment for llgi_data (on f_obs's index set), if the model
  has groups of chains related by tNCS: the observations are kept on the E
  scale without tNCS (nacelle's TEPS is not used), and renormalised by
  t0 = 1 + rho A, rho fitted to the data alone, until sigmaA and tau are
  first fitted (fit_and_renormalise). Returns the new llgi_data (Feff in
  it is the new f_obs), or None if there are no tNCS groups or the
  treatment is switched off. """
  if(params is None):
    params = llgi_tncs_params.extract()
  if(params.enabled is False):
    return None
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
    mod.f_total, linear=params.tau_linear_in_resolution)
  t0 = 1 + rho.at(dss) * mod.a
  d = renormalise_e(e_base, sig_base, llgi_data.resn.data(), t0, centric)
  if(log is not None):
    print("tNCS: %.3f of the model's scattering in tNCS groups; data-only "
      "modulation rho = %.3f%s" % (mod.f_total, rho.p0,
      "" if rho.p1 == 0 else " + %.3f d*^2/d*^2_max" % rho.p1), file=log)
  return _framed_llgi_data(llgi_data, f_obs, d, u=t0, t=t0, a=mod.a,
    e_obs_sq_base=f_obs.array(data=e_base),
    sig_e_obs_sq_base=f_obs.array(data=sig_base),
    tncs_state=group_args(groups=groups, tau=rho, rho=rho,
      f_total=mod.f_total, params=params))

def fit_and_renormalise(fmodel, scatfrac_params=None, e_params=None):
  """ Per macrocycle, for fmodel's current model and scales: A(h); sigmaA
  (E scale, R-free set) and tau (params.tau_reflections: all reflections,
  or the R-free set), alternately, in the u-frame; ScatFrac (F scale,
  working set) in the u-frame; then the data
  renormalised by t = sigmaA^2 (1 + A) + (1 - a_eff^2) u, with Feff the new
  f_obs. Returns a group_args like llgi_e_sigmaa.
  estimate_sigmaa_e_then_scatfrac_f's, plus .tau. """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
  llgi_data = fmodel.llgi_data()
  state = llgi_data.tncs_state
  params = state.params
  if(scatfrac_params is None):
    scatfrac_params = llgi_scatfrac.llgi_scatfrac_params.extract()
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  dss = f_obs.d_star_sq().data()
  eps = f_obs.epsilons().data().as_double()
  n = f_obs.size()
  resn = llgi_data.resn.data()
  e_base = llgi_data.e_obs_sq_base.data()
  sig_base = llgi_data.sig_e_obs_sq_base.data()
  rice_kappa = llgi_data.hybrid_params.rice_kappa \
    if getattr(llgi_data, "hybrid_params", None) is not None else 0.1
  a = pair_modulation(f_obs, state.groups, fmodel.xray_structure).a
  e_model = flex.abs(llgi_e_sigmaa.build_e_model(
    fmodel.f_model_no_aniso_scale().data(), eps, dss).e_model)
  tau = state.tau
  warnings = []
  for cycle in range(max(1, params.n_cycles)):
    u = 1 + tau.at(dss) * a
    d = renormalise_e(e_base, sig_base, resn, u, centric)
    hybrid_u = _hybrid(d, centric, rice_kappa)
    sigmaa_result = llgi_e_sigmaa.estimate_e_sigmaa_for_fmodel(fmodel,
      dobs=d.dobs, feff=d.eeff * resn, resn=resn, params=e_params,
      hybrid=hybrid_u, tncs_u=u)
    sigmaa = sigmaa_result.sigmaa
    tau_selection = flex.bool(n, True)
    if(getattr(params, "tau_reflections", "all") == "free"):
      tau_selection = fmodel.r_free_flags().data()
    tau = fit_tau(e_base, sig_base, flex.double(n, 1.), centric, e_model,
      sigmaa, a, dss, state.f_total, tau_selection,
      tau_start=[tau.p0, tau.p1], linear=params.tau_linear_in_resolution)
  if(sigmaa_result.lbfgs_error is not None):
    warnings.append("sigmaA fit: L-BFGS stopped early: %s" %
      sigmaa_result.lbfgs_error)
  u = 1 + tau.at(dss) * a
  d = renormalise_e(e_base, sig_base, resn, u, centric)
  scale_factor = fmodel.scale_ml_wrapper()
  scatfrac_result = llgi_scatfrac.estimate_llgi_scatfrac(
    f_eff=d.feff, working_selection=~fmodel.r_free_flags().data(),
    f_calc=fmodel.f_model().data(), dobs=d.dobs, sigmaa=sigmaa,
    teps=flex.double(n, 1.), resn=resn, centric_flags=centric,
    d_star_sq=dss, scale_factor=scale_factor, params=scatfrac_params,
    hybrid=_hybrid(d, centric, rice_kappa), tncs_u=u)
  if(scatfrac_result.lbfgs_error is not None):
    warnings.append("ScatFrac fit: L-BFGS stopped early: %s" %
      scatfrac_result.lbfgs_error)
  scatfrac = scatfrac_result.scatfrac
  form = getattr(scatfrac_params, "form", "likelihood")
  a_eff = llgi_scatfrac.effective_sigmaa(sigmaa, scatfrac, scale_factor,
    form)
  t = expected_intensity_modulation(sigmaa * sigmaa, a_eff, a, u)
  d = renormalise_e(e_base, sig_base, resn, t, centric)
  new = _framed_llgi_data(llgi_data, f_obs, d, u=u, t=t, a=a,
    sigmaa=f_obs.array(data=sigmaa),
    sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=f_obs.array(data=scatfrac), scatfrac_form=form,
    tncs_state=group_args(groups=state.groups, tau=tau, rho=state.rho,
      f_total=state.f_total, params=params),
    _sigmaa_f_model=fmodel.f_model().data().deep_copy())
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=f_obs.customized_copy(data=d.feff))
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac, scatfrac_form=form, target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac,
    n_scatfrac_at_floor=scatfrac_result.n_at_floor,
    warnings=warnings, tau=tau)

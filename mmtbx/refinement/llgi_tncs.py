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
of the true copies relative to the consensus (gamma <= 1, or the
cctbx normalisation's gamma_max; v the rms
difference per direction; one curve shared by all groups). The error
variance is then modulated by u(h) = (T - sigmaA^2 M)/(1 - sigmaA^2), and
the expected intensity by t(h) = sigmaA^2 M(h) + (1 - sigmaA^2) u(h) = T,
sigmaA being the coefficient of the unit-normalised model E-value. With
coherence = tau, u(h) = 1 + tau(s) A(h) on the model's own pair terms
instead. The LLGI functors take u and t per reflection (cctbx/xray/targets/llgi.h, tncs_frame); the data are normalised
by t (renormalise). Parameter fits use the u-frame (t = u), in which the
frame is valid for any sigmaA and which, with the exact likelihood, gives
the same result as any other frame.

With the cctbx normalisation (llgi_data.normalisation_source = "cctbx",
mmtbx.refinement.llgi_normalisation), RESN and TEPS_data = 1 +
Gamma_data(s) A_cons(h) are refitted to the intensities alone at each
stage (the data-only stage), and the frame t is TEPS_data or t_model
(params.frame), with u clamped where t <= (1 - sigmaA^2) u (setup,
fit_and_renormalise).
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

def _polymer_atoms(chain, selection=None, out=None):
  """ {(resseq, icode, atom name): (i_seq, resname)} for the first
  conformer of the amino-acid and nucleotide residues of a chain, without
  hydrogens; only the atoms with selection[i_seq] if selection (a flex.bool)
  is given; added to out if given """
  if(out is None): out = {}
  for rg in chain.residue_groups():
    ag = rg.atom_groups()[0]
    cls = iotbx.pdb.common_residue_names_get_class(ag.resname)
    if cls not in ("common_amino_acid", "common_rna_dna"):
      continue
    for a in ag.atoms():
      if a.element.strip().upper() in ("H", "D"):
        continue
      if(selection is not None and not selection[a.i_seq]):
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

  @classmethod
  def from_selections(cls, pdb_hierarchy, crystal_symmetry, groups,
                      angle_max=10.):
    """ Groups from atom selections instead of whole chains (domain-level
    tNCS). groups: a list of groups, each a list of atom-selection strings,
    one per copy (the first is the reference copy). The copies' polymer
    atoms are matched by residue number, insertion code and atom name; each
    copy's operator is the symmetry operator that brings it closest to a
    pure translation of the reference copy, and a copy related to the
    reference by more than angle_max degrees is an error. The result has the
    attributes of tncs_groups, with the selection strings as chain_ids and
    i_seqs indexing the atoms of pdb_hierarchy (counted from 0, whatever
    their i_seq values: those are reset on a copy of pdb_hierarchy, which is
    left unchanged, if they are not 0..n-1). """
    from libtbx.utils import Sorry
    n_atoms = pdb_hierarchy.atoms().size()
    if(list(pdb_hierarchy.atoms().extract_i_seq()) != list(range(n_atoms))):
      pdb_hierarchy = pdb_hierarchy.deep_copy()
      pdb_hierarchy.atoms().reset_i_seq()
    cache = pdb_hierarchy.atom_selection_cache()
    sites_frac = crystal_symmetry.unit_cell().fractionalize(
      pdb_hierarchy.atoms().extract_xyz())
    self = cls.__new__(cls)
    self.crystal_symmetry = crystal_symmetry
    self.n_atoms = n_atoms
    self.groups = []
    for i_group, selections in enumerate(groups):
      selections = list(selections)
      if(len(selections) < 2):
        raise Sorry("tNCS group %d: at least two copies are needed, not %d"
          % (i_group + 1, len(selections)))
      copies = []
      for sel_str in selections:
        sel = cache.selection(sel_str)
        atoms = {}
        for ch in pdb_hierarchy.models()[0].chains():
          _polymer_atoms(ch, selection=sel, out=atoms)
        if(len(atoms) == 0):
          raise Sorry("tNCS group %d: copy '%s' selects no polymer atoms"
            % (i_group + 1, sel_str))
        copies.append(atoms)
      common = set(copies[0])
      for atoms in copies[1:]:
        common &= set(atoms)
      common = sorted(k for k in common
        if all(atoms[k][1] == copies[0][k][1] for atoms in copies))
      if(len(common) < 3):
        raise Sorry("tNCS group %d: the copies have only %d atoms in common"
          % (i_group + 1, len(common)))
      i_seqs = [flex.size_t([atoms[k][0] for k in common]) for atoms in copies]
      ops, angles = [sgtbx.rt_mx()], [0.]
      ref_sites = sites_frac.select(i_seqs[0])
      for sel_str, sel in zip(selections[1:], i_seqs[1:]):
        angle, op = _best_op(ref_sites, sites_frac.select(sel),
                             crystal_symmetry)
        if(angle > angle_max):
          raise Sorry("tNCS group %d: copy '%s' is related to '%s' by %.1f "
            "degrees, more than angle_max" % (i_group + 1, sel_str,
            selections[0], angle))
        ops.append(op); angles.append(angle)
      self.groups.append(group_args(
        chain_ids=selections, i_seqs=i_seqs,
        chain_i_seqs=[flex.size_t(sorted(v[0] for v in atoms.values()))
                      for atoms in copies],
        ops=ops, angles=angles))
    return self

  def __len__(self):
    return len(self.groups)

  def show(self, log):
    if not self.groups:
      print("tNCS: no groups of chains related by tNCS", file=log)
      return
    def q(cid):   # selections with blanks in quotes
      return "'%s'" % cid if " " in cid else cid
    for g in self.groups:
      print("tNCS group: chains %s, %d common atoms" % (
        " ".join(q(c) for c in g.chain_ids), g.i_seqs[0].size()), file=log)
      for cid, op, angle in zip(g.chain_ids[1:], g.ops[1:], g.angles[1:]):
        print("  %s -> %s by %s: rotation %.2f degrees" % (
          q(cid), q(g.chain_ids[0]), op.as_xyz(), angle), file=log)

def groups_from_params(pdb_hierarchy, crystal_symmetry, params):
  """ The model's tNCS groups: from the atom selections of params.group
  (llgi_tncs_params; tncs_groups.from_selections) if any is given,
  otherwise found automatically from the chains (tncs_groups, with
  params.angle_max). """
  selected = [[c for c in g.copy if c is not None]
              for g in (getattr(params, "group", None) or [])]
  selected = [g for g in selected if g]
  if(selected):
    return tncs_groups.from_selections(pdb_hierarchy, crystal_symmetry,
      selected, angle_max=params.angle_max)
  return tncs_groups(pdb_hierarchy, crystal_symmetry,
    angle_max=params.angle_max)

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
  d = (sites_frac_m - sites_frac_n).as_double()
  d = flex.vec3_double(d - flex.floor(d))
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
      bases=("model", "consensus"), per_group=False):
  """ pair_modulation for several bases at once (the weights are computed
  once): group_args(a ({basis: flex.double}), fractions, f_total,
  a_groups). With per_group, a_groups is {basis: list of flex.double}, the
  pair modulation of each group (in groups.groups' order; a is their sum),
  otherwise None. """
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
  a_groups = dict((basis, []) for basis in bases) if per_group else None
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
    if(per_group):
      for basis in bases:
        a_groups[basis].append(flex.double(dss.size(), 0))
    for m, n in itertools.combinations(range(len(g.i_seqs)), 2):
      b_sum = b_iso.select(g.i_seqs[m]) + b_iso.select(g.i_seqs[n])
      weight = 2 * flex.sqrt(f_s[m] * f_s[n]) * _weight_sums(dss, z_sq, b_sum) \
        / flex.sqrt(w_common[m] * w_common[n])   # including the coherence
      for basis in bases:
        term = weight * pair_term(miller_set, placed[basis][m],
          placed[basis][n], z_sq, b_sum)
        a[basis] += term
        if(per_group): a_groups[basis][-1] += term
  return group_args(a=a, fractions=fractions,
    f_total=max(flex.max(f_in_groups), sum(sum(f) for f in fractions)),
    a_groups=a_groups)

def select_pair_terms(pm, source, target):
  """ pm, a pair_modulations result on source's indices, on target's
  indices (both miller sets): a, and a_groups if present, selected through
  miller.match_indices; fractions and f_total are kept. None if an index of
  target is not in source. """
  matches = miller.match_indices(source.indices(), target.indices())
  if(matches.singles(1).size() != 0):
    return None
  pairs = matches.pairs()
  i_source = pairs.column(0)
  # into the target's order (the pairs need not be in it)
  assert pairs.size() == target.size()
  order = flex.sort_permutation(pairs.column(1))
  i_source = i_source.select(order)
  a = dict((basis, v.select(i_source)) for basis, v in pm.a.items())
  a_groups = None
  if(pm.a_groups is not None):
    a_groups = dict((basis, [v.select(i_source) for v in vs])
      for basis, vs in pm.a_groups.items())
  return group_args(a=a, fractions=pm.fractions, f_total=pm.f_total,
    a_groups=a_groups)

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
      a_model, a_cons, d_star_sq, selection, start=None, gamma_max=1.0):
  """ Gamma(s) (coherence_curve) by the exact E-scale likelihood over the
  selected reflections, with sigmaA fixed, as fit_tau but with u from
  error_variance_modulation; 0 <= gamma <= gamma_max (1 for identical
  copies; the cctbx normalisation allows a little more, its gamma_max),
  0 <= v <= 3 A. The default start is gamma = 0.95, or start's gamma
  at most gamma_max - 0.01. """
  from scipy.optimize import minimize
  isel = selection.iselection()
  io, si = e_obs_sq.select(isel), sig_e_obs_sq.select(isel)
  cen, em = centric_flags.select(isel), e_model.select(isel)
  sa, am = sigmaa.select(isel), a_model.select(isel)
  ac, ss = a_cons.select(isel), d_star_sq.select(isel)
  n = isel.size()
  one = flex.double(n, 1.)
  def f(x):
    if not (0 <= x[0] <= gamma_max and 0 <= x[1] <= 3): return 1e10
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
  x0 = [0.95, 0.2] if start is None else [min(start.gamma,
    gamma_max - 0.01), max(start.v, 0.05)]
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

U_CLAMP_WARNING = ("tNCS: u clamped for %d reflections where t <= "
  "(1 - sigmaA^2) u")

def clamp_u_to_frame(sigmaa, u, t):
  """ u lowered to (1 - U_MIN) t/(1 - sigmaA^2) where t <= (1 - sigmaA^2) u,
  so that the frame of the LLGI functors (tncs_frame: a'^2 = 1 - (1 -
  sigmaA^2) u/t) is valid, a'^2 >= U_MIN there; only where 0 < sigmaA < 1
  (elsewhere the frame is invalid whatever u is). Returns (u, a new
  array; the number of reflections clamped). """
  sa2 = flex.pow2(sigmaa)
  clamp = (sigmaa > 0) & (sigmaa < 1) & (t <= (1 - sa2) * u)
  u = u.deep_copy()
  isel = clamp.iselection()
  u.set_selected(isel, (1 - U_MIN) * t.select(isel) / (1 - sa2.select(isel)))
  return u, isel.size()

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
  group
    .multiple = True
    .short_caption = tNCS group of domain copies
    .help = "A group of copies related by tNCS, given by atom selections " \\
            "(domain-level tNCS, e.g. for domains that have moved " \\
            "relative to each other between the copies). Each copy is a " \\
            "selection of the same residues (matched by residue number, " \\
            "insertion code and atom name); the first is the reference " \\
            "copy. If any group is given, the groups of chains are not " \\
            "searched for."
  {
    copy = None
      .type = atom_selection
      .multiple = True
      .short_caption = Atom selection of one copy
      .help = "Atoms of one copy of the group; give at least two."
  }
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
  frame = *model data
    .type = choice
    .short_caption = Expected-intensity modulation t of the tNCS frame
    .help = "The modulation t by which the data are normalised (the frame " \\
            "of the LLGI target), with the cctbx normalisation only " \\
            "(llgi_data.normalisation = cctbx). model: t = sigmaA^2 M + " \\
            "(1 - sigmaA^2) u from the model's pair terms and the fitted " \\
            "coherence; data: t = TEPS_data from the data-only " \\
            "normalisation fit, with u lowered where t <= (1 - sigmaA^2) u " \\
            "(reported). The exact likelihood does not depend on t; t " \\
            "affects the Rice/hybrid approximation and Feff."
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

def is_cctbx_source(llgi_data):
  """ llgi_data is from the cctbx normalisation (phenix.refinement.
  llgi_data.get_llgi_data_cctbx: .normalisation_source = "cctbx", with the
  intensities .i_obs_base) """
  return (getattr(llgi_data, "normalisation_source", None) == "cctbx" and
    getattr(llgi_data, "i_obs_base", None) is not None)

def _normalisation_fit(i_obs_base, a_groups, normalisation_params,
      selection=None, log=None):
  """ the data-only stage: mmtbx.refinement.llgi_normalisation.fit of the
  intensities with the consensus pair terms of the tNCS groups """
  import mmtbx.refinement.llgi_normalisation as llgi_normalisation
  return llgi_normalisation.fit(i_obs=i_obs_base, tncs_a_groups=a_groups,
    params=normalisation_params, log=log, selection=selection)

def _resn_scale(resn, sig_i_obs, sig_e_obs_sq):
  """ k in resn = k RESN_fit: llgi_data's RESN (which mmtbx.f_model.manager.
  apply_scale_k1_to_f_obs puts on the model's scale) over the RESN of the
  normalisation that gave sig_e_obs_sq = SIGI/RESN_fit^2, the median of the
  ratios """
  import numpy
  ratio = resn / flex.sqrt(sig_i_obs / sig_e_obs_sq)
  return float(numpy.median(ratio.as_numpy_array()))

def _base_from_fit(i_obs_base, fit, resn_scale):
  """ (RESN on llgi_data's scale, E_obs^2 and sigma(E_obs^2) without tNCS,
  I/RESN^2 with the fit's RESN) """
  resn = fit.resn.data()
  resn_sq = flex.pow2(resn)
  return (resn * resn_scale, i_obs_base.data() / resn_sq,
    i_obs_base.sigmas() / resn_sq)

def setup(llgi_data, f_obs, pdb_hierarchy, xray_structure, params=None,
      log=None, normalisation_params=None, pair_terms=None):
  """ The tNCS treatment for llgi_data (on f_obs's index set), if the model
  has groups of chains related by tNCS: the observations are kept on the E
  scale without tNCS (nacelle's TEPS is not used), and renormalised by
  t0 = 1 + rho A, rho fitted to the data alone, until sigmaA and the
  coherence (or tau) are first fitted (fit_and_renormalise). With the
  cctbx normalisation (is_cctbx_source), the data-only stage instead:
  llgi_normalisation.fit of llgi_data.i_obs_base with the consensus pair
  terms of the groups (normalisation_params: extracted
  llgi_normalisation.normalisation_params, or None for defaults; kept for
  fit_and_renormalise), RESN from it, the observations on the E scale
  I/RESN^2, and the frame started at u = t = TEPS_data; the fit's
  selection is kept for later stages. pair_terms (cctbx source only): the
  pair_modulations result (with "consensus", per_group=True) on f_obs's
  indices, if already known (select_pair_terms); otherwise computed here.
  Returns the new llgi_data (Feff in it is the new f_obs), or None if there
  are no tNCS groups or the treatment is switched off. """
  if(params is None):
    params = llgi_tncs_params.extract()
  if(params.enabled is False):
    return None
  assert getattr(llgi_data, "tncs_state", None) is None  # t would apply twice
  if(getattr(llgi_data, "e_obs_sq", None) is None):
    from libtbx.utils import Sorry
    raise Sorry("The LLGI treatment of translational NCS needs the "
      "intensities (I, SIGI); the LLGI data have none.")
  groups = groups_from_params(pdb_hierarchy, f_obs.crystal_symmetry(), params)
  if(log is not None): groups.show(log)
  if(len(groups) == 0):
    return None
  if(is_cctbx_source(llgi_data)):
    return _setup_cctbx(llgi_data, f_obs, groups, xray_structure, params,
      normalisation_params, log, pair_terms)
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
      f_total=mod.f_total, params=params, teps_data=None, t_model=None,
      normalisation=None, normalisation_params=None, monitor=None))

def _setup_cctbx(llgi_data, f_obs, groups, xray_structure, params,
      normalisation_params, log, pair_terms=None):
  """ setup with the cctbx normalisation (only the consensus pair terms
  are needed, and not computed if pair_terms is given) """
  if(normalisation_params is None):
    import mmtbx.refinement.llgi_normalisation as llgi_normalisation
    normalisation_params = llgi_normalisation.normalisation_params.extract()
  i_obs = llgi_data.i_obs_base
  assert i_obs.indices().all_eq(f_obs.indices())
  centric = f_obs.centric_flags().data()
  dss = f_obs.d_star_sq().data()
  pm = pair_terms
  if(pm is None):
    pm = pair_modulations(f_obs, groups, xray_structure, ("consensus",),
      per_group=True)
  for a in pm.a_groups["consensus"]:
    assert a.size() == f_obs.size()
  if(log is not None):
    print("tNCS: cctbx normalisation with the consensus pair terms of the "
      "tNCS groups (data-only stage)", file=log)
  fit = _normalisation_fit(i_obs, pm.a_groups["consensus"],
    normalisation_params, log=log)
  scale = _resn_scale(llgi_data.resn.data(), i_obs.sigmas(),
    llgi_data.sig_e_obs_sq.data() * llgi_data.teps.data())
  resn, e_base, sig_base = _base_from_fit(i_obs, fit, scale)
  teps_data = fit.teps.data()
  d = renormalise_e(e_base, sig_base, resn, teps_data, centric)
  gamma = fit.coherence.gamma
  if(log is not None):
    print("tNCS: %.3f of the model's scattering in tNCS groups; data-only "
      "%s; frame %s" % (pm.f_total, fit.coherence.summary(), params.frame),
      file=log)
  return _framed_llgi_data(llgi_data, f_obs, d, u=teps_data.deep_copy(),
    t=teps_data.deep_copy(), resn=f_obs.array(data=resn),
    e_obs_sq_base=f_obs.array(data=e_base),
    sig_e_obs_sq_base=f_obs.array(data=sig_base),
    normalisation_selection=f_obs.array(data=fit.selection),
    tncs_state=group_args(groups=groups,
      tau=linear_in_d_star_sq(flex.max(dss), min(max(gamma, 0.), 0.99)),
      coherence=fit.coherence, f_total=pm.f_total, params=params,
      teps_data=teps_data, t_model=None, normalisation=fit,
      normalisation_params=normalisation_params, monitor=None))

def _rice_kappa(llgi_data):
  return llgi_data.hybrid_params.rice_kappa \
    if getattr(llgi_data, "hybrid_params", None) is not None else 0.1

def _model_anisotropy(e_params, llgi_data):
  # llgi_e_sigmaa's model_anisotropy (resolved), as in the sigmaA fit
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  import mmtbx.refinement.llgi_model_anisotropy as llgi_model_anisotropy
  if(e_params is None):
    e_params = llgi_e_sigmaa.llgi_e_sigmaa_params.extract()
  return llgi_model_anisotropy.resolve(
    getattr(e_params, "model_anisotropy", "off"), llgi_data)

def _fit_scatfrac(fmodel, d, u, sigmaa, scatfrac_params, rice_kappa,
      warnings, model_anisotropy="off", resn=None):
  """ ScatFrac (F scale, working set) in the u-frame, for the data d
  (renormalise_e by u); model_anisotropy as llgi_e_sigmaa's; resn: RESN
  of d (None: llgi_data's) """
  import mmtbx.refinement.llgi_scatfrac as llgi_scatfrac
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  if(scatfrac_params is None):
    scatfrac_params = llgi_scatfrac.llgi_scatfrac_params.extract()
  for w in llgi_e_sigmaa.isotropic_scatfrac_warnings(scatfrac_params,
      model_anisotropy):
    if(w not in warnings): warnings.append(w)
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  if(resn is None):
    resn = fmodel.llgi_data().resn.data()
  result = llgi_scatfrac.estimate_llgi_scatfrac(
    f_eff=d.feff, working_selection=~fmodel.r_free_flags().data(),
    f_calc=fmodel.f_model().data(), dobs=d.dobs, sigmaa=sigmaa,
    teps=flex.double(f_obs.size(), 1.), resn=resn,
    centric_flags=centric, d_star_sq=f_obs.d_star_sq().data(),
    scale_factor=fmodel.scale_ml_wrapper(), params=scatfrac_params,
    hybrid=_hybrid(d, centric, rice_kappa), tncs_u=u, miller_set=f_obs,
    model_anisotropy=model_anisotropy)
  if(result.lbfgs_error is not None):
    warnings.append("ScatFrac fit: L-BFGS stopped early: %s" %
      result.lbfgs_error)
  return result

def refit_sigmaa_scatfrac(fmodel, scatfrac_params=None, e_params=None):
  """ sigmaA and ScatFrac for fmodel's current model in the current tNCS
  frame (u and t, and so the data, are kept; A(h) is not recomputed): for
  refits on a copy of fmodel whose model is changed temporarily, e.g. with
  atoms removed or B-factors shaken (mmtbx.refinement.weights.
  adp_gradients). sigmaA is fitted in the stored u-frame, and u is then
  clamped for the new sigmaA and the stored t (clamp_u_to_frame; the number
  is the result's .n_u_clamped, with a warning), whatever the normalisation
  source: a lower sigmaA can invalidate the frame of t_model as well.
  Returns a group_args as fit_and_renormalise. """
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
    rice_kappa, warnings, model_anisotropy=_model_anisotropy(e_params,
      fmodel.llgi_data()))
  # as in fit_and_renormalise: ScatFrac with the u of the fits, then u
  # clamped for the frame (t is kept, so the data are; the cached hybrid
  # does not depend on u)
  u, n_clamped = clamp_u_to_frame(sigmaa, u, llgi_data.tncs_t.data())
  if(n_clamped > 0):
    warnings.append(U_CLAMP_WARNING % n_clamped)
  fmodel.set_llgi_data(replace_llgi_data(llgi_data,
    tncs_u=f_obs.array(data=u),
    sigmaa=f_obs.array(data=sigmaa), sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=f_obs.array(data=scatfrac_result.scatfrac),
    _sigmaa_f_model=fmodel.f_model().data().deep_copy()))
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac_result.scatfrac, target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac, warnings=warnings, tau=None,
    coherence=None, n_u_clamped=n_clamped, tncs_summary=None)

# ---------------------------------------------------------------------------
# Monitor of the data-based against the model-based TEPS

MONITOR_GAP_GROWTH = 0.002   # nats/reflection
MONITOR_STANDS_OUT = 3.      # times the median absolute coefficient
MONITOR_MIN_COEFFICIENT = 0.05
MONITOR_N_BINS = 5
MONITOR_MIN_CELL = 10        # reflections for a meaningful correlation
MONITOR_REPORT_CHANGE = 0.1  # coefficient change that reports a stand-out
                             # group again

GAP_WARNING = ("tNCS monitor: free-set gap grew from %.4f to %.4f "
  "nats/reflection")
STANDS_OUT_WARNING = "tNCS monitor: group %d stands out (coefficient %.3f)"

def _teps_class(teps):
  """ "weak" (TEPS < 0.5), "neutral" (0.5-1.5) or "strong" (> 1.5) """
  return ["weak" if t < 0.5 else ("strong" if t > 1.5 else "neutral")
    for t in teps]

def _correlation(x, y):
  """ Pearson correlation; 1 if x and y are both constant and equal, 0 if
  either is constant otherwise """
  import numpy
  sx, sy = numpy.std(x), numpy.std(y)
  if(sx < 1e-12 or sy < 1e-12):
    return 1. if (sx < 1e-12 and sy < 1e-12 and
      abs(numpy.mean(x) - numpy.mean(y)) < 1e-9) else 0.
  return float(numpy.corrcoef(x, y)[0, 1])

def _monitor_table(teps_data, t_model, d_spacings):
  """ rows (d_max, d_min, class, n, corr_ln, mean_ratio, rms_ln_ratio) for
  MONITOR_N_BINS resolution bins of equal numbers of reflections, each split
  by the class of TEPS_data; corr_ln: correlation of ln TEPS_data with
  ln t_model, mean_ratio: mean of TEPS_data/t_model, rms_ln_ratio: rms of
  ln(TEPS_data/t_model) """
  import numpy
  teps = teps_data.as_numpy_array()
  t = t_model.as_numpy_array()
  d = d_spacings.as_numpy_array()
  classes = numpy.array(_teps_class(teps))
  order = numpy.argsort(-d)    # low resolution first
  rows = []
  for chunk in numpy.array_split(order, MONITOR_N_BINS):
    if(chunk.size == 0): continue
    d_max, d_min = float(d[chunk].max()), float(d[chunk].min())
    for c in ("weak", "neutral", "strong"):
      sel = chunk[classes[chunk] == c]
      if(sel.size == 0): continue
      ln_ratio = numpy.log(teps[sel] / t[sel])
      rows.append((d_max, d_min, c, int(sel.size),
        _correlation(numpy.log(teps[sel]), numpy.log(t[sel])),
        float(numpy.mean(teps[sel] / t[sel])),
        float(numpy.sqrt(numpy.mean(ln_ratio ** 2)))))
  return rows

def _group_coefficients(teps_data, t_model, a_groups):
  """ least-squares coefficients of ln(TEPS_data/t_model) on each A_g """
  import numpy
  y = numpy.log(teps_data.as_numpy_array() / t_model.as_numpy_array())
  x = numpy.array([a.as_numpy_array() for a in a_groups]).T
  return [float(c) for c in numpy.linalg.lstsq(x, y, rcond=None)[0]]

def _free_log_likelihood(i_obs, resn, t, free_flags):
  """ mean over the free reflections with SIGI > 0 (all such, if there are
  no free ones) of the Wilson log-likelihood of the intensities with
  Sigma = RESN^2 t (nats per reflection) """
  sig = i_obs.sigmas()
  sel = (sig > 0) & free_flags
  if(sel.count(True) == 0): sel = sig > 0
  isel = sel.iselection()
  r = xray_ext.wilson_intensity(i_obs=i_obs.data().select(isel),
    sig_i_obs=sig.select(isel), sigma=(flex.pow2(resn) * t).select(isel),
    centric_flags=i_obs.centric_flags().data().select(isel))
  return flex.mean(r.log_likelihood)

def monitor(teps_data, t_model, f_obs, free_flags, a_groups, i_obs, resn,
      coherence_data, coherence_fc, previous=None):
  """ Compares the data-only TEPS (teps_data) with the model-based
  modulation t_model: group_args(table, ll_free_data, ll_free_model,
  group_coefficients, coherence_data, coherence_fc, reported, warnings).
  table: rows (d_max, d_min, class, n, corr_ln, mean_ratio, rms_ln_ratio)
  (_monitor_table); ll_free_data and ll_free_model: nats per free
  reflection of the intensities i_obs (miller array with sigmas, on the
  scale of resn: Sigma = resn^2 t) with t = teps_data and t_model;
  group_coefficients: least-squares coefficients of ln(TEPS_data/t_model)
  on the A_g in a_groups (list of flex.double, one per group, f_obs's
  order); coherence_data and coherence_fc: the data-only and the phased
  coherence_curve. previous: the monitor result of the last macrocycle.
  Warnings: the gap ll_free_data - ll_free_model grew by more than
  MONITOR_GAP_GROWTH since previous; a group's coefficient exceeds
  MONITOR_STANDS_OUT times the median absolute coefficient (the lower
  median, so that a single stand-out among two groups can be seen) and
  MONITOR_MIN_COEFFICIENT. A stand-out group is warned about once:
  reported is a dict {group index (0-based): coefficient at its last
  report} of the groups that stand out now. A group is warned about if it
  is not in previous.reported (a previous without reported counts as
  empty), or if its coefficient differs from the reported one by more than
  MONITOR_REPORT_CHANGE; otherwise its reported coefficient is kept. A group
  that no longer stands out is dropped from reported, so that it is
  reported again if it stands out later. """
  warnings = []
  ll_data = _free_log_likelihood(i_obs, resn, teps_data, free_flags)
  ll_model = _free_log_likelihood(i_obs, resn, t_model, free_flags)
  if(previous is not None):
    gap_before = previous.ll_free_data - previous.ll_free_model
    gap = ll_data - ll_model
    if(gap - gap_before > MONITOR_GAP_GROWTH):
      warnings.append(GAP_WARNING % (gap_before, gap))
  coefficients = _group_coefficients(teps_data, t_model, a_groups)
  magnitudes = sorted(abs(c) for c in coefficients)
  median = magnitudes[(len(magnitudes) - 1) // 2]
  before = getattr(previous, "reported", None) or {}
  reported = {}
  for i, c in enumerate(coefficients):
    if(abs(c) > MONITOR_STANDS_OUT * median and abs(c) > MONITOR_MIN_COEFFICIENT):
      if(i not in before or abs(c - before[i]) > MONITOR_REPORT_CHANGE):
        warnings.append(STANDS_OUT_WARNING % (i + 1, c))
        reported[i] = c
      else:
        reported[i] = before[i]
  return group_args(
    table=_monitor_table(teps_data, t_model, f_obs.d_spacings().data()),
    ll_free_data=ll_data, ll_free_model=ll_model,
    group_coefficients=coefficients, coherence_data=coherence_data,
    coherence_fc=coherence_fc, reported=reported, warnings=warnings)

def show_monitor(result, log):
  """ prints the coherence curves (gamma, v) and the table of monitor's
  result """
  print("tNCS monitor: data-only %s; phased %s" % (
    result.coherence_data.summary(), result.coherence_fc.summary()), file=log)
  print("  free-set log-likelihood per reflection: data-based TEPS %.4f, "
    "model-based %.4f (gap %.4f)" % (result.ll_free_data, result.ll_free_model,
    result.ll_free_data - result.ll_free_model), file=log)
  print("  group coefficients: %s" % " ".join("%.3f" % c
    for c in result.group_coefficients), file=log)
  print("  %7s %7s %-8s %7s %8s %10s %12s" % ("d_max", "d_min", "class", "n",
    "corr_ln", "mean_ratio", "rms_ln_ratio"), file=log)
  for row in result.table:
    print("  %7.2f %7.2f %-8s %7d %8.3f %10.3f %12.3f" % row, file=log)

def _current_groups(state, crystal_symmetry, xray_structure, pdb_hierarchy):
  """ (groups, note): state.groups, or, if the model has gained or lost
  atoms since they were found (e.g. waters, ions), the groups found again
  in pdb_hierarchy, which must match xray_structure, and a note saying so
  (None otherwise) """
  from libtbx.utils import Sorry
  n_atoms = xray_structure.scatterers().size()
  if(state.groups.n_atoms == n_atoms):
    return state.groups, None
  if(pdb_hierarchy is None or pdb_hierarchy.atoms().size() != n_atoms):
    raise Sorry("tNCS: the model now has %d atoms, not the %d for which the "
      "tNCS groups were found, and no model with the current atoms was "
      "given to find them again." % (n_atoms, state.groups.n_atoms))
  groups = groups_from_params(pdb_hierarchy, crystal_symmetry, state.params)
  if(len(groups) == 0):
    raise Sorry("tNCS: no groups of chains related by tNCS in the model with "
      "%d atoms (they were found in the model with %d atoms)." % (n_atoms,
      state.groups.n_atoms))
  return groups, "groups found again for the model's %d atoms (previously " \
    "%d)" % (n_atoms, state.groups.n_atoms)

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
  atoms has changed.

  With the cctbx normalisation (setup), first the data-only stage:
  llgi_normalisation.fit with the current consensus pair terms, on the
  reflections selected at setup (n_would_change is reported), and RESN
  and the observations on the E scale (I/RESN^2) from it; then the fits as
  above; then t = TEPS_data (params.frame = data) or t_model = sigmaA^2
  (1 + A) + (1 - sigmaA^2) u (model), u clamped where t <= (1 - sigmaA^2) u
  (clamp_u_to_frame, reported), and the data renormalised by t.
  tncs_state then has .teps_data, .t_model and .normalisation (the fit),
  and the result .n_would_change and .n_u_clamped (None otherwise). """
  import mmtbx.refinement.llgi_e_sigmaa as llgi_e_sigmaa
  llgi_data = fmodel.llgi_data()
  state = llgi_data.tncs_state
  params = state.params
  warnings = []
  groups, regrouped = _current_groups(state, fmodel.f_obs().crystal_symmetry(),
    fmodel.xray_structure, pdb_hierarchy)
  f_obs = fmodel.f_obs()
  centric = f_obs.centric_flags().data()
  dss = f_obs.d_star_sq().data()
  n = f_obs.size()
  resn = llgi_data.resn.data()
  e_base = llgi_data.e_obs_sq_base.data()
  sig_base = llgi_data.sig_e_obs_sq_base.data()
  rice_kappa = _rice_kappa(llgi_data)
  physical = params.coherence == "physical"
  cctbx_source = getattr(state, "normalisation", None) is not None
  pm = pair_modulations(f_obs, groups, fmodel.xray_structure,
    ("model", "consensus") if (physical or cctbx_source) else ("model",),
    per_group=cctbx_source)
  a, a_cons = pm.a["model"], pm.a.get("consensus")
  normalisation = None
  if(cctbx_source):
    # the data-only stage at the current geometry, on the reflections
    # selected at setup
    i_obs = llgi_data.i_obs_base
    normalisation = _normalisation_fit(i_obs, pm.a_groups["consensus"],
      state.normalisation_params,
      selection=llgi_data.normalisation_selection.data())
    resn, e_base, sig_base = _base_from_fit(i_obs, normalisation,
      _resn_scale(resn, i_obs.sigmas(), sig_base))
    warnings.extend(normalisation.warnings)
  tau, coherence = state.tau, state.coherence
  # the bound of the coherence: the normalisation's gamma_max with the cctbx
  # source (state.normalisation_params, kept by _setup_cctbx), 1 with nacelle's
  gamma_max = state.normalisation_params.gamma_max if cctbx_source else 1.0
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
        a, a_cons, dss, fit_selection, start=coherence, gamma_max=gamma_max)
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
    rice_kappa, warnings, model_anisotropy=_model_anisotropy(e_params,
      fmodel.llgi_data()), resn=resn)
  scatfrac = scatfrac_result.scatfrac
  t_model = expected_intensity_modulation(sigmaa, a, u)
  t = t_model
  teps_data = n_clamped = monitor_result = None
  changes = {}
  if(cctbx_source):
    teps_data = normalisation.teps.data()
    if(params.frame == "data"):
      t = teps_data
    u, n_clamped = clamp_u_to_frame(sigmaa, u, t)
    if(n_clamped > 0):
      warnings.append(U_CLAMP_WARNING % n_clamped)
    monitor_result = monitor(teps_data, t_model, f_obs,
      fmodel.r_free_flags().data(), pm.a_groups["consensus"], i_obs,
      normalisation.resn.data(), normalisation.coherence, coherence,
      previous=getattr(state, "monitor", None))
    warnings.extend(monitor_result.warnings)
    changes = dict(resn=f_obs.array(data=resn),
      e_obs_sq_base=f_obs.array(data=e_base),
      sig_e_obs_sq_base=f_obs.array(data=sig_base))
  d = renormalise_e(e_base, sig_base, resn, t, centric)
  new = _framed_llgi_data(llgi_data, f_obs, d, u=u, t=t,
    sigmaa=f_obs.array(data=sigmaa),
    sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=f_obs.array(data=scatfrac),
    tncs_state=group_args(groups=groups, tau=tau,
      coherence=coherence, f_total=f_total, params=params,
      teps_data=teps_data, t_model=t_model, normalisation=normalisation,
      normalisation_params=getattr(state, "normalisation_params", None),
      monitor=monitor_result),
    _sigmaa_f_model=fmodel.f_model().data().deep_copy(), **changes)
  fmodel.set_llgi_data(new)
  fmodel.update(f_obs=f_obs.customized_copy(data=d.feff))
  return group_args(
    sigmaa=sigmaa, sigmaa_curve=sigmaa_result.evaluate_at,
    scatfrac=scatfrac, target=scatfrac_result.target,
    scatfrac_inf=scatfrac_result.scatfrac_inf,
    b_scatfrac=scatfrac_result.b_scatfrac,
    warnings=warnings, tau=None if physical else tau,
    coherence=coherence if physical else None,
    n_would_change=None if normalisation is None else
      normalisation.n_would_change,
    n_u_clamped=n_clamped,
    tncs_summary=(coherence.summary() if physical else "tau = %.3f%s" % (
      tau.p0, "" if tau.p1 == 0 else " + %.3f d*^2/d*^2_max" % tau.p1)) +
      ("" if normalisation is None else "; data-only %s, %d reflection(s) "
        "would change outlier status (selection kept); frame %s" % (
        normalisation.coherence.summary(), normalisation.n_would_change,
        params.frame)) +
      ("" if regrouped is None else "; " + regrouped))

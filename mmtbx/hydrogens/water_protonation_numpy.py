"""H-bond-aware placement of the two hydrogens on bare water oxygens.

For every bare water oxygen (any common water residue: HOH, DOD, H2O, WAT,
OH2, ...) the two H are placed pointing at H-bond acceptors, clear of the
whole structure (including H placed on other waters) and out of the
hemisphere of nearby metal cations. Geometry only: no map, no monomer
library. O-H is 0.984 A (neutron) or 0.957 A (X-ray) and H-O-H is
104.5 deg. Given a crystal symmetry, waters at a lattice contact also see
the neighbouring asymmetric units.

The public entry point :func:`place_water_hydrogens` modifies a hierarchy
in place; :class:`mmtbx.programs.water_protonation.Program` wraps it as the
``mmtbx.development.naiad`` command line.

Needs ``numpy`` and ``scipy`` (KDTree).
"""

from __future__ import absolute_import, division, print_function

import math
import random

import iotbx.pdb
from cctbx.crystal import super_cell
from iotbx.pdb.utils import check_for_missing_elements
from libtbx import group_args
import numpy as np
from scitbx import matrix
from scitbx.array_family import flex
from scipy.spatial import KDTree


# ---------------------------------------------------------------------------
# Constants for the H-bond-aware water-H placer (``place_water_hydrogens``).
# ---------------------------------------------------------------------------

_WATER_OH_XRAY = 0.957        # canonical X-ray O-H bond length (A)
_WATER_OH_NEUTRON = 0.984     # canonical neutron O-H bond length (A)
_WATER_HOH_DEG = 104.5        # canonical H-O-H angle (deg)
_WATER_ACCEPTOR_RADIUS = 3.5  # max distance to search for H-bond acceptors (A)
_WATER_ACCEPTOR_ELEMENTS = frozenset({"O", "N", "F", "S", "CL"})
_WATER_NH_BOND = 1.3          # max N-H distance for the "N carries an H" test (A)
# Lone-pair-directed placement (opt-in): H1 aims at an acceptor's lone-pair
# lobe rather than its nucleus. Its distances and angles:
_WATER_BOND_HEAVY = 1.9       # max heavy-heavy bond distance for lobe geometry (A)
_WATER_HBOND_HA = 1.8         # nominal H...acceptor distance for the lobe target (A)
_WATER_SP2_LOBE_DEG = 60.0    # half-angle of the two sp2 carbonyl/-late lobes
_WATER_CONE_SAMPLES = 36      # angular samples around the O-H1 cone
# Element-aware "clash-free" thresholds: a candidate H must clear every
# heavy atom by _WATER_MIN_CLEARANCE and every hydrogen (and cation) by the
# larger _WATER_MIN_H_CLEARANCE.
_WATER_MIN_CLEARANCE = 1.5
_WATER_MIN_H_CLEARANCE = 2.0
_WATER_CLEARANCE_RADIUS = 3.0  # neighbour search radius for clearance (A)
# Distance under which a symmetry equivalent is the same atom, not a copy.
_WATER_SYM_EQUIV_TOL = 0.5

# Relaxation sweeps after the greedy pass, each re-placing every water
# against the final positions of all the others.
_WATER_REFINE_SWEEPS = 5
# Early-stop tolerance: keep refining only while a sweep removes at least
# this many close (<2.0 A) H-H contacts.
_WATER_REFINE_TOL = 1
# Basin-hopping (opt-in): each round restarts from the best state, randomly
# re-orients the still-clashing waters, and relaxes. Seeded.
_WATER_BASIN_SEED = 0
_WATER_BASIN_RELAX = 2        # relaxation sweeps after each random kick
# Cations within _WATER_CATION_RADIUS of the water O are repellers, keeping
# both H in the hemisphere away from them, rather than acceptors.
_WATER_CATION_ELEMENTS = frozenset({
  "LI", "NA", "K", "RB", "CS",                                    # alkali
  "BE", "MG", "CA", "SR", "BA",                                   # alkaline earth
  "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",      # 3d
  "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",      # 4d
  "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",            # 5d
  "AL", "GA", "IN", "SN", "SB", "TL", "PB", "BI",                 # p-block
  "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",                 # lanthanides
  "TB", "DY", "HO", "ER", "TM", "YB", "LU",
  "TH", "U",                                                      # actinides
})
_WATER_CATION_RADIUS = 3.0
_WATER_METAL_COORD_RADIUS = 2.6  # first-shell M-O bond, reporting only
_WATER_H1_ACC_ALIGN = 0.7        # min cos(O-H1, acceptor) to count as donated to


def _has_deuterium(hier):
  """True if the model carries any D atom."""
  return (hier.atoms().extract_element(strip=True) == "D").count(True) > 0


def count_environment_hydrogens(hier):
  """Number of H/D carried by atoms outside water residues."""
  return sum(_n_hd(ag) for ag in hier.atom_groups()
             if not _is_water(ag.resname))


def _hd_flags(ag):
  """``(atoms array, per-atom H/D flags)`` for one atom_group."""
  # Callers index this array, never ``for a in ats``: af_shared_atom has no
  # __iter__, so a Python loop ends on a Boost.Python IndexError, ~66 us.
  ats = ag.atoms()
  e = ats.extract_element(strip=True)
  return ats, (e == "H") | (e == "D")


def _n_hd(ag):
  """Number of H/D in one atom_group."""
  e = ag.atoms().extract_element(strip=True)
  return (e == "H").count(True) + (e == "D").count(True)


def _is_water(resname):
  """True for any common water alias (HOH, DOD, H2O, WAT, OH2, ...)."""
  return (iotbx.pdb.common_residue_names_get_class(resname.strip().upper())
          == "common_water")


def _fibonacci_sphere(n):
  """``n`` roughly-uniform unit vectors over the sphere (Fibonacci spiral)."""
  golden = math.pi * (3.0 - math.sqrt(5.0))
  pts = []
  for i in range(n):
    y = 1.0 - 2.0 * (i + 0.5) / n
    r = math.sqrt(max(0.0, 1.0 - y * y))
    th = golden * i
    pts.append((math.cos(th) * r, y, math.sin(th) * r))
  return tuple(pts)


_WATER_FALLBACK_DIRECTIONS = _fibonacci_sphere(64)

# Shared empty neighbour-slot list.
_EMPTY_SLOTS = np.zeros(0, dtype=np.intp)


def _as_np(cols):
  """Stack scitbx col vectors into an (n, 3) float array."""
  if not cols:
    return np.zeros((0, 3))
  return np.array([c.elems for c in cols], dtype=float)


# The H2 cone is sampled at fixed azimuths, so its cosines and sines are
# constants, as is the normalized fallback sphere.
_CONE_COS = np.array([math.cos(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                      for k in range(_WATER_CONE_SAMPLES)])
_CONE_SIN = np.array([math.sin(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                      for k in range(_WATER_CONE_SAMPLES)])
_FALLBACK_NP = _as_np([matrix.col(v).normalize()
                       for v in _WATER_FALLBACK_DIRECTIONS])


def _ortho_frame(d1):
  """Two unit vectors ``p, q`` making ``(d1, p, q)`` an orthonormal frame.

  Mirrors ``scitbx.matrix.col.ortho()`` operation for operation.
  """
  x, y, z = d1[0], d1[1], d1[2]
  a, b, c = abs(x), abs(y), abs(z)
  if c <= a and c <= b:
    p = np.array((-y, x, 0.0))
  elif b <= a and b <= c:
    p = np.array((-z, 0.0, x))
  else:
    p = np.array((0.0, -z, y))
  p = p / math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
  return p, np.array((y * p[2] - p[1] * z,
                      z * p[0] - p[2] * x,
                      x * p[1] - p[0] * y))


def _rand_unit(rng):
  """A random unit vector, from Gaussian deviates of ``rng``."""
  while True:
    v = np.array((rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)))
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n > 1e-6:
      return v / n


def _strip_water_hydrogens(hier):
  """Remove every H/D from water residues.

  Returns the elements removed, as ``{atom_group memory_id: [element, ...]}``
  over the waters that carried any.
  """
  stripped = {}
  for ag in hier.atom_groups():
    if not _is_water(ag.resname):
      continue
    ats, hd = _hd_flags(ag)
    removed = [ats[int(k)] for k in hd.iselection()]
    if not removed:
      continue
    stripped[ag.memory_id()] = [a.element.strip().upper() for a in removed]
    for a in removed:
      ag.remove_atom(a)
  return stripped


def _new_h_atom(name, element, xyz, occ, b, hetero=False):
  """Build a hierarchy atom for a placed water H."""
  atom = iotbx.pdb.hierarchy.atom()
  atom.name = name
  atom.element = element
  atom.xyz = xyz
  atom.occ = occ
  atom.b = b
  atom.segid = " " * 4
  atom.hetero = hetero
  return atom


def _free_proton_name(existing_names, element):
  """First of ``<element>1`` / ``<element>2`` not in ``existing_names``.

  Returns a 4-character PDB atom name, or None if both are taken.
  """
  for di in (1, 2):
    name = f"{element}{di}"
    if name not in existing_names:
      return f" {name} "
  return None


def _sp2_plane_normal(atoms, static_tree, c, exclude):
  """Unit normal of the sp2 plane around atom ``c``.

  Built from ``c``'s heavy substituents other than ``exclude``, the bonded
  acceptor O that supplies one in-plane vector. None if underdetermined.
  """
  C = matrix.col(atoms[c].xyz)
  oc = matrix.col(atoms[exclude].xyz) - C   # C -> O
  for k in static_tree.query_ball_point(atoms[c].xyz, _WATER_BOND_HEAVY):
    if k == c or k == exclude:
      continue
    if atoms[k].element_is_hydrogen():
      continue
    v = matrix.col(atoms[k].xyz) - C
    if v.length() > _WATER_BOND_HEAVY:
      continue
    nrm = oc.cross(v)
    if nrm.length() > 1e-3:
      return nrm.normalize()
  return None


def _acceptor_lobes(atoms, static_tree, donor_n):
  """Lone-pair lobe unit vectors per acceptor atom, from bonded geometry.

  - terminal O (1 bond, carbonyl/carboxylate): two in-plane sp2 lobes
    ``2 * _WATER_SP2_LOBE_DEG`` apart, straddling the direction away from the
    bonded atom.
  - otherwise: one lobe opposite the sum of the bond directions.

  ``donor_n`` holds the indices of N that carry an H (donors, not acceptors).
  Returns acceptor index -> lobe list, empty where the geometry is
  underdetermined.
  """
  lobes = {}
  for i, a in enumerate(atoms):
    el = a.element.strip().upper()
    if el not in _WATER_ACCEPTOR_ELEMENTS:
      continue
    if el == "N" and i in donor_n:
      continue  # protonated N is a donor, not an acceptor
    A = matrix.col(a.xyz)
    nbrs = []
    for j in static_tree.query_ball_point(a.xyz, _WATER_BOND_HEAVY):
      if j == i:
        continue
      d = (matrix.col(atoms[j].xyz) - A).length()
      lim = (_WATER_NH_BOND if atoms[j].element_is_hydrogen()
             else _WATER_BOND_HEAVY)
      # An atom on top of this one contributes no bond direction.
      if 1e-3 < d <= lim:
        nbrs.append(j)
    bond_dirs = [(matrix.col(atoms[j].xyz) - A).normalize() for j in nbrs]
    if not bond_dirs:
      lobes[i] = []
    elif el == "O" and len(nbrs) == 1:
      away = bond_dirs[0] * -1.0
      n = _sp2_plane_normal(atoms, static_tree, nbrs[0], i)
      if n is None:
        lobes[i] = [away]
      else:
        lobes[i] = [
          away.rotate_around_origin(axis=n, angle=_WATER_SP2_LOBE_DEG, deg=True),
          away.rotate_around_origin(axis=n, angle=-_WATER_SP2_LOBE_DEG, deg=True)]
    else:
      bsum = matrix.col((0.0, 0.0, 0.0))
      for b in bond_dirs:
        bsum = bsum + b
      lobes[i] = [(bsum * -1.0).normalize()] if bsum.length() > 1e-6 else []
  return lobes


def _symmetry_environment(hier, sites_cart, crystal_symmetry, radius,
                          min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Copies of the atoms crystal symmetry places within ``radius`` of a water.

  Each copy carries its whole residue group, so an acceptor keeps the bonded
  neighbours its lobe geometry needs and an N keeps the H that marks it a
  donor. Copies are environment only: the waters to protonate come from the
  asymmetric unit, and a mate's placed protons are tracked separately (see
  :func:`_water_image_neighbours`). An equivalent closer than
  ``min_distance_sym_equiv`` to its own site counts as that same atom, which
  the asymmetric unit already holds.

  An atom's image lands within ``radius`` of a water exactly when the atom
  lies within ``radius`` of that water's inverse image, so the operators go on
  the few water oxygens instead of on every atom. A pair table cannot be
  seeded on the waters and costs several times as much.

  Returns ``(hierarchies, atoms, xyz)``: the sub-hierarchies, which own the
  atom objects and must be kept alive, their atoms in coordinate order, and
  their sites. All three are empty when symmetry places nothing in range.
  """
  o_sel = flex.size_t([a.i_seq for a in hier.atoms()
                       if _is_water(a.parent().resname)
                       and a.element.strip().upper() == "O"])
  if not o_sel.size():
    return [], [], flex.vec3_double()
  o_sites = sites_cart.select(o_sel)
  unit_cell = crystal_symmetry.unit_cell()
  tree = KDTree(sites_cart.as_numpy_array())
  # Deposited coordinates need not lie inside one cell, so the translations to
  # try span the fractional spread of the sites, not just the radius.
  margin = [radius * x for x in unit_cell.reciprocal_parameters()[:3]]
  fr_water = unit_cell.fractionalize(o_sites).parts()
  lo_water = [flex.min(c) for c in fr_water]
  hi_water = [flex.max(c) for c in fr_water]
  tol_sq = min_distance_sym_equiv ** 2
  rg_of = {}
  for rg in hier.residue_groups():
    idx = [a.i_seq for a in rg.atoms()]
    for i in idx:
      rg_of[i] = idx
  by_key = {}
  for op in crystal_symmetry.space_group():
    rot = unit_cell.matrix_cart(op.r())
    trn = unit_cell.orthogonalize(op.t().as_double())
    rot_inv = matrix.sqr(rot).transpose().elems
    fwd = unit_cell.fractionalize(rot * sites_cart + trn).parts()
    spans = []
    for axis in range(3):
      lo = lo_water[axis] - flex.max(fwd[axis]) - margin[axis]
      hi = hi_water[axis] - flex.min(fwd[axis]) + margin[axis]
      spans.append(range(int(math.floor(lo)), int(math.ceil(hi)) + 1))
    for i in spans[0]:
      for j in spans[1]:
        for k in spans[2]:
          shift = unit_cell.orthogonalize((i, j, k))
          off = (trn[0] + shift[0], trn[1] + shift[1], trn[2] + shift[2])
          near = tree.query_ball_point(
            (rot_inv * (o_sites - off)).as_numpy_array(), radius,
            return_sorted=False)
          hit = set()
          for row in near:
            hit.update(row)
          if not hit:
            continue
          idx = flex.size_t(sorted(hit))
          own = sites_cart.select(idx)
          dx, dy, dz = ((rot * own + off) - own).parts()
          keep = (flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)) >= tol_sq
          grown = set()
          for j_seq in idx.select(keep):
            grown.update(rg_of[int(j_seq)])
          if grown:
            by_key.setdefault((str(op), i, j, k), (rot, off, set()))[2].update(
              grown)
  hiers = []
  atoms = []
  xyz = flex.vec3_double()
  for key in sorted(by_key):
    rot, off, grp = by_key[key]
    # copy_atoms: without it set_xyz moves the model's own atoms.
    sub = hier.select(flex.size_t(sorted(grp)), copy_atoms=True)
    sub_atoms = sub.atoms()
    sub_atoms.set_xyz(rot * sub_atoms.extract_xyz() + off)
    hiers.append(sub)
    atoms.extend(list(sub_atoms))
    xyz.extend(sub_atoms.extract_xyz())
  return hiers, atoms, xyz


def _water_image_neighbours(sites, crystal_symmetry, radius,
                            min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Images of ``sites`` that crystal symmetry brings within ``radius`` of one.

  Enumerating the operators beats a pair table here: the sites are few, a
  table cannot be seeded on them, and the table would have to be built at this
  radius rather than the smaller one the atom environment needs.

  Returns ``(by_site, transforms, own)``: per site, the sorted ``(block, other
  site)`` pairs whose image is in range; per block the cartesian ``(rotation,
  translation)`` that produces it, block 0 being an unused placeholder so a
  block index is never zero; and per site the blocks whose operator brings the
  site's own image into range. An operator that fixes a site is dropped
  throughout, the equivalent a pair table drops for being coincident with its
  site.
  """
  n = sites.size()
  by_site = [[] for _ in range(n)]
  own = [[] for _ in range(n)]
  transforms = [None]
  if not n:
    return by_site, transforms, own
  unit_cell = crystal_symmetry.unit_cell()
  # Deposited coordinates need not lie inside one cell, so the translations
  # to try span the fractional spread of the sites, not just the radius.
  margin = [radius * x for x in unit_cell.reciprocal_parameters()[:3]]
  real = unit_cell.fractionalize(sites).parts()
  lo_real = [flex.min(c) for c in real]
  hi_real = [flex.max(c) for c in real]
  tree = KDTree(sites.as_numpy_array())
  tol_sq = min_distance_sym_equiv ** 2
  found = []
  found_own = []
  for op in crystal_symmetry.space_group():
    base = super_cell.sym_equiv_sites_cart(sites, unit_cell, op)
    frac = unit_cell.fractionalize(base).parts()
    spans = []
    for axis in range(3):
      lo = lo_real[axis] - flex.max(frac[axis]) - margin[axis]
      hi = hi_real[axis] - flex.min(frac[axis]) + margin[axis]
      spans.append(range(int(math.floor(lo)), int(math.ceil(hi)) + 1))
    for i in spans[0]:
      for j in spans[1]:
        for k in spans[2]:
          shift = unit_cell.orthogonalize((i, j, k))
          pts = base + shift
          dx, dy, dz = (pts - sites).parts()
          fixed = (flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)) < tol_sq
          for wj, near in enumerate(tree.query_ball_point(
              pts.as_numpy_array(), radius, return_sorted=False)):
            if fixed[wj]:
              continue
            for wi in near:
              if wi == wj:
                found_own.append((wj, (str(op), i, j, k)))
              else:
                found.append((int(wi), wj, (str(op), i, j, k)))
  if not (found or found_own):
    return by_site, transforms, own
  # Number the blocks by operator rather than by discovery order, so the slot
  # order a site sees does not depend on how the spatial query enumerated it.
  blocks = {}
  by_str = {str(op): op for op in crystal_symmetry.space_group()}
  keys = {f[2] for f in found} | {f[1] for f in found_own}
  for key in sorted(keys):
    blocks[key] = len(blocks) + 1
    op_str, i, j, k = key
    op = by_str[op_str]
    transforms.append((
      matrix.sqr(unit_cell.matrix_cart(op.r())),
      matrix.col(unit_cell.orthogonalize(op.t().as_double()))
      + matrix.col(unit_cell.orthogonalize((i, j, k)))))
  for wi, wj, key in found:
    by_site[wi].append((blocks[key], wj))
  for wj, key in found_own:
    own[wj].append(blocks[key])
  for pairs in by_site:
    pairs.sort()
  for pairs in own:
    pairs.sort()
  return by_site, transforms, own


class _WaterHydrogenPlacer(object):
  """Engine behind :func:`place_water_hydrogens`.

  Holds the placement state (static-atom KDTree, donor/acceptor bookkeeping,
  per-water neighbour blocks, placed protons) and the placement, refinement
  and basin-hopping logic.

  Only the placed water H move: the model and the water O do not, and every
  candidate H lies exactly ``oh_length`` from its O. Each water's static
  neighbour block and water-neighbour list are therefore built once in
  :meth:`run` and reused by the greedy pass, every relaxation sweep and every
  basin round.

  The constructor records the parameters; :meth:`run` does the work, in
  place. Parameter semantics: :func:`place_water_hydrogens`.
  """

  def __init__(self, hier, oh_length=None, element=None,
               n_refine=_WATER_REFINE_SWEEPS, refine_tol=_WATER_REFINE_TOL,
               n_basin=0, existing_h="keep", lone_pair_directed=False,
               on_state=None, crystal_symmetry=None,
               min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
    self.hier = hier
    self.oh_length = oh_length
    self.element = element
    self.n_refine = n_refine
    self.refine_tol = refine_tol
    self.n_basin = n_basin
    self.existing_h = existing_h
    self.lone_pair_directed = lone_pair_directed
    self.on_state = on_state
    self.crystal_symmetry = crystal_symmetry
    self.min_distance_sym_equiv = min_distance_sym_equiv

    # Placed-H coordinates, one slot per proton, filled in run(); placed_np
    # is the same data as an (n, 3) array.
    self.placed_coords = []
    self.placed_np = None
    # Sweep active-set stamps (see _dirty), sized in run(): off one
    # monotonic counter, the tick at which each water was last placed and the
    # tick at which its H last moved.
    self.w_placed_at = []
    self.w_moved_at = []
    self.tick = 0
    # (water index, [(atom, slot, di), ...], fixed_d1)
    self.records = []
    # (residue id, (metal element, distance) or None, action)
    self.partial_waters = []

  def _static_block(self, wi, cands):
    """Candidates against water ``wi``'s static neighbours.

    Returns ``(d, ok)``: the (k, m) squared-distance block and whether each
    candidate clears every one of those neighbours. ``d`` is None when the
    water has no static neighbour at all.
    """
    P = self.w_sxyz[wi]
    if not len(P):
      return None, np.ones(len(cands), dtype=bool)
    dx = cands[:, 0, None] - P[None, :, 0]
    dy = cands[:, 1, None] - P[None, :, 1]
    dz = cands[:, 2, None] - P[None, :, 2]
    d = dx * dx + dy * dy + dz * dz
    return d, ~(d < self.w_sthr[wi][None, :]).any(axis=1)

  def _static_clear(self, wi, cands, key):
    """Static half of the clearance test, reduced and cached per water.

    ``cands`` must be an invariant candidate set (the acceptor directions or
    the fallback sphere) and ``key`` the geometry-dict slot to cache under.
    Returns ``(ok, mins)``: the per-candidate flag and the per-candidate
    minimum squared distance, the latter None when the water has no static
    neighbour.
    """
    g = self.w_geom[wi]
    st = g[key]
    if st is None:
      d, ok = self._static_block(wi, cands)
      st = g[key] = (ok, None if d is None else d.min(axis=1))
    return st

  def _self_image_clear(self, wi, cands, own_fixed=()):
    """Candidates against the images of water ``wi``'s own protons.

    A water whose own image is in reach cannot be scored against a standing
    point set: moving a candidate moves its image with it. Each operator is
    applied to the whole candidate array instead, which gives every
    candidate's distance to its own image, and to the image of each proton in
    ``own_fixed``, this water's protons already settled this pass.

    The transform runs through flex: numpy sums the three terms of the matrix
    product in another order, and on a dense rotation the two part company in
    the last ulp, which is enough to flip a threshold test.

    Returns ``(ok, mins)``, both None when no operator brings this water's own
    image into range.
    """
    ops = self.w_self_ops[wi]
    if not ops:
      return None, None
    ok = np.ones(len(cands), dtype=bool)
    best = None
    thr = _WATER_MIN_H_CLEARANCE ** 2
    for rot, trn in ops:
      img = (rot * flex.vec3_double(cands) + trn).as_numpy_array()
      pairs = [(cands, img)]
      for pt in own_fixed:
        one = flex.vec3_double(np.asarray(pt, dtype=float).reshape(1, 3))
        pairs.append((cands, (rot * one + trn).as_numpy_array()[0]))
        pairs.append((img, np.asarray(pt, dtype=float)))
      for left, right in pairs:
        dx = left[..., 0] - right[..., 0]
        dy = left[..., 1] - right[..., 1]
        dz = left[..., 2] - right[..., 2]
        d = dx * dx + dy * dy + dz * dz
        ok &= ~(d < thr)
        best = d if best is None else np.minimum(best, d)
    return ok, best

  def _clear(self, wi, cands, nbr_slots, static_key=None, own_fixed=()):
    """Clearance test over the candidate H positions of one water.

    Every candidate lies exactly ``oh_length`` from the water O, so water
    ``wi``'s cached static block (a ball of
    ``oh_length + _WATER_CLEARANCE_RADIUS`` about the O, own atoms already
    dropped) covers every candidate's own clearance ball. The test is one dense
    (candidate x neighbour) block.

    ``nbr_slots`` are the placed-H slots that may lie near this water, its own
    excluded. ``static_key`` names the geometry-dict slot holding the cached
    static half, for a candidate set that does not move; None for the cone,
    which is rebuilt around each pass's O-H1 axis. ``own_fixed`` are this
    water's protons already settled this pass, which only matter to a water
    whose own image is in reach (see :meth:`_self_image_clear`).

    Returns ``(min_dist, ok)`` per candidate: the distance to the nearest
    non-own atom within ``_WATER_CLEARANCE_RADIUS`` (the search radius if
    nothing is near, a soft tie-breaker) and whether the candidate clears every
    heavy atom by ``_WATER_MIN_CLEARANCE`` and every hydrogen (static or
    placed) by ``_WATER_MIN_H_CLEARANCE``.
    """
    # Squared distances throughout: the thresholds and the cap are exact
    # squares and sqrt is monotone.
    best = np.full(len(cands), _WATER_CLEARANCE_RADIUS ** 2)
    if static_key is not None:
      ok, mins = self._static_clear(wi, cands, static_key)
      ok = ok.copy()
      if mins is not None:
        np.minimum(best, mins, out=best)
    else:
      d, ok = self._static_block(wi, cands)
      if d is not None:
        np.minimum(best, d.min(axis=1), out=best)
    if len(nbr_slots):
      Q = self.pool_np[nbr_slots]
      dx = cands[:, 0, None] - Q[None, :, 0]
      dy = cands[:, 1, None] - Q[None, :, 1]
      dz = cands[:, 2, None] - Q[None, :, 2]
      d = dx * dx + dy * dy + dz * dz
      np.minimum(best, d.min(axis=1), out=best)
      ok &= ~(d < _WATER_MIN_H_CLEARANCE ** 2).any(axis=1)
    s_ok, s_mins = self._self_image_clear(wi, cands, own_fixed)
    if s_ok is not None:
      ok &= s_ok
      np.minimum(best, s_mins, out=best)
    return np.sqrt(best), ok

  @staticmethod
  def _cat_ok(cat, pts, o):
    """Whether each candidate is clear of every nearby cation.

    True where ``(H - O).(cation - O) <= 0`` for every O->cation direction in
    ``cat``, with ``pts`` the candidate H positions and ``o`` the water O.
    """
    if not len(cat):
      return np.ones(len(pts), dtype=bool)
    dx = pts[:, 0, None] - o[0]
    dy = pts[:, 1, None] - o[1]
    dz = pts[:, 2, None] - o[2]
    dp = dx * cat[None, :, 0] + dy * cat[None, :, 1] + dz * cat[None, :, 2]
    return (dp <= 0.0).all(axis=1)

  def _stats(self):
    """:func:`_water_clash_stats` from the engine's own arrays."""
    ns = len(self.placed_coords)
    if ns:
      X = np.concatenate((self.placed_np[:ns], self.wh_xyz))
      W = np.concatenate((self.slot_wid[:ns], self.wh_wid))
    else:
      X, W = self.wh_xyz, self.wh_wid
    return _contact_stats(
      len(X), _water_h_contacts(flex.vec3_double(X), flex.size_t(W))[2])

  def _nearest_cation(self, o_xyz, own_idx):
    """Closest metal cation coordinating a water O, if any.

    Returns ``(element, distance)`` or None, over
    ``_WATER_METAL_COORD_RADIUS`` (a first-shell bond) and excluding the
    water's own atoms ``own_idx``. Reporting only; placement ignores it.
    """
    o = matrix.col(o_xyz)
    best = None
    for i in self.static_tree.query_ball_point(tuple(o_xyz),
                                               _WATER_METAL_COORD_RADIUS):
      if i in own_idx:
        continue
      el = self.atoms[i].element.strip().upper()
      if el not in _WATER_CATION_ELEMENTS:
        continue
      d = (matrix.col(self.atoms[i].xyz) - o).length()
      if best is None or d < best[1]:
        best = (el, d)
    return best

  def _water_geom(self, wi):
    """Fixed per-water geometry, computed once.

    Returns the dict of ``o`` (O coordinates), ``acceptors`` (atom indices,
    nearest first), ``acc_dirs``/``acc_pts`` (O-H unit directions and the H
    positions built from them), ``cat`` (O->cation unit directions), ``sph``
    (the lazily built fallback-sphere H positions) and
    ``acc_static``/``sph_static`` (the cached static half of the clearance test
    for those two candidate sets, see :meth:`_static_clear`).
    """
    g = self.w_geom[wi]
    if g is not None:
      return g
    atoms = self.atoms
    o_xyz = self.w_o_col[wi]
    own_idx = self.w_own[wi]
    # Acceptors (static O/N/F/S/Cl) within range, nearest first.
    acceptors = []
    for i in self.w_acc_raw[wi]:
      if i in own_idx:
        continue
      if i in self.donor_n:
        continue  # protonated N: a donor, not a usable acceptor
      a = atoms[i]
      if a.element.strip().upper() not in _WATER_ACCEPTOR_ELEMENTS:
        continue
      d = (matrix.col(a.xyz) - o_xyz).length()
      if d < 1e-3:
        continue  # atom coincident with O (alt-conf / overlap): no direction
      acceptors.append((d, i))
    acceptors.sort(key=lambda t: t[0])
    acceptors = [i for _, i in acceptors]

    def accept_dir(i):
      """Unit direction from the water O toward acceptor ``i``.

      A lone-pair lobe facing the water when lone-pair-directed placement is
      on, else the acceptor nucleus.
      """
      a_xyz = matrix.col(atoms[i].xyz)
      if self.lone_pair_directed:
        toward = (o_xyz - a_xyz).normalize()
        best = None
        for lobe in self.acc_lobes.get(i, ()):
          if best is None or lobe.dot(toward) > best.dot(toward):
            best = lobe
        if best is not None and best.dot(toward) > 0.0:
          return (a_xyz + best * _WATER_HBOND_HA - o_xyz).normalize()
      return (a_xyz - o_xyz).normalize()

    # Nearby metal cations, as O->cation unit directions (see _cat_ok).
    cation_dirs = []
    for i in self.w_cat_raw[wi]:
      if i in own_idx:
        continue
      a = atoms[i]
      if a.element.strip().upper() not in _WATER_CATION_ELEMENTS:
        continue
      v = matrix.col(a.xyz) - o_xyz
      if v.length() < 1e-3:
        continue
      cation_dirs.append(v.normalize())

    o = np.array(o_xyz.elems, dtype=float)
    acc_dirs = _as_np([accept_dir(i) for i in acceptors])
    g = {"o": o,
         "acceptors": acceptors,
         "acc_dirs": acc_dirs,
         "acc_pts": o + self.oh_length * acc_dirs,
         "cat": _as_np(cation_dirs),
         "sph": None,
         "acc_static": None,
         "sph_static": None}
    self.w_geom[wi] = g
    return g

  def _place_one(self, wi, nbr_slots, fixed_d1=None):
    """The two H positions for one water O, clash-aware.

    ``nbr_slots`` are the placed-H slots that may lie near this water, its own
    excluded. ``fixed_d1`` holds the unit O-H1 direction fixed instead of
    searching for one, for a water that already carries a proton: the returned
    ``h1_xyz`` then only restates it and just ``h2_xyz`` is new.

    Returns ``(h1_xyz, h2_xyz)``.
    """
    g = self._water_geom(wi)
    o = g["o"]
    acceptors = g["acceptors"]
    acc_dirs = g["acc_dirs"]
    acc_pts = g["acc_pts"]
    cat = g["cat"]
    na = len(acceptors)
    if na:
      acc_best, acc_ok = self._clear(wi, acc_pts, nbr_slots, "acc_static")
      acc_cat = self._cat_ok(cat, acc_pts, o)

    # H1: nearest acceptor giving a placement that is away from cations and
    # clash-free; else the best direction over the acceptors plus a dense
    # fallback sphere, ranked (away-from-cation, clash-free, clearance).
    # ``h1_k`` is H1's position in the acceptor list, or -1.
    if fixed_d1 is not None:
      # Find the acceptor the deposited H1 donates to, so H2 aims elsewhere.
      d1 = fixed_d1
      h1_k = -1
      best_align = _WATER_H1_ACC_ALIGN
      for k in range(na):
        align = (acc_dirs[k, 0] * d1[0] + acc_dirs[k, 1] * d1[1]
                 + acc_dirs[k, 2] * d1[2])
        if align > best_align:
          best_align = align
          h1_k = k
    else:
      d1 = None
      h1_k = -1
      for k in range(na):
        if acc_ok[k] and acc_cat[k]:  # ok and away from cations
          d1 = acc_dirs[k]
          h1_k = k
          break
      if d1 is None:
        if g["sph"] is None:
          g["sph"] = o + self.oh_length * _FALLBACK_NP
        sph_pts = g["sph"]
        s_best, s_ok = self._clear(wi, sph_pts, nbr_slots, "sph_static")
        s_cat = self._cat_ok(cat, sph_pts, o)
        if na:
          cand_dirs = np.concatenate((acc_dirs, _FALLBACK_NP))
          cand_best = np.concatenate((acc_best, s_best))
          cand_ok = np.concatenate((acc_ok, s_ok))
          cand_cat = np.concatenate((acc_cat, s_cat))
        else:
          cand_dirs, cand_best, cand_ok, cand_cat = (
            _FALLBACK_NP, s_best, s_ok, s_cat)
        # (cation-ok, clash-free) as one rank, then clearance, first wins.
        rank = cand_cat.astype(np.int8) * 2 + cand_ok
        top = rank == rank.max()
        d1 = cand_dirs[int(np.argmax(np.where(top, cand_best, -np.inf)))]
    h1_xyz = o + self.oh_length * d1

    p, q = _ortho_frame(d1)

    # H2: rank the cone angles by (away-from-cation, clash-free, best
    # alignment over the acceptors H1 did not take, clearance).
    cone_dirs = (self.cos_hoh * d1
                 + self.sin_hoh * (_CONE_COS[:, None] * p
                                   + _CONE_SIN[:, None] * q))
    cone_pts = o + self.oh_length * cone_dirs
    c_best, c_ok = self._clear(wi, cone_pts, nbr_slots,
                               own_fixed=(h1_xyz,))
    c_cat = self._cat_ok(cat, cone_pts, o)

    if na - (h1_k >= 0):
      dots = (cone_dirs[:, 0, None] * acc_dirs[None, :, 0]
              + cone_dirs[:, 1, None] * acc_dirs[None, :, 1]
              + cone_dirs[:, 2, None] * acc_dirs[None, :, 2])
      if h1_k >= 0:
        dots[:, h1_k] = -np.inf
      align = dots.max(axis=1)
    else:
      align = np.zeros(_WATER_CONE_SAMPLES)
    rank = c_cat.astype(np.int8) * 2 + c_ok
    top = rank == rank.max()
    a3 = np.where(c_cat & c_ok, align, 0.0)
    top &= a3 == np.where(top, a3, -np.inf).max()
    return h1_xyz, cone_pts[int(np.argmax(np.where(top, c_best, -np.inf)))]

  def _pool_write(self, wi, slot, xyz):
    """Put one proton of water ``wi`` at ``xyz``, its images with it."""
    self.placed_np[slot] = xyz
    if self.pool_np is self.placed_np:
      return
    self.pool_np[slot] = xyz
    v = matrix.col(xyz)
    for b, rot, trn in self.w_blocks[wi]:
      self.pool_np[b * self.n_slots + slot] = (rot * v + trn).elems

  def _store(self, wi, slots, h1, h2):
    """Write one water's new H positions to the model and the arrays.

    True if any of them moved, which is what :meth:`_dirty` reads.
    """
    moved = False
    for atom, slot, di in slots:
      xyz = tuple((h1 if di == 1 else h2).tolist())
      if xyz != self.placed_coords[slot]:
        moved = True
      self.placed_coords[slot] = xyz
      self._pool_write(wi, slot, xyz)
      atom.set_xyz(xyz)
    return moved

  def _dirty(self, wi):
    """Whether re-placing water ``wi`` could move its H.

    A water is placed against its static surroundings, which never move, and
    the placed H of the waters in ``w_wnbr[wi]`` and the images in
    ``w_inbr[wi]``, an image moving exactly when its own water does. If none
    of those H has moved since this water was last placed, the placement
    re-derives the two positions it already holds. The monotonic tick is bumped once per water per
    sweep, so a neighbour that moves earlier in the same sweep still counts.
    """
    t = self.w_placed_at[wi]
    if t < 0:
      return True   # never placed against the finished set of water H
    moved_at = self.w_moved_at
    if moved_at[wi] > t:
      return True   # moved since it was last placed, i.e. kicked
    for wj in self.w_wnbr[wi]:
      if moved_at[wj] > t:
        return True
    for _b, wj in self.w_inbr[wi]:
      if moved_at[wj] > t:
        return True
    return False

  def _apply_sweep(self):
    """Run one relaxation sweep.

    Re-places every recorded water whose neighbourhood has changed (see
    :meth:`_dirty`) against the current positions of all the others, its own H
    excluded, updating ``placed_coords`` and the atoms in place. A water placed
    later in the sweep sees the H the earlier ones just moved.
    """
    for wi, slots, fixed_d1 in self.records:
      self.tick += 1
      if self._dirty(wi):
        if self._store(wi, slots,
                       *self._place_one(wi, self.w_nbr_slots[wi], fixed_d1)):
          self.w_moved_at[wi] = self.tick
      self.w_placed_at[wi] = self.tick

  def _snapshot(self):
    """Capture the placed-H state as ``(coords, placed_at, moved_at)``.

    Whether a water may be skipped is a function of the coordinates and the
    stamps together, so the three are captured and restored as one.
    """
    return (list(self.placed_coords), list(self.w_placed_at),
            list(self.w_moved_at))

  def _restore(self, snap):
    """Reset all placed H to a snapshot from :meth:`_snapshot`."""
    coords, placed_at, moved_at = snap
    for wi, slots, _fixed in self.records:
      for atom, slot, di in slots:
        self.placed_coords[slot] = coords[slot]
        self._pool_write(wi, slot, coords[slot])
        atom.set_xyz(coords[slot])
    self.w_placed_at = list(placed_at)
    self.w_moved_at = list(moved_at)

  def _clashing_records(self):
    """Indices into ``records`` of waters with a placed H that clashes."""
    bad = []
    for ri, (wi, slots, _fixed) in enumerate(self.records):
      pts = self.placed_np[[slot for _, slot, _ in slots]]
      if not self._clear(wi, pts, self.w_nbr_slots[wi],
                         own_fixed=tuple(pts))[1].all():
        bad.append(ri)
    return bad

  def _kick(self, ri, rng):
    """Re-orient the water ``records[ri]`` randomly, from seeded ``rng``.

    H1 goes along a random axis and H2 to a random azimuth on the H-O-H cone.
    A water being completed has a fixed O-H1, so only the azimuth is random.
    """
    wi, slots, fixed_d1 = self.records[ri]
    o = self._water_geom(wi)["o"]
    d1 = fixed_d1 if fixed_d1 is not None else _rand_unit(rng)
    p, q = _ortho_frame(d1)
    theta = rng.uniform(0.0, 2.0 * math.pi)
    d2 = self.cos_hoh * d1 + self.sin_hoh * (math.cos(theta) * p
                                             + math.sin(theta) * q)
    self.tick += 1
    if self._store(wi, slots, o + self.oh_length * d1,
                   o + self.oh_length * d2):
      self.w_moved_at[wi] = self.tick

  def run(self):
    """Place H on every bare water, refine, optionally basin-hop, keep best.

    Modifies the hierarchy in place. Returns the kept-state label (see
    :func:`place_water_hydrogens`).
    """
    hier = self.hier
    assert hier.models_size() <= 1, (
      f"place_water_hydrogens takes one model, not {hier.models_size()}")
    # Every test below reads the element column.
    check_for_missing_elements(hier)
    # Resolve before stripping, which would remove the D this keys on.
    if self.oh_length is None:
      self.oh_length = (_WATER_OH_NEUTRON if _has_deuterium(hier)
                        else _WATER_OH_XRAY)
    # What each water carried is read off the strip; the walk below sees the
    # same protons in every other mode.
    single_h = None
    stripped = {}
    if self.existing_h == "reorient":
      stripped = _strip_water_hydrogens(hier)
      single_h = {k for k, els in stripped.items() if len(els) == 1}

    sel = hier.atoms()
    atoms = list(sel)
    if not atoms:
      return None
    sel.reset_i_seq()

    # Neighbouring asymmetric units, as environment atoms appended after the
    # model's own. They keep the asymmetric unit's indices 0..n-1 valid as
    # both tree and i_seq, and the water walk below reads the hierarchy, so
    # only the model's waters are protonated.
    self.sym_hier = []
    sym_atoms = []
    sym_xyz = flex.vec3_double()
    if self.crystal_symmetry is not None:
      self.sym_hier, sym_atoms, sym_xyz = _symmetry_environment(
        hier, sel.extract_xyz(), self.crystal_symmetry,
        max(self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
            _WATER_ACCEPTOR_RADIUS), self.min_distance_sym_equiv)
    atoms = atoms + sym_atoms
    self.atoms = atoms

    # Static neighbours (protein, ligands, water O, pre-existing H) never
    # move; the placed water H are tracked by slot in placed_coords/placed_np.
    self.static_np = sel.extract_xyz().as_numpy_array()
    if sym_atoms:
      self.static_np = np.vstack([self.static_np, sym_xyz.as_numpy_array()])
    self.static_tree = KDTree(self.static_np)
    _el = sel.extract_element(strip=True)
    _el.extend(flex.std_string([a.element.strip() for a in sym_atoms]))
    self.static_is_h = ((_el == "H") | (_el == "D")).as_numpy_array()
    self.static_thr = np.where(self.static_is_h, _WATER_MIN_H_CLEARANCE ** 2,
                               _WATER_MIN_CLEARANCE ** 2)

    # N atoms that already carry an H are donors, not acceptors (amide,
    # ammonium, guanidinium, protonated His ring N, ...). O always accepts,
    # so only N is filtered.
    self.donor_n = set()
    h_idx = np.nonzero(self.static_is_h)[0]
    if len(h_idx):
      for i, nbrs in zip(h_idx, self.static_tree.query_ball_point(
          self.static_np[h_idx], _WATER_NH_BOND,
          return_sorted=False)):
        for j in nbrs:
          if j != i and atoms[j].element.strip().upper() == "N":
            self.donor_n.add(j)

    # Lone-pair lobe directions per acceptor (opt-in; empty when off).
    self.acc_lobes = _acceptor_lobes(atoms, self.static_tree, self.donor_n) \
        if self.lone_pair_directed else {}

    self.cos_hoh = math.cos(math.radians(_WATER_HOH_DEG))
    self.sin_hoh = math.sin(math.radians(_WATER_HOH_DEG))

    # Gather the waters to protonate, in one walk per water residue.
    waters = []
    wh_xyz = []
    wh_wid = []
    for wgid, ag in enumerate(g for g in hier.atom_groups()
                              if _is_water(g.resname)):
      o = None
      existing = []
      names = set()
      own_idx = set()
      ats, hd = _hd_flags(ag)
      els = ats.extract_element(strip=True)
      for k in range(ats.size()):
        a = ats[k]
        names.add(a.name.strip())
        own_idx.add(a.i_seq)
        if hd[k]:
          existing.append(a)
          wh_xyz.append(a.xyz)
          wh_wid.append(wgid)
        elif o is None and els[k].upper() == "O":
          o = a
      if o is None:
        continue
      fixed_d1 = None
      skip = len(existing) >= 2          # already protonated
      if existing and not skip:
        if self.existing_h != "complete":
          skip = True
        else:
          v = matrix.col(existing[0].xyz) - matrix.col(o.xyz)
          if v.length() < 1e-3:
            skip = True  # H coincident with O: no direction for a cone
          else:
            fixed_d1 = np.array(v.normalize().elems, dtype=float)
      is_single = (ag.memory_id() in single_h if single_h is not None
                   else len(existing) == 1)
      if skip and not is_single:
        continue                         # nothing to place and nothing to say
      if is_single:
        action = ("stripped" if self.existing_h == "reorient"
                  else "completed" if fixed_d1 is not None else "kept")
        self.partial_waters.append(
          (_water_id(ag), self._nearest_cation(o.xyz, own_idx), action))
      if skip:
        continue
      waters.append((ag, o, own_idx, fixed_d1, existing, names, wgid))
    self.wh_xyz = np.array(wh_xyz, dtype=float) if wh_xyz else np.zeros((0, 3))
    self.wh_wid = np.array(wh_wid, dtype=np.int64)

    # Per-water constants: every neighbour list the placement needs, built
    # once here and reused by the greedy pass, every relaxation sweep and
    # every basin round.
    n = len(waters)
    self.w_geom = [None] * n
    self.placed_coords = []
    self.placed_np = np.zeros((2 * n, 3))
    # Protons the clearance test may draw on: the placed ones, plus one
    # transformed block per symmetry operator that brings another water's
    # protons into reach. Block b slot s lives at b * n_slots + s, so an image
    # proton is just another slot. Without images the pool is the placed array
    # itself and nothing extra is paid.
    self.n_slots = 2 * n
    self.pool_np = self.placed_np
    self.img_tf = [None]
    self.w_inbr = [[] for _ in range(n)]
    self.w_blocks = [()] * n
    # Cartesian (rotation, translation) per operator that brings the water's
    # own image into range.
    self.w_self_ops = [()] * n
    self.slot_wid = np.zeros(2 * n, dtype=np.int64)
    self.records = []   # (water index, [(atom, slot, di), ...], fixed_d1)
    self.w_wnbr = []
    # Every water is dirty for the first sweep.
    self.w_placed_at = [-1] * n
    self.w_moved_at = [0] * n
    self.tick = 0
    if n:
      o_pts = np.array([w[1].xyz for w in waters], dtype=float)
      # Order most-crowded first. ``crowd`` is the neighbour count within
      # the clash radius, the water's own atoms excluded.
      crowd = [sum(1 for j in nb if j not in waters[k][2]) for k, nb in
               enumerate(self.static_tree.query_ball_point(
                 o_pts, _WATER_CLEARANCE_RADIUS,
                 return_sorted=False))]
      order = sorted(range(n), key=lambda k: crowd[k], reverse=True)
      waters = [waters[k] for k in order]
      o_pts = o_pts[order]
      self.w_own = [w[2] for w in waters]
      self.w_o_col = [matrix.col(w[1].xyz) for w in waters]
      # query_ball_point sorts a batch's indices but not a single point's, and
      # the acceptor order is that sort's tie-break: the batch must not sort.
      self.w_acc_raw = self.static_tree.query_ball_point(
        o_pts, _WATER_ACCEPTOR_RADIUS, return_sorted=False)
      self.w_cat_raw = self.static_tree.query_ball_point(
        o_pts, _WATER_CATION_RADIUS, return_sorted=False)
      # One static-neighbour block per water: every candidate H lies on the
      # O-H sphere about the O, so a single ball of oh_length + clearance
      # covers every candidate's own clearance ball.
      self.w_sxyz = []
      self.w_sthr = []
      for wi, nb in enumerate(self.static_tree.query_ball_point(
          o_pts, self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
          return_sorted=False)):
        own = self.w_own[wi]
        idx = np.array([j for j in nb if j not in own], dtype=np.intp)
        self.w_sxyz.append(self.static_np[idx])
        self.w_sthr.append(self.static_thr[idx])
      # A placed H sits within oh_length of its own O, so only waters whose
      # O lie within clearance + 2 oh_length can hold one near this water's
      # candidates.
      r_wh = _WATER_CLEARANCE_RADIUS + 2.0 * self.oh_length + 0.01
      self.w_wnbr = [[j for j in nb if j != wi] for wi, nb in
                     enumerate(KDTree(o_pts).query_ball_point(
                       o_pts, r_wh, return_sorted=False))]
      # The same test against the images of these waters.
      if self.crystal_symmetry is not None:
        self.w_inbr, self.img_tf, own_blocks = _water_image_neighbours(
          flex.vec3_double(o_pts), self.crystal_symmetry, r_wh,
          self.min_distance_sym_equiv)
        # matrix_cart's raw tuple is what multiplies a vec3_double array;
        # the matrix.sqr wrapper _pool_write uses does not.
        self.w_self_ops = [
          tuple((self.img_tf[b][0].elems, self.img_tf[b][1].elems)
                for b in own_blocks[wi]) for wi in range(n)]
        if len(self.img_tf) > 1:
          in_block = [set() for _ in range(n)]
          for wi in range(n):
            for b, wj in self.w_inbr[wi]:
              in_block[wj].add(b)
          self.w_blocks = [tuple((b,) + self.img_tf[b] for b in sorted(bs))
                           for bs in in_block]
          self.pool_np = np.zeros((len(self.img_tf) * self.n_slots, 3))

    # Initial greedy pass over the ordered waters, each avoiding the H
    # already placed on earlier ones. ``records`` keeps per-water (index,
    # placed-H slots, fixed_d1) for the refinement sweeps.
    slots_of = [()] * n

    def nbr_slots(wi):
      """Pool slots of every water and image neighbouring water ``wi``."""
      nbr = [s for wj in self.w_wnbr[wi] for s in slots_of[wj]]
      for b, wj in self.w_inbr[wi]:
        off = b * self.n_slots
        nbr.extend(off + s for s in slots_of[wj])
      return np.array(nbr, dtype=np.intp) if nbr else _EMPTY_SLOTS

    for wi in range(n):
      ag, o, own_idx, fixed_d1, existing, existing_names, wgid = waters[wi]
      if self.element is not None:
        proton_element = self.element
      elif existing:
        proton_element = existing[0].element.strip().upper()
      elif ag.memory_id() in stripped:
        proton_element = stripped[ag.memory_id()][0]
      else:
        proton_element = "D" if ag.resname.strip().upper() == "DOD" else "H"
      h1, h2 = self._place_one(wi, nbr_slots(wi), fixed_d1)

      # Completing builds only the cone proton; the name takes the free slot.
      slots = []
      for di in ((2,) if fixed_d1 is not None else (1, 2)):
        proton_name = _free_proton_name(existing_names, proton_element)
        if proton_name is None:
          continue
        existing_names.add(proton_name.strip())
        xyz = tuple((h1 if di == 1 else h2).tolist())
        atom = _new_h_atom(proton_name, proton_element, xyz, o.occ, o.b,
                           o.hetero)
        ag.append_atom(atom)
        slot = len(self.placed_coords)
        slots.append((atom, slot, di))
        self.placed_coords.append(xyz)
        self._pool_write(wi, slot, xyz)
        self.slot_wid[slot] = wgid
      if slots:
        self.records.append((wi, slots, fixed_d1))
        slots_of[wi] = tuple(s for _, s, _ in slots)

    # The placed-H neighbourhood of every water, now that all of them exist.
    self.w_nbr_slots = [nbr_slots(wi) for wi in range(n)]

    # Refinement sweeps: re-place each water against the final set of all
    # placed H, its own excluded. With refine_tol > 0 this stops early once a
    # sweep reduces the close (<2.0 A) contact count by fewer than refine_tol
    # (n_refine is then a cap); refine_tol == 0 runs all n_refine sweeps. The
    # best state is kept and restored at the end.
    stats = self._stats() \
        if (self.on_state or self.n_refine or self.n_basin) else None
    if self.on_state is not None:
      self.on_state("initial", stats)
    prev_n20 = stats[1] if stats is not None else None
    best_n20 = prev_n20
    best_coords = self._snapshot() if (self.n_refine or self.n_basin) else None
    best_label = "initial"
    for i in range(self.n_refine):
      self._apply_sweep()
      stats = self._stats()
      if self.on_state is not None:
        self.on_state(f"sweep {i + 1}", stats)
      if best_n20 is None or stats[1] < best_n20:
        best_n20 = stats[1]
        best_coords = self._snapshot()
        best_label = f"sweep {i + 1}"
      if self.refine_tol and prev_n20 is not None and prev_n20 - stats[1] < self.refine_tol:
        break  # gain below tolerance: converged
      prev_n20 = stats[1]

    # Basin-hopping (optional): each round restarts from the best state,
    # kicks the still-clashing waters to a random orientation, relaxes, and
    # keeps the result if it improved. Deterministic (seeded).
    if self.n_basin and best_coords is not None:
      rng = random.Random(_WATER_BASIN_SEED)
      for it in range(self.n_basin):
        self._restore(best_coords)
        offenders = self._clashing_records()
        if not offenders:
          break  # nothing left to relax
        for ri in offenders:
          self._kick(ri, rng)
        for s in range(_WATER_BASIN_RELAX):
          self._apply_sweep()
          stats = self._stats()
          label = f"basin {it + 1}.{s + 1}"
          if self.on_state is not None:
            self.on_state(label, stats)
          if best_n20 is None or stats[1] < best_n20:
            best_n20 = stats[1]
            best_coords = self._snapshot()
            best_label = label

    if best_coords is not None:
      self._restore(best_coords)
    return best_label if (self.n_refine or self.n_basin) else None


def place_water_hydrogens(hier, oh_length=None, element=None,
                          n_refine=_WATER_REFINE_SWEEPS,
                          refine_tol=_WATER_REFINE_TOL, n_basin=0,
                          existing_h="keep", lone_pair_directed=False,
                          on_state=None, crystal_symmetry=None,
                          min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Place the two H on every bare water, H-bond-aware.

  For each water residue missing H (any common water alias: HOH, DOD, H2O,
  WAT, OH2, ...): H1 along O -> the nearest acceptor giving a clash-free H
  (``_WATER_ACCEPTOR_RADIUS``, ``_WATER_ACCEPTOR_ELEMENTS``, N carrying an H
  excluded as donors), else the max-clearance direction over a dense sphere;
  H2 on the ``_WATER_HOH_DEG`` cone about O-H1, at a clash-free angle toward
  a second acceptor when there is one, else the clearest. Candidates are
  scored against every other atom and every other water's placed H; waters go
  most-crowded first and new H inherit the parent O's occupancy and B
  factor.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Single-model hierarchy; modified in place. Every atom must carry an
      element symbol.
  oh_length : float or None, optional
      O-H bond length in A, positive. None (default) picks
      ``_WATER_OH_NEUTRON`` (0.984) if the model contains D, else 0.957.
  element : str or None, optional
      Element of the placed atoms, ``"H"`` or ``"D"``, forced on every water;
      None (default) takes the element of the water's own H/D, those
      ``existing_h="reorient"`` strips included, else ``"D"`` for DOD and
      ``"H"`` for HOH. Names follow it: ``H1``/``H2`` or ``D1``/``D2``.
  n_refine : int, optional
      Maximum relaxation sweeps after the greedy pass, each re-placing every
      water against the final environment, best state kept; 0 disables it.
  refine_tol : int, optional
      Stop refining once a sweep removes fewer than this many close (<2.0 A)
      H-H contacts, making ``n_refine`` a cap; 0 runs all ``n_refine`` sweeps.
  n_basin : int, optional
      Basin-hopping rounds after refinement (default 0 = off): each restarts
      from the best state, randomly re-orients the still-clashing waters,
      relaxes, and keeps any improvement (deterministic, seeded).
  existing_h : str, optional
      What to do with a water that already carries H: ``"keep"`` (default,
      leave it untouched, single-H waters included), ``"complete"`` (a water
      with exactly one H gets its partner on that proton's cone, inheriting
      its element) or ``"reorient"`` (strip all water H and re-place both).
  lone_pair_directed : bool, optional
      If True, each O-H aims at an acceptor's lone-pair lobe (estimated from
      its bonded-neighbour geometry) rather than its nucleus.
  on_state : callable, optional
      Called as ``on_state(label, stats)`` at each state reached:
      ``"initial"`` after the greedy pass, ``"sweep N"`` after each refinement
      sweep, ``"basin N.M"`` during basin-hopping; ``stats`` is the
      ``_water_clash_stats`` tuple.
  crystal_symmetry : cctbx.crystal.symmetry or None, optional
      Honour crystal packing: atoms that symmetry places within reach of a
      water join its environment, so H at a lattice contact avoid the
      neighbouring asymmetric units instead of pointing into them. None
      (default) treats the model as isolated. A symmetry mate contributes its
      O and its placed protons both.
  min_distance_sym_equiv : float, optional
      Distance in A under which a symmetry equivalent counts as coincident
      with its own site, and so as that same atom rather than a second copy
      (default 0.5). A water refined a little off a symmetry element needs a
      larger value to be recognised as sitting on it.

  Returns
  -------
  libtbx.group_args
      ``kept_label``: label of the kept state (``"initial"``, ``"sweep N"`` or
      ``"basin N.M"``) when refinement or basin-hopping ran, else None.
      ``partial_waters``: one ``(residue_id, metal, action)`` per water that
      carried exactly one H on input, ``metal`` being ``(element, distance)``
      for a coordinating cation or None and ``action`` one of ``"kept"``,
      ``"completed"``, ``"stripped"``.
  """
  placer = _WaterHydrogenPlacer(
    hier, oh_length=oh_length, element=element, n_refine=n_refine,
    refine_tol=refine_tol, n_basin=n_basin,
    existing_h=existing_h,
    lone_pair_directed=lone_pair_directed,
    on_state=on_state, crystal_symmetry=crystal_symmetry,
    min_distance_sym_equiv=min_distance_sym_equiv)
  kept_label = placer.run()
  return group_args(kept_label=kept_label,
                    partial_waters=placer.partial_waters)


def _water_h_contacts(xyz, wid):
  """Inter-water H-H contacts within 2.0 A.

  ``xyz`` are water H sites and ``wid`` the water each belongs to. Returns
  ``(i, j, d)``: the two H of each contact, H on the same water excluded, and
  their distance.
  """
  empty = (flex.size_t(), flex.size_t(), flex.double())
  if xyz.size() < 2:
    return empty
  pairs = KDTree(xyz.as_numpy_array()).query_pairs(2.0, output_type="ndarray")
  if not len(pairs):
    return empty
  # KDTree hands back an (n, 2) ndarray; flex flattens it row-major, so the
  # pair members are the even and odd entries.
  p = flex.size_t(pairs)
  i = p.select(flex.size_t_range(0, p.size(), 2))
  j = p.select(flex.size_t_range(1, p.size(), 2))
  keep = (wid.select(i) != wid.select(j)).iselection()
  i = i.select(keep)
  j = j.select(keep)
  dx, dy, dz = (xyz.select(i) - xyz.select(j)).parts()
  return i, j, flex.sqrt(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz))


def _contact_stats(n, d):
  """``(n, n_lt_20, n_lt_18, n_lt_15, closest)`` for ``n`` water H with
  contact distances ``d``; ``closest`` is None without a contact."""
  if not d.size():
    return n, 0, 0, 0, None
  return (n, d.size(), (d < 1.8).count(True), (d < 1.5).count(True),
          flex.min(d))


def _water_h_sites(hier):
  """``(xyz, wid, atoms)`` for every water H/D in ``hier``, ``wid`` numbering
  the water each belongs to."""
  xyz = flex.vec3_double()
  wid = flex.size_t()
  atoms = []
  for w, ag in enumerate(g for g in hier.atom_groups()
                         if _is_water(g.resname)):
    ats, hd = _hd_flags(ag)
    sel = hd.iselection()
    xyz.extend(ats.extract_xyz().select(sel))
    wid.extend(flex.size_t(sel.size(), w))
    atoms.extend(ats[int(k)] for k in sel)
  return xyz, wid, atoms


def _water_clash_stats(hier):
  """Count H-H contacts between the H of different waters.

  Returns ``(n_placed, n_lt_20, n_lt_18, n_lt_15, closest)``: the number of
  water H, the counts of inter-water H-H contacts below 2.0/1.8/1.5 A, and
  the closest such distance (None if no pair is within 2.0 A).
  """
  xyz, wid, _ = _water_h_sites(hier)
  return _contact_stats(xyz.size(), _water_h_contacts(xyz, wid)[2])


def _clash_row(label, stats, log):
  """Print one row of the per-sweep clash table for ``stats`` to ``log``."""
  _, n20, n18, n15, worst = stats
  w = f"{worst:.2f}" if worst is not None else ">2.0"
  print(f"  {label:<9} <2.0={n20:<5} <1.8={n18:<5} <1.5={n15:<5} closest={w}",
        file=log)


def _atom_id(a):
  """Compact atom identity, e.g. ``"HOH A 863 H2"``."""
  L = a.fetch_labels()
  return (f"{L.resname.strip()} {L.chain_id.strip()} "
          f"{L.resseq.strip()} {L.name.strip()}")


def _water_id(ag):
  """Compact water identity, e.g. ``"HOH A 863"`` (altloc in parentheses)."""
  rg = ag.parent()
  alt = ag.altloc.strip()
  return (f"{ag.resname.strip()} {rg.parent().id.strip()} {rg.resseq.strip()}"
          + (f" ({alt})" if alt else ""))


def _worst_water_clashes(hier):
  """Inter-water H-H contacts within 2.0 A, closest first.

  One ``(distance, id_a, id_b)`` per contact.
  """
  xyz, wid, atoms = _water_h_sites(hier)
  i, j, d = _water_h_contacts(xyz, wid)
  return [(d[k], _atom_id(atoms[i[k]]), _atom_id(atoms[j[k]]))
          for k in flex.sort_permutation(d)]


def _detect_neutron(pdb_in, hier):
  """Classify a model as neutron- or X-ray-like for O-H distance selection.

  Prefers the deposited experiment type (``EXPDTA`` in PDB,
  ``_exptl.method`` in mmCIF), falling back to the presence of D atoms in
  ``hier``. Returns ``(is_neutron, source)``, ``source`` being a short
  human-readable reason.
  """
  exp = pdb_in.get_experiment_type()
  if not exp.is_empty():
    if exp.is_neutron():            # includes joint X-ray/neutron
      return True, f"experiment type {exp!r}"
    if exp.is_xray() or exp.is_electron_microscopy():
      return False, f"experiment type {exp!r}"
  if _has_deuterium(hier):
    return True, "D atoms present (no conclusive experiment metadata)"
  return False, "no neutron signal (assuming X-ray)"

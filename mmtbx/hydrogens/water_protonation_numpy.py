"""H-bond-aware placement of the two hydrogens on bare water oxygens.

A map-free, H-bond-network-aware placer for water hydrogens. For every
bare water oxygen (any common water residue: HOH, DOD, H2O, WAT, OH2, ...)
it places the two H so they (a) point at real H-bond acceptors, (b) avoid
clashes against the whole structure (including H placed on other waters),
and (c) keep off metal cations -- placed purely from geometry, with no
experimental data/map and no monomer library. It enforces O-H = 0.984 A
(neutron) or 0.957 A (X-ray) and H-O-H = 104.5 deg.

The public entry point is :func:`place_water_hydrogens`, which modifies a
hierarchy in place. The :class:`mmtbx.programs.water_protonation.Program`
wraps it with DataManager I/O, experiment-aware defaults, and a clash
report; the command-line dispatcher is ``mmtbx.naiad``.

Needs ``numpy`` and ``scipy`` (KDTree).
"""

from __future__ import absolute_import, division, print_function

import math
import random

import iotbx.pdb
from libtbx import group_args
import numpy as np
from scitbx import matrix
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
# Lone-pair-directed placement (opt-in). When on, H1 aims at an acceptor's
# lone-pair lobe (derived from its bonded-neighbour geometry) rather than
# its nucleus, giving better D-H...A angles. Helper distances/angles:
_WATER_BOND_HEAVY = 1.9       # max heavy-heavy bond distance for lobe geometry (A)
_WATER_HBOND_HA = 1.8         # nominal H...acceptor distance for the lobe target (A)
_WATER_SP2_LOBE_DEG = 60.0    # half-angle of the two sp2 carbonyl/-late lobes
_WATER_CONE_SAMPLES = 36      # angular samples around the O-H1 cone
# Element-aware "clash-free" thresholds. A candidate H must clear every
# heavy atom by _WATER_MIN_CLEARANCE and every hydrogen (and cation) by
# the larger _WATER_MIN_H_CLEARANCE. ~1.5 A still admits a genuine
# H...acceptor contact (the acceptor heavy atom sits ~1.7-1.8 A from the
# H); the larger H...H bound keeps two protons from facing each other
# (donor-donor) -- a uniform raise to 2.0 would instead reject real
# H-bonds, hence the split by element.
_WATER_MIN_CLEARANCE = 1.5
_WATER_MIN_H_CLEARANCE = 2.0
_WATER_CLEARANCE_RADIUS = 3.0  # neighbour search radius for clearance (A)
# Relaxation sweeps after the greedy pass. Each re-places every water
# against the final positions of all the others, relaxing the water-water
# clashes the order-dependent greedy pass leaves in dense clusters. Each
# sweep ~halves the remaining clash count at ~one placement-pass cost.
_WATER_REFINE_SWEEPS = 5
# Early-stop tolerance: keep refining only while a sweep removes at least
# this many close (<2.0 A) H-H contacts; below it, treat as converged.
# 1 = stop at a true plateau; larger = stop sooner on diminishing returns.
_WATER_REFINE_TOL = 1
# Basin-hopping (opt-in): each round restarts from the best state, randomly
# re-orients the still-clashing waters, and relaxes -- to escape the local
# minima refinement settles into. Seeded for reproducibility.
_WATER_BASIN_SEED = 0
_WATER_BASIN_RELAX = 2        # relaxation sweeps after each random kick
# Metal cations carry a positive charge, so a water that coordinates one
# (via its O lone pair) should keep both H in the hemisphere *away* from
# the metal -- pointing H+ at M(n+) is electrostatically unfavourable.
# Cations within _WATER_CATION_RADIUS of the water O are treated as
# repellers (hemisphere constraint), not acceptors.
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
  """True if the model carries any D atom.

  Only neutron and joint X-ray/neutron refinements model deuterium, so this
  is the structure-only signal for which canonical O-H distance applies. It
  is what ``_detect_neutron`` falls back to when the deposited experiment
  type is absent or inconclusive.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy.

  Returns
  -------
  bool
      True if any atom has element D.
  """
  return any(a.element.strip().upper() == "D" for a in hier.atoms())


def count_environment_hydrogens(hier):
  """Number of H/D carried by atoms outside water residues.

  Water orientation is decided against the surrounding atoms, so hydrogens
  on the rest of the model are what the clash test sees, and a bonded H is
  what marks an N a donor rather than an acceptor. This counts them so a
  caller can tell whether that information is available.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy.

  Returns
  -------
  int
      Count of H/D on non-water residues.
  """
  return sum(1 for ag in hier.atom_groups() if not _is_water(ag.resname)
             for a in ag.atoms() if a.element_is_hydrogen())


def _is_water(resname):
  """True for any common water alias (HOH, DOD, H2O, WAT, OH2, ...)."""
  return (iotbx.pdb.common_residue_names_get_class(resname.strip().upper())
          == "common_water")


def _fibonacci_sphere(n):
  """Roughly-uniform unit vectors over the sphere (Fibonacci spiral).

  Used as generic fallback directions when no acceptor gives a clash-free
  H -- dense enough to find an open "away from everything" pocket.

  Parameters
  ----------
  n : int
      Number of points to generate.

  Returns
  -------
  tuple of tuple of float
      ``n`` unit vectors as ``(x, y, z)`` tuples.
  """
  golden = math.pi * (3.0 - math.sqrt(5.0))
  pts = []
  for i in range(n):
    y = 1.0 - 2.0 * (i + 0.5) / n
    r = math.sqrt(max(0.0, 1.0 - y * y))
    th = golden * i
    pts.append((math.cos(th) * r, y, math.sin(th) * r))
  return tuple(pts)


_WATER_FALLBACK_DIRECTIONS = _fibonacci_sphere(64)

# Empty neighbour-slot list, shared so the common "no placed H nearby" case
# allocates nothing.
_EMPTY_SLOTS = np.zeros(0, dtype=np.intp)


def _as_np(cols):
  """Stack scitbx col vectors into an (n, 3) float array."""
  if not cols:
    return np.zeros((0, 3))
  return np.array([c.elems for c in cols], dtype=float)


# The H2 cone is sampled at fixed azimuths, so its cosines and sines are
# constants; likewise the normalized fallback sphere.
_CONE_COS = np.array([math.cos(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                      for k in range(_WATER_CONE_SAMPLES)])
_CONE_SIN = np.array([math.sin(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                      for k in range(_WATER_CONE_SAMPLES)])
_FALLBACK_NP = _as_np([matrix.col(v).normalize()
                       for v in _WATER_FALLBACK_DIRECTIONS])


def _ortho_frame(d1):
  """Build an orthonormal frame around a direction.

  Mirrors ``scitbx.matrix.col.ortho()`` operation for operation.

  Parameters
  ----------
  d1 : numpy.ndarray
      Unit direction (the O-H1 axis), shape (3,).

  Returns
  -------
  tuple of numpy.ndarray
      Two unit vectors ``p, q`` with ``(d1, p, q)`` mutually perpendicular
      -- the frame for sampling the H2 cone around ``d1``.
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
  """A random unit vector (rejection-free Gaussian normalization).

  Parameters
  ----------
  rng : random.Random
      Seeded RNG (for reproducibility).

  Returns
  -------
  numpy.ndarray
      A uniformly random unit vector.
  """
  while True:
    v = np.array((rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)))
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n > 1e-6:
      return v / n


def _strip_water_hydrogens(hier):
  """Remove every H/D from water residues, leaving bare O.

  The stripped waters are re-placed from scratch by the caller.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Hierarchy to modify in place.

  Returns
  -------
  int
      Number of H/D atoms removed.
  """
  n = 0
  for ag in hier.atom_groups():
    if not _is_water(ag.resname):
      continue
    for a in list(ag.atoms()):
      if a.element_is_hydrogen():
        ag.remove_atom(a)
        n += 1
  return n


def _new_h_atom(name, element, xyz, occ, b, hetero=False):
  """Build a hierarchy atom for a placed water H.

  The atom is given a blank segid, like the parent O.

  Parameters
  ----------
  name : str
      Atom name (e.g. ``" H1 "``).
  element : str
      Element symbol (``"H"`` or ``"D"``).
  xyz : tuple of float
      Cartesian coordinates.
  occ : float
      Occupancy (inherited from the parent O).
  b : float
      B factor (inherited from the parent O).
  hetero : bool, optional
      HETATM flag (inherited from the parent O, so the placed H are
      written as the same record type as the water they belong to).

  Returns
  -------
  iotbx.pdb.hierarchy.atom
      The new atom.
  """
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
  """First of ``<element>1`` / ``<element>2`` not already used in a residue.

  A water being completed already carries one proton, which may be named
  either H1 or H2, so the new one takes whichever slot is free.

  Parameters
  ----------
  existing_names : set of str
      Stripped atom names already present in the residue.
  element : str
      Element symbol (``"H"`` or ``"D"``).

  Returns
  -------
  str or None
      A 4-character PDB atom name, or None if both slots are taken.
  """
  for di in (1, 2):
    name = f"{element}{di}"
    if name not in existing_names:
      return f" {name} "
  return None


def _sp2_plane_normal(atoms, static_tree, c, exclude):
  """Unit normal of the sp2 plane around an atom.

  Computed from the heavy substituents of atom index ``c`` (the atom a
  terminal acceptor O is bonded to) other than ``exclude``.

  Parameters
  ----------
  atoms : list of iotbx.pdb.hierarchy.atom
      All atoms (positionally indexed).
  static_tree : scipy.spatial.KDTree
      KDTree over the atom coordinates.
  c : int
      Index of the atom whose substituent plane is wanted.
  exclude : int
      Index of the bonded acceptor O to exclude (supplies one in-plane
      vector).

  Returns
  -------
  scitbx.matrix.col or None
      Unit normal of the plane, or None if underdetermined.
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
  """Estimate lone-pair lobe directions for every acceptor atom.

  The lobes are derived from each acceptor's bonded-neighbour geometry
  (library-free):

  - terminal O (1 bond, carbonyl/carboxylate): two in-plane sp2 lobes
    ~``2 * _WATER_SP2_LOBE_DEG`` apart, straddling the direction away from
    the bonded atom, in that atom's substituent plane.
  - otherwise: a single lobe along the open hemisphere (opposite the sum
    of the bond directions).

  Acceptors whose geometry can't be resolved get an empty list, and the
  caller falls back to aiming at the nucleus.

  Parameters
  ----------
  atoms : list of iotbx.pdb.hierarchy.atom
      All atoms (positionally indexed).
  static_tree : scipy.spatial.KDTree
      KDTree over the atom coordinates.
  donor_n : set of int
      Indices of N atoms that carry an H (donors, excluded as acceptors).

  Returns
  -------
  dict
      Maps each acceptor-atom index to a list of lone-pair lobe unit
      vectors (``scitbx.matrix.col``); the list is empty when the geometry
      is underdetermined.
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
      # An atom on top of this one has no bond direction to contribute.
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


class _WaterHydrogenPlacer(object):
  """The H-bond-aware water-hydrogen placer (engine behind
  :func:`place_water_hydrogens`).

  Holds the shared placement state -- the static-atom KDTree, the
  donor/acceptor bookkeeping, the geometry constants, the per-water
  neighbour blocks and the placed protons -- and the per-water placement,
  refinement and basin-hopping logic. :func:`place_water_hydrogens` is a
  thin wrapper that constructs this and calls :meth:`run`.

  Everything the placement scores against is static except the placed water
  H themselves: the model never moves and neither do the water O, and every
  candidate H lies exactly ``oh_length`` from its O. Each water therefore
  gets one fixed block of static neighbours and one fixed list of
  neighbouring waters, built once in :meth:`run` and reused by the greedy
  pass, every relaxation sweep and every basin round.

  The constructor only records the parameters; all the work (and all the
  derived state) happens in :meth:`run`, which modifies the hierarchy in
  place. See :func:`place_water_hydrogens` for the parameter semantics.
  """

  def __init__(self, hier, oh_length=None, element=None,
               n_refine=_WATER_REFINE_SWEEPS, refine_tol=_WATER_REFINE_TOL,
               n_basin=0, existing_h="keep", lone_pair_directed=False,
               on_state=None):
    self.hier = hier
    self.oh_length = oh_length
    self.element = element
    self.n_refine = n_refine
    self.refine_tol = refine_tol
    self.n_basin = n_basin
    self.existing_h = existing_h
    self.lone_pair_directed = lone_pair_directed
    self.on_state = on_state

    # Placed-H coordinates, one slot per placed proton, set up in run():
    # placed_np is the same data as an (n, 3) array.
    self.placed_coords = []
    self.placed_np = None
    # Active-set bookkeeping for the relaxation sweeps (see _dirty), sized
    # in run(): the tick at which each water was last placed and the tick
    # at which its H last moved, off one monotonic counter.
    self.w_placed_at = []
    self.w_moved_at = []
    self.tick = 0
    # (water index, [(atom, slot, di), ...], fixed_d1)
    self.records = []
    # (residue id, (metal element, distance) or None, action)
    self.partial_waters = []

  def _static_block(self, wi, cands):
    """Candidates against water ``wi``'s static neighbours.

    Parameters
    ----------
    wi : int
        Index of the water being placed.
    cands : numpy.ndarray
        Candidate H positions, shape (k, 3).

    Returns
    -------
    tuple
        ``(d, ok)``: the (k, m) squared-distance block and whether each
        candidate clears every one of those neighbours. ``d`` is None when
        the water has no static neighbour at all.
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

    The static block of a candidate set that never moves is the same on
    every pass, so it is computed once and kept as the two things the test
    reads out of it: the per-candidate flag and the per-candidate minimum.
    Both are exact reductions of the block, so the cached form is the block,
    only smaller (two values per candidate instead of one per neighbour).

    Parameters
    ----------
    wi : int
        Index of the water being placed.
    cands : numpy.ndarray
        Candidate H positions, shape (k, 3); an invariant set (the acceptor
        directions or the fallback sphere).
    key : str
        Slot of the water's geometry dict to cache under.

    Returns
    -------
    tuple
        ``(ok, mins)``: the per-candidate flag and the per-candidate
        minimum squared distance, or None for the latter when the water has
        no static neighbour and there is nothing to reduce.
    """
    g = self.w_geom[wi]
    st = g[key]
    if st is None:
      d, ok = self._static_block(wi, cands)
      st = g[key] = (ok, None if d is None else d.min(axis=1))
    return st

  def _clear(self, wi, cands, nbr_slots, static_key=None):
    """Clearance test over the candidate H positions of one water.

    Every candidate lies exactly ``oh_length`` from the water O, so the
    cached static block of water ``wi`` (a ball of ``oh_length +
    _WATER_CLEARANCE_RADIUS`` about the O, own atoms already dropped)
    covers every candidate's own clearance ball; atoms in the block but
    outside a given candidate's ball are harmless, since they can neither
    beat the 3.0 A cap nor trip a threshold. The whole test is then one
    dense (candidate x neighbour) block instead of a query and a Python
    loop per candidate.

    Parameters
    ----------
    wi : int
        Index of the water being placed.
    cands : numpy.ndarray
        Candidate H positions, shape (k, 3).
    nbr_slots : numpy.ndarray
        Placed-H slots that may lie near this water (its own excluded).
    static_key : str or None, optional
        Geometry-dict slot holding the cached static half of the test, for
        a candidate set that does not move (the acceptor directions, the
        fallback sphere). None for the cone, which is rebuilt around each
        pass's O-H1 axis and so has nothing to reuse.

    Returns
    -------
    tuple of numpy.ndarray
        ``(min_dist, ok)``, one entry per candidate: the distance to the
        nearest non-own atom within ``_WATER_CLEARANCE_RADIUS`` (the search
        radius if nothing is near, a soft tie-breaker) and whether the
        candidate clears every heavy atom by ``_WATER_MIN_CLEARANCE`` and
        every hydrogen (static or placed) by ``_WATER_MIN_H_CLEARANCE``.
    """
    # Squared distances throughout: the thresholds and the cap are exact
    # squares and sqrt is monotone, so both the comparisons and the
    # per-candidate minimum are unchanged.
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
      Q = self.placed_np[nbr_slots]
      dx = cands[:, 0, None] - Q[None, :, 0]
      dy = cands[:, 1, None] - Q[None, :, 1]
      dz = cands[:, 2, None] - Q[None, :, 2]
      d = dx * dx + dy * dy + dz * dz
      np.minimum(best, d.min(axis=1), out=best)
      ok &= ~(d < _WATER_MIN_H_CLEARANCE ** 2).any(axis=1)
    return np.sqrt(best), ok

  @staticmethod
  def _cat_ok(cat, pts, o):
    """Whether each candidate is clear of every nearby cation.

    True where the candidate is in the hemisphere away from every nearby
    cation, i.e. ``(H - O).(cation - O) <= 0`` for all of them.

    Parameters
    ----------
    cat : numpy.ndarray
        O->cation unit directions, shape (m, 3).
    pts : numpy.ndarray
        Candidate H positions, shape (k, 3).
    o : numpy.ndarray
        Water O coordinates, shape (3,).

    Returns
    -------
    numpy.ndarray
        Boolean, one entry per candidate.
    """
    if not len(cat):
      return np.ones(len(pts), dtype=bool)
    dx = pts[:, 0, None] - o[0]
    dy = pts[:, 1, None] - o[1]
    dz = pts[:, 2, None] - o[2]
    dp = dx * cat[None, :, 0] + dy * cat[None, :, 1] + dz * cat[None, :, 2]
    return (dp <= 0.0).all(axis=1)

  def _stats(self):
    """:func:`_water_clash_stats` from the engine's own arrays.

    The placer already knows every H it placed and which water each belongs
    to, and the water H it did not place never move, so the whole-hierarchy
    walk the public helper does is unnecessary between sweeps.

    Returns
    -------
    tuple
        The ``_water_clash_stats`` tuple.
    """
    ns = len(self.placed_coords)
    if ns:
      X = np.concatenate((self.placed_np[:ns], self.wh_xyz))
      W = np.concatenate((self.slot_wid[:ns], self.wh_wid))
    else:
      X, W = self.wh_xyz, self.wh_wid
    n = len(X)
    if n < 2:
      return n, 0, 0, 0, None
    pairs = KDTree(X).query_pairs(2.0, output_type="ndarray")
    if len(pairs):
      pairs = pairs[W[pairs[:, 0]] != W[pairs[:, 1]]]  # drop same-water H
    if not len(pairs):
      return n, 0, 0, 0, None
    D = X[pairs[:, 0]] - X[pairs[:, 1]]
    d = np.sqrt(D[:, 0] * D[:, 0] + D[:, 1] * D[:, 1] + D[:, 2] * D[:, 2])
    return (n, int(len(d)), int((d < 1.8).sum()), int((d < 1.5).sum()),
            float(d.min()))

  def _nearest_cation(self, o_xyz, own_idx):
    """Closest metal cation coordinating a water O, if any.

    Uses ``_WATER_METAL_COORD_RADIUS`` (a first-shell bond), not the looser
    ``_WATER_CATION_RADIUS`` used for proton repulsion. Reporting only; it
    does not affect placement.

    Parameters
    ----------
    o_xyz : tuple of float
        Coordinates of the water O.
    own_idx : set of int
        Static-atom indices belonging to this water (excluded).

    Returns
    -------
    tuple or None
        ``(element, distance)`` of the closest coordinating cation, or None.
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

    The acceptor list, the O->acceptor directions and the candidate points
    built from them, and the O->cation directions, all depend only on the
    static atoms and the water O -- none of which move. They are therefore
    built on first use and reused by the greedy pass, every relaxation
    sweep and every basin round, instead of being re-derived per pass.

    Parameters
    ----------
    wi : int
        Index of the water.

    Returns
    -------
    dict
        ``o`` (O coordinates), ``acceptors`` (atom indices, nearest first),
        ``acc_dirs``/``acc_pts`` (O-H unit directions and H positions for
        them), ``cat`` (O->cation unit directions), ``sph`` (the lazily
        built fallback-sphere H positions) and ``acc_static``/``sph_static``
        (the cached static half of the clearance test for those two
        candidate sets, see :meth:`_static_clear`).
    """
    g = self.w_geom[wi]
    if g is not None:
      return g
    atoms = self.atoms
    o_xyz = self.w_o_col[wi]
    own_idx = self.w_own[wi]
    # Acceptors (static O/N/F/S/Cl) within range, nearest first. Keep them
    # all so H1 can skip a near acceptor whose own H is in the way.
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

      Aimed at a lone-pair lobe when lone-pair-directed placement is on and
      a lobe faces the water, else straight at the acceptor nucleus.
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

    # Nearby metal cations: O->cation unit directions. A placed H should
    # stay out of the metal's hemisphere (see _cat_ok).
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
    """Compute the two H positions for one water O.

    Clash-aware against the current environment (the water's own static
    atoms and own placed H are excluded).

    Parameters
    ----------
    wi : int
        Index of the water being placed.
    nbr_slots : numpy.ndarray
        Placed-H slots that may lie near this water (its own excluded).
    fixed_d1 : numpy.ndarray or None, optional
        Unit O-H1 direction to hold fixed instead of searching for one,
        used when completing a water that already carries one proton. The
        returned ``h1_xyz`` then just restates that direction at the
        canonical bond length and the caller ignores it; only ``h2_xyz``
        is new.

    Returns
    -------
    tuple of numpy.ndarray
        ``(h1_xyz, h2_xyz)``, the placed H1 and H2 positions.
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
    # clash-free; else the best direction over acceptors + a dense fallback
    # sphere, ranked (away-from-cation, clash-free, clearance). ``h1_k`` is
    # the position in the acceptor list H1 took, or -1.
    if fixed_d1 is not None:
      # H1 is the deposited proton; find the acceptor it donates to so H2
      # aims elsewhere.
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

    # Orthonormal frame for the H2 cone (d1, p, q mutually perpendicular)
    p, q = _ortho_frame(d1)

    # H2: rank the cone angles by (away-from-cation, clash-free, best
    # alignment over the acceptors H1 did not take, clearance).
    cone_dirs = (self.cos_hoh * d1
                 + self.sin_hoh * (_CONE_COS[:, None] * p
                                   + _CONE_SIN[:, None] * q))
    cone_pts = o + self.oh_length * cone_dirs
    c_best, c_ok = self._clear(wi, cone_pts, nbr_slots)
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

  def _store(self, slots, h1, h2):
    """Write one water's new H positions to the model and the arrays.

    Returns True if any of them is not where it already was -- what the
    sweep's active set is stamped from (see :meth:`_dirty`).
    """
    moved = False
    for atom, slot, di in slots:
      xyz = tuple((h1 if di == 1 else h2).tolist())
      if xyz != self.placed_coords[slot]:
        moved = True
      self.placed_coords[slot] = xyz
      self.placed_np[slot] = xyz
      atom.set_xyz(xyz)
    return moved

  def _dirty(self, wi):
    """Whether re-placing water ``wi`` could move its H.

    A water is placed against its static surroundings, which never move,
    and the placed H of its neighbouring waters -- ``w_wnbr[wi]``, the same
    set ``w_nbr_slots[wi]`` is built from. If none of those H has moved
    since this water was last placed, the placement re-derives exactly the
    two positions it already holds, so the sweep can skip it. Waters are
    stamped with a monotonic tick, bumped once per water per sweep, so a
    neighbour that moves earlier in the *same* sweep still counts.

    Parameters
    ----------
    wi : int
        Index of the water.

    Returns
    -------
    bool
        True if the water must be re-placed.
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
    return False

  def _apply_sweep(self):
    """Run one relaxation sweep.

    Re-places every recorded water whose neighbourhood has changed against
    the current positions of all the others (its own H excluded), updating
    ``placed_coords`` and the atoms in place. A water placed later in the
    sweep therefore sees the H the earlier ones just moved, at the position
    they now hold. The waters whose neighbours have not moved since they
    were last placed would re-derive the positions they already hold, so
    they are skipped: the result is the same as sweeping all of them, and
    later sweeps touch only the shrinking active set (see :meth:`_dirty`).
    """
    for wi, slots, fixed_d1 in self.records:
      self.tick += 1
      if self._dirty(wi):
        if self._store(slots,
                       *self._place_one(wi, self.w_nbr_slots[wi], fixed_d1)):
          self.w_moved_at[wi] = self.tick
      self.w_placed_at[wi] = self.tick

  def _snapshot(self):
    """Capture the placed-H state: coordinates and both stamp arrays.

    Whether a water may be skipped is a function of the coordinates and the
    stamps together, so the two are captured and restored as one.

    Returns
    -------
    tuple
        ``(coords, placed_at, moved_at)``.
    """
    return (list(self.placed_coords), list(self.w_placed_at),
            list(self.w_moved_at))

  def _restore(self, snap):
    """Reset all placed H to a snapshot from :meth:`_snapshot`.

    Parameters
    ----------
    snap : tuple
        ``(coords, placed_at, moved_at)``.
    """
    coords, placed_at, moved_at = snap
    for _wi, slots, _fixed in self.records:
      for atom, slot, di in slots:
        self.placed_coords[slot] = coords[slot]
        self.placed_np[slot] = coords[slot]
        atom.set_xyz(coords[slot])
    self.w_placed_at = list(placed_at)
    self.w_moved_at = list(moved_at)

  def _clashing_records(self):
    """Find records whose placed H still clash.

    Returns
    -------
    list of int
        Indices into ``records`` of waters with a placed H that fails the
        clash gate.
    """
    bad = []
    for ri, (wi, slots, _fixed) in enumerate(self.records):
      pts = self.placed_np[[slot for _, slot, _ in slots]]
      if not self._clear(wi, pts, self.w_nbr_slots[wi])[1].all():
        bad.append(ri)
    return bad

  def _kick(self, ri, rng):
    """Re-orient one water to a random orientation.

    H1 is placed along a random axis and H2 on the H-O-H cone at a random
    azimuth. Used to kick a water out of its local minimum during
    basin-hopping. A water being completed carries a fixed O-H1 (its
    deposited proton), so only the cone azimuth is randomized there.

    Parameters
    ----------
    ri : int
        Index into ``records`` of the water to kick.
    rng : random.Random
        Seeded RNG (for reproducibility).
    """
    wi, slots, fixed_d1 = self.records[ri]
    o = self._water_geom(wi)["o"]
    d1 = fixed_d1 if fixed_d1 is not None else _rand_unit(rng)
    p, q = _ortho_frame(d1)
    theta = rng.uniform(0.0, 2.0 * math.pi)
    d2 = self.cos_hoh * d1 + self.sin_hoh * (math.cos(theta) * p
                                             + math.sin(theta) * q)
    self.tick += 1
    if self._store(slots, o + self.oh_length * d1, o + self.oh_length * d2):
      self.w_moved_at[wi] = self.tick

  def run(self):
    """Place H on every bare water, refine, optionally basin-hop, keep best.

    Modifies the hierarchy in place. Returns the kept-state label (see
    :func:`place_water_hydrogens`).
    """
    hier = self.hier
    # Resolve before stripping, which would remove the D this keys on.
    if self.oh_length is None:
      self.oh_length = (_WATER_OH_NEUTRON if _has_deuterium(hier)
                        else _WATER_OH_XRAY)
    # The single-H set is only needed from *before* the strip, so the extra
    # hierarchy pass is done only in the mode that strips; otherwise the
    # walk below sees the same protons and counts them itself.
    single_h = None
    if self.existing_h == "reorient":
      single_h = {ag.memory_id() for ag in hier.atom_groups()
                  if _is_water(ag.resname)
                  and sum(1 for a in ag.atoms()
                          if a.element_is_hydrogen()) == 1}
      _strip_water_hydrogens(hier)

    sel = hier.atoms()
    atoms = list(sel)
    if not atoms:
      return None
    self.atoms = atoms

    # Static neighbours (protein, ligands, water O, any pre-existing H) are
    # fixed, so the KDTree over them is built once. The water H we place are
    # tracked separately, by slot, in placed_coords/placed_np.
    self.static_np = sel.extract_xyz().as_numpy_array()
    self.static_tree = KDTree(self.static_np)
    self.static_is_h = np.array([a.element_is_hydrogen() for a in atoms],
                                dtype=bool)
    self.static_thr = np.where(self.static_is_h, _WATER_MIN_H_CLEARANCE ** 2,
                               _WATER_MIN_CLEARANCE ** 2)
    # Key by memory_id() (stable C++ identity), not id(): iotbx returns a
    # fresh Python wrapper on each atom access, so id() is not stable across
    # the hier.atoms() and ag.atoms() calls used to build the per-water
    # own_idx below.
    pos_by_id = {a.memory_id(): i for i, a in enumerate(atoms)}

    # N atoms that already carry an H are donors, not acceptors -- their lone
    # pair is occupied/conjugated (amide, ammonium, guanidinium, protonated
    # His ring N, ...). An N's bonded H thus encodes its donor/acceptor role,
    # charge and tautomer, giving a name-free, library-free donor test. O
    # always accepts (it keeps lone pairs even when donating), so only N is
    # filtered. Inert on an unprotonated model (no H -> no N flagged -> every
    # N stays an acceptor).
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

    # Gather the waters to protonate. One walk of each water residue serves
    # the O lookup, the existing-H list, the free-name set, the own-atom
    # index set and the clash-stats bookkeeping.
    waters = []
    wh_xyz = []
    wh_wid = []
    for wgid, ag in enumerate(g for g in hier.atom_groups()
                              if _is_water(g.resname)):
      o = None
      existing = []
      names = set()
      own_idx = set()
      for a in ag.atoms():
        names.add(a.name.strip())
        i = pos_by_id.get(a.memory_id())
        if i is not None:
          own_idx.add(i)
        if a.element_is_hydrogen():
          existing.append(a)
          wh_xyz.append(a.xyz)
          wh_wid.append(wgid)
        elif o is None and a.element.strip().upper() == "O":
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

    # Per-water constants. Neither the static atoms nor the water O move, so
    # every neighbour list the placement needs is built once, here, and
    # reused by the greedy pass, every relaxation sweep and every basin
    # round.
    n = len(waters)
    self.w_geom = [None] * n
    self.placed_coords = []
    self.placed_np = np.zeros((2 * n, 3))
    self.slot_wid = np.zeros(2 * n, dtype=np.int64)
    self.records = []   # (water index, [(atom, slot, di), ...], fixed_d1)
    self.w_wnbr = []
    # Every water is dirty for the first sweep: the greedy pass placed each
    # one against the H of the earlier waters only.
    self.w_placed_at = [-1] * n
    self.w_moved_at = [0] * n
    self.tick = 0
    if n:
      o_pts = np.array([w[1].xyz for w in waters], dtype=float)
      # Order most-crowded first: a water boxed in by many neighbours is the
      # hardest to satisfy, so placing it while few water H are fixed gives
      # it the most freedom; roomy waters have options left over and adapt
      # around it. ``crowd`` is the neighbour count within the clash radius
      # (excluding the water's own atoms).
      crowd = [sum(1 for j in nb if j not in waters[k][2]) for k, nb in
               enumerate(self.static_tree.query_ball_point(
                 o_pts, _WATER_CLEARANCE_RADIUS,
                 return_sorted=False))]
      order = sorted(range(n), key=lambda k: crowd[k], reverse=True)
      waters = [waters[k] for k in order]
      o_pts = o_pts[order]
      self.w_own = [w[2] for w in waters]
      self.w_o_col = [matrix.col(w[1].xyz) for w in waters]
      # query_ball_point leaves a single point's indices unsorted but sorts
      # a batch's; the acceptor order is the tie-break of the nearest-first
      # sort, so the batch must not sort either.
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
      # O lie within clearance + 2 oh_length can ever hold one near this
      # water's candidates. A neighbour map over the (fixed) water O
      # therefore replaces the growing KDTree over the placed H outright.
      self.w_wnbr = [[j for j in nb if j != wi] for wi, nb in
                     enumerate(KDTree(o_pts).query_ball_point(
                       o_pts,
                       _WATER_CLEARANCE_RADIUS + 2.0 * self.oh_length + 0.01,
                       return_sorted=False))]

    # Initial greedy pass over the ordered waters, each avoiding the H
    # already placed on earlier ones. ``records`` keeps per-water (index,
    # placed-H slots, fixed_d1) so the refinement sweeps can re-place each
    # water against the *final* environment, breaking the order-dependence
    # of the greedy pass.
    slots_of = [()] * n

    def nbr_slots(wi):
      """Placed-H slots of every water neighbouring water ``wi``."""
      nbr = [s for wj in self.w_wnbr[wi] for s in slots_of[wj]]
      return np.array(nbr, dtype=np.intp) if nbr else _EMPTY_SLOTS

    for wi in range(n):
      ag, o, own_idx, fixed_d1, existing, existing_names, wgid = waters[wi]
      if self.element is not None:
        proton_element = self.element
      elif existing:
        proton_element = existing[0].element.strip().upper()
      else:
        proton_element = "D" if ag.resname.strip().upper() == "DOD" else "H"
      h1, h2 = self._place_one(wi, nbr_slots(wi), fixed_d1)

      # Completing builds only the cone proton; names take whichever of
      # H1/H2 is free.
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
        self.placed_np[slot] = xyz
        self.slot_wid[slot] = wgid
      if slots:
        self.records.append((wi, slots, fixed_d1))
        slots_of[wi] = tuple(s for _, s, _ in slots)

    # The placed-H neighbourhood of every water, now that all of them exist.
    self.w_nbr_slots = [nbr_slots(wi) for wi in range(n)]

    # Refinement sweeps: re-place each water against the *final* set of all
    # placed H (its own excluded), so a water no longer ignores neighbours
    # placed after it. This relaxes the water-water clashes the greedy pass
    # leaves in dense clusters. With refine_tol > 0 it stops early once a
    # sweep reduces the close (<2.0 A) contact count by fewer than refine_tol
    # (n_refine is then a cap); refine_tol == 0 runs all n_refine sweeps.
    # Refinement isn't strictly monotonic, so the best state is kept and
    # restored at the end -- the result is never worse than any state seen.
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
        break  # gain below tolerance -- converged
      prev_n20 = stats[1]

    # Basin-hopping (optional): each round restarts from the best state, kicks
    # the still-clashing waters to a random orientation, relaxes, and keeps
    # the result if it improved. Deterministic (seeded). Helps only where a
    # better orientation exists -- it re-aims the H, so a clash rooted in a
    # misplaced water O won't resolve.
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
                          on_state=None):
  """Place the two H on every bare water, H-bond-aware.

  For each water residue missing H (any common water alias -- HOH, DOD,
  H2O, WAT, OH2, ...):

  - H1 along O -> nearest H-bond acceptor (within
    ``_WATER_ACCEPTOR_RADIUS``, element in ``_WATER_ACCEPTOR_ELEMENTS``,
    excluding N that already carry an H -- those are donors, not
    acceptors) that yields a *clash-free* H -- acceptors whose own H would
    collide with the placed proton are skipped; falls back to the
    max-clearance direction (over a dense sphere) if none is clash-free.
  - H2 on the cone of half-angle ``_WATER_HOH_DEG`` around the O-H1 axis.
    Among clash-free cone angles it points toward a second acceptor when
    one exists, else maximizes clearance; if no angle is clash-free it
    falls back to the clearest one.

  Clash avoidance is global: every candidate position is scored against
  all other atoms *and all other waters' placed H*, so neighbouring waters
  don't aim H at each other. Positions stay at least
  ``_WATER_MIN_CLEARANCE`` from everything else where the local packing
  allows it; in genuinely crowded pockets the clearest available position
  is used.

  Waters are placed most-crowded first (the hardest to satisfy get the
  most freedom), then up to ``n_refine`` refinement sweeps re-place each
  water against the final environment to relax water-water clashes; the
  best sweep is always kept, so the result is never worse than any state
  seen. New H inherit the parent O's occupancy and B factor.

  Modifies *hier* in place. By default it is idempotent: water residues
  already carrying H are left untouched. A water carrying exactly one H is
  left alone too, since omitting the second proton is a common way of
  writing hydroxide; ``existing_h="complete"`` builds the missing partner
  on the cone of the deposited proton instead.

  ``oh_length`` defaults to picking the canonical distance from the model:
  neutron if it carries any D atom, X-ray otherwise. ``mmtbx.naiad`` reads
  the deposited experiment type instead, which is stronger evidence, and
  passes the result explicitly.

  This is a thin wrapper around :class:`_WaterHydrogenPlacer`; the
  parameters and return value are documented here.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy; modified in place.
  oh_length : float or None, optional
      O-H bond length in A, positive. None (default) picks
      ``_WATER_OH_NEUTRON`` (0.984) if the model contains D, else
      ``_WATER_OH_XRAY`` (0.957).
  element : str or None, optional
      Element of the placed atoms, ``"H"`` or ``"D"``, forced on every
      water; None (default) picks per residue -- ``"D"`` for DOD, ``"H"``
      for HOH. Atom names follow the element (``H1``/``H2`` or
      ``D1``/``D2``).
  n_refine : int, optional
      Maximum relaxation sweeps after the greedy pass; 0 disables
      refinement.
  refine_tol : int, optional
      Stop refining once a sweep removes fewer than this many close
      (<2.0 A) H-H contacts (so ``n_refine`` is a cap). 1 = stop at a true
      plateau, larger = stop sooner on diminishing returns, 0 = run all
      ``n_refine`` sweeps.
  n_basin : int, optional
      Basin-hopping rounds after refinement (default 0 = off). Each
      restarts from the best state, randomly re-orients the waters still in
      a clash, relaxes, and keeps the result if it improved (deterministic,
      seeded). Helps only where a better orientation exists.
  existing_h : str, optional
      What to do with a water that already carries H: ``"keep"``
      (default, leave it untouched), ``"complete"`` (a water with exactly
      one H gets its partner built on that proton's cone, inheriting its
      element), or ``"reorient"`` (strip all water H and re-place both).
  lone_pair_directed : bool, optional
      If True, each O-H aims at an acceptor's lone-pair lobe (estimated
      from the acceptor's bonded-neighbour geometry) rather than its
      nucleus, for better D-H...A angles.
  on_state : callable, optional
      If given, called as ``on_state(label, stats)`` once the placement
      reaches each state -- ``"initial"`` after the greedy pass, then
      ``"sweep N"`` after each refinement sweep (and ``"basin N.M"`` during
      basin-hopping) -- where ``stats`` is the ``_water_clash_stats``
      tuple. Lets a caller stream a progress table row by row.

  Returns
  -------
  libtbx.group_args
      ``kept_label``: the label of the kept state (``"initial"``,
      ``"sweep N"`` or ``"basin N.M"``) when refinement or basin-hopping
      ran, else None. ``partial_waters``: one
      ``(residue_id, metal, action)`` per water that carried exactly one H
      on input, where ``metal`` is ``(element, distance)`` for a
      coordinating cation or None, and ``action`` is ``"kept"``,
      ``"completed"`` or ``"stripped"`` after ``existing_h``.
  """
  placer = _WaterHydrogenPlacer(
    hier, oh_length=oh_length, element=element, n_refine=n_refine,
    refine_tol=refine_tol, n_basin=n_basin,
    existing_h=existing_h,
    lone_pair_directed=lone_pair_directed,
    on_state=on_state)
  kept_label = placer.run()
  return group_args(kept_label=kept_label,
                    partial_waters=placer.partial_waters)


def _water_clash_stats(hier):
  """Count water-H vs water-H contacts between different waters.

  Heavy-atom contacts are not counted -- the placer's clash gate keeps
  those clear by construction.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy.

  Returns
  -------
  tuple
      ``(n_placed, n_lt_20, n_lt_18, n_lt_15, closest)``: the number of
      placed water H, the counts of inter-water H-H contacts below
      2.0/1.8/1.5 A, and the closest such distance (None if no pair is
      within 2.0 A).
  """
  wh = []
  for ag in hier.atom_groups():
    if not _is_water(ag.resname):
      continue
    wid = ag.memory_id()
    for a in ag.atoms():
      if a.element_is_hydrogen():
        wh.append((tuple(a.xyz), wid))
  if len(wh) < 2:
    return len(wh), 0, 0, 0, None
  tree = KDTree([x for x, _ in wh])
  n20 = n18 = n15 = 0
  worst = None
  for i, (x, wid) in enumerate(wh):
    xc = matrix.col(x)
    for j in tree.query_ball_point(x, 2.0):
      if j <= i or wh[j][1] == wid:  # skip self-pair and same-water H
        continue
      d = (matrix.col(wh[j][0]) - xc).length()
      if worst is None or d < worst:
        worst = d
      n20 += 1
      n18 += d < 1.8
      n15 += d < 1.5
  return len(wh), n20, n18, n15, worst


def _clash_row(label, stats, log):
  """Print one row of the per-sweep clash table.

  Parameters
  ----------
  label : str
      Row label (state name, e.g. ``"sweep 2"``).
  stats : tuple
      A ``_water_clash_stats`` tuple.
  log : file object
      Destination stream.
  """
  _, n20, n18, n15, worst = stats
  w = f"{worst:.2f}" if worst is not None else ">2.0"
  print(f"  {label:<9} <2.0={n20:<5} <1.8={n18:<5} <1.5={n15:<5} closest={w}",
        file=log)


def _atom_id(a):
  """Compact identity string for an atom, e.g. ``"HOH A 863 H2"``.

  Parameters
  ----------
  a : iotbx.pdb.hierarchy.atom
      Atom to identify.

  Returns
  -------
  str
      ``"<resname> <chain> <resseq> <name>"``.
  """
  rg = a.parent().parent()
  return (f"{a.parent().resname.strip()} {rg.parent().id.strip()} "
          f"{rg.resseq.strip()} {a.name.strip()}")


def _water_id(ag):
  """Compact residue identity for a water, e.g. ``"HOH A 863"``.

  Parameters
  ----------
  ag : iotbx.pdb.hierarchy.atom_group
      The water atom group.

  Returns
  -------
  str
      ``"<resname> <chain> <resseq>"``, with the altloc appended in
      parentheses when the group has one.
  """
  rg = ag.parent()
  alt = ag.altloc.strip()
  return (f"{ag.resname.strip()} {rg.parent().id.strip()} {rg.resseq.strip()}"
          + (f" ({alt})" if alt else ""))


def _worst_water_clashes(hier):
  """Closest inter-water H-H contacts -- the offenders behind the counts.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy.

  Returns
  -------
  list of tuple
      ``(distance, id_a, id_b)`` tuples for every water-H vs water-H contact
      between different waters within 2.0 A, closest first.
  """
  wh = []
  for ag in hier.atom_groups():
    if not _is_water(ag.resname):
      continue
    wid = ag.memory_id()
    for a in ag.atoms():
      if a.element_is_hydrogen():
        wh.append((tuple(a.xyz), wid, a))
  if len(wh) < 2:
    return []
  tree = KDTree([x for x, _, _ in wh])
  pairs = []
  for i, (x, wid, a) in enumerate(wh):
    xc = matrix.col(x)
    for j in tree.query_ball_point(x, 2.0):
      if j <= i or wh[j][1] == wid:
        continue
      pairs.append(((matrix.col(wh[j][0]) - xc).length(), a, wh[j][2]))
  pairs.sort(key=lambda t: t[0])
  return [(d, _atom_id(a), _atom_id(b)) for d, a, b in pairs]


def _detect_neutron(pdb_in, hier):
  """Classify a model as neutron- or X-ray-like for O-H distance selection.

  Prefers the deposited experiment type (``EXPDTA`` in PDB, ``_exptl.method``
  in mmCIF); when that is absent or inconclusive, falls back to the presence
  of D atoms (only neutron / joint refinements model deuterium).

  Parameters
  ----------
  pdb_in : iotbx.pdb.input or iotbx.pdb.mmcif.cif_input
      The parsed input (carries the experiment-type metadata).
  hier : iotbx.pdb.hierarchy.root
      The model hierarchy (for the D-atom fallback).

  Returns
  -------
  tuple
      ``(is_neutron, source)``: whether neutron O-H distances should be
      used, and a short human-readable string explaining the decision.
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

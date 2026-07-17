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

Needs ``scipy`` (KDTree).
"""

from __future__ import absolute_import, division, print_function

import math
import random

import iotbx.pdb
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
_WATER_CONE_SAMPLES = 18      # angular samples around the O-H1 cone
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
# The KDTree over already-placed water H is expensive to rebuild from
# scratch every water (O(n) each -> O(n^2) overall). Instead the bulk of
# the placed H live in a "prefix" tree rebuilt only once this many new H
# have accumulated; the few placed since then live in a small "pending"
# tree rebuilt each water. Their union is still every earlier water's H.
_WATER_TREE_REBUILD = 64
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
  "LI", "NA", "K", "RB", "CS",
  "MG", "CA", "SR", "BA",
  "MN", "FE", "CO", "NI", "CU", "ZN", "CD", "HG",
  "AL",
})
_WATER_CATION_RADIUS = 3.0


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
      if a.element.strip().upper() in ("H", "D"):
        ag.remove_atom(a)
        n += 1
  return n


def _new_h_atom(name, element, xyz, occ, b):
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
  return atom


def _ortho_frame(d1):
  """Build an orthonormal frame around a direction.

  Parameters
  ----------
  d1 : scitbx.matrix.col
      Unit direction (the O-H1 axis).

  Returns
  -------
  tuple of scitbx.matrix.col
      Two unit vectors ``p, q`` with ``(d1, p, q)`` mutually perpendicular
      -- the frame for sampling the H2 cone around ``d1``.
  """
  helper = matrix.col((1.0, 0.0, 0.0))
  if abs(d1.dot(helper)) > 0.95:
    helper = matrix.col((0.0, 1.0, 0.0))
  p = (helper - d1 * d1.dot(helper)).normalize()
  return p, d1.cross(p)


def _tilt_dir(dir_a, dir_b, hoh_rad):
  """H1 direction near one acceptor but tilted so H2 can reach a second.

  The returned direction is near acceptor ``dir_a`` but tilted *away* from
  acceptor ``dir_b`` in their shared plane, so an H2 on H1's H-O-H cone can
  reach ``dir_b`` -- the balanced two-acceptor placement. When the
  acceptors subtend angle phi, each H ends up ``(104.5 - phi) / 2`` off its
  acceptor, so both H-bonds are near-ideal.

  Parameters
  ----------
  dir_a : scitbx.matrix.col
      Unit direction toward the first acceptor.
  dir_b : scitbx.matrix.col
      Unit direction toward the second acceptor.
  hoh_rad : float
      Target H-O-H angle in radians.

  Returns
  -------
  scitbx.matrix.col
      The tilted H1 unit direction, or ``dir_a`` if the two acceptors are
      collinear.
  """
  w = dir_b - dir_a * dir_a.dot(dir_b)
  if w.length() < 1e-6:
    return dir_a
  w = w.normalize()
  delta = (hoh_rad - dir_a.angle(dir_b)) / 2.0
  return (dir_a * math.cos(delta) - w * math.sin(delta)).normalize()


def _rand_unit(rng):
  """A random unit vector (rejection-free Gaussian normalization).

  Parameters
  ----------
  rng : random.Random
      Seeded RNG (for reproducibility).

  Returns
  -------
  scitbx.matrix.col
      A uniformly random unit vector.
  """
  while True:
    v = matrix.col((rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1)))
    if v.length() > 1e-6:
      return v.normalize()


def _kick_water(record, placed_coords, oh_length, cos_hoh, sin_hoh, rng):
  """Re-orient one water to a random orientation.

  H1 is placed along a random axis and H2 on the H-O-H cone at a random
  azimuth, updating the water's placed-H coordinates and atoms in place.
  Used to kick a water out of its local minimum during basin-hopping.

  Parameters
  ----------
  record : tuple
      Per-water record ``(o_xyz, own_idx, slots)`` (see
      ``place_water_hydrogens``); only ``o_xyz`` and ``slots`` are used.
  placed_coords : list of tuple
      Global list of placed-H coordinates, updated in place.
  oh_length : float
      O-H bond length.
  cos_hoh, sin_hoh : float
      Cosine and sine of the H-O-H angle.
  rng : random.Random
      Seeded RNG.
  """
  o_xyz, _own_idx, slots = record
  d1 = _rand_unit(rng)
  helper = matrix.col((1.0, 0.0, 0.0))
  if abs(d1.dot(helper)) > 0.95:
    helper = matrix.col((0.0, 1.0, 0.0))
  p = (helper - d1 * d1.dot(helper)).normalize()
  q = d1.cross(p)
  theta = rng.uniform(0.0, 2.0 * math.pi)
  d2 = cos_hoh * d1 + sin_hoh * (math.cos(theta) * p + math.sin(theta) * q)
  dirs = {1: d1, 2: d2}
  for atom, slot, di in slots:
    xyz = tuple(o_xyz + oh_length * dirs[di])
    placed_coords[slot] = xyz
    atom.set_xyz(xyz)


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
    if atoms[k].element.strip().upper() in ("H", "D"):
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
      ej = atoms[j].element.strip().upper()
      d = (matrix.col(atoms[j].xyz) - A).length()
      if (d <= _WATER_NH_BOND) if ej in ("H", "D") else (d <= _WATER_BOND_HEAVY):
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
  donor/acceptor bookkeeping, the geometry constants, and the growing set
  of placed protons (a prefix + pending KDTree pair) -- and the per-water
  placement, refinement and basin-hopping logic. :func:`place_water_hydrogens`
  is a thin wrapper that constructs this and calls :meth:`run`.

  The constructor only records the parameters; all the work (and all the
  derived state) happens in :meth:`run`, which modifies the hierarchy in
  place. See :func:`place_water_hydrogens` for the parameter semantics.
  """

  def __init__(self, hier, oh_length=_WATER_OH_NEUTRON, element=None,
               n_refine=_WATER_REFINE_SWEEPS, refine_tol=_WATER_REFINE_TOL,
               n_basin=0, reorient_existing=False, lone_pair_directed=False,
               joint=False, on_state=None):
    self.hier = hier
    self.oh_length = oh_length
    self.element = element
    self.n_refine = n_refine
    self.refine_tol = refine_tol
    self.n_basin = n_basin
    self.reorient_existing = reorient_existing
    self.lone_pair_directed = lone_pair_directed
    self.joint = joint
    self.on_state = on_state

    # Placed-H state (the two KDTrees and their coordinate backing store),
    # set up in run(). placed_tree covers placed_coords[:n_built], the
    # rebuilt "prefix"; pending_tree covers placed_coords[n_built:], the
    # small "pending" set; their union is every placed water H.
    self.placed_coords = []
    self.placed_tree = None
    self.pending_tree = None
    self.n_built = 0
    # Per-water records (o_xyz, own_idx, [(atom, slot, di), ...]); the
    # refinement/basin passes re-place each against the final environment.
    self.records = []

  def _clearances(self, cands, own_idx, own_slots):
    """Batched clearance test over candidate H positions.

    One KDTree query per tree covers all candidates, rather than one query
    per candidate.

    Parameters
    ----------
    cands : list of scitbx.matrix.col
        Candidate H positions to score.
    own_idx : set of int
        Static-atom indices belonging to the water being placed (excluded).
    own_slots : set of int
        Placed-H slot indices belonging to the water being placed
        (excluded).

    Returns
    -------
    list of tuple
        One ``(min_dist, ok)`` per candidate, where ``min_dist`` is the
        distance to the nearest non-own atom within
        ``_WATER_CLEARANCE_RADIUS`` (the search radius if nothing is near,
        a soft tie-breaker) and ``ok`` is True iff the candidate clears
        every heavy atom by ``_WATER_MIN_CLEARANCE`` and every hydrogen
        (static or placed) by the larger ``_WATER_MIN_H_CLEARANCE``. The
        element split lets an H sit at H-bonding distance from an acceptor
        heavy atom while keeping it away from other protons.
    """
    if not cands:
      return []
    pts = [tuple(c) for c in cands]
    best = [_WATER_CLEARANCE_RADIUS] * len(cands)
    ok = [True] * len(cands)
    for ci, nbrs in enumerate(
        self.static_tree.query_ball_point(pts, _WATER_CLEARANCE_RADIUS)):
      c = cands[ci]
      for j in nbrs:
        if j in own_idx:
          continue
        d = (matrix.col(self.static_coords[j]) - c).length()
        if d < best[ci]:
          best[ci] = d
        thr = _WATER_MIN_H_CLEARANCE if self.static_is_h[j] else _WATER_MIN_CLEARANCE
        if d < thr:
          ok[ci] = False
    # placed H live in two trees (prefix + pending); base maps a tree-local
    # index back to its global placed_coords slot. Placed atoms are all H.
    for tree, base in ((self.placed_tree, 0), (self.pending_tree, self.n_built)):
      if tree is None:
        continue
      for ci, nbrs in enumerate(
          tree.query_ball_point(pts, _WATER_CLEARANCE_RADIUS)):
        c = cands[ci]
        for j in nbrs:
          slot = base + j
          if slot in own_slots:
            continue
          d = (matrix.col(self.placed_coords[slot]) - c).length()
          if d < best[ci]:
            best[ci] = d
          if d < _WATER_MIN_H_CLEARANCE:
            ok[ci] = False
    return list(zip(best, ok))

  def _place_one(self, o_xyz, own_idx, own_slots):
    """Compute the two H positions for one water O.

    Clash-aware against the current environment (the water's own static
    atoms and own placed H are excluded).

    Parameters
    ----------
    o_xyz : scitbx.matrix.col
        Coordinates of the water O.
    own_idx : set of int
        Static-atom indices belonging to this water (excluded).
    own_slots : set of int
        Placed-H slot indices belonging to this water (excluded).

    Returns
    -------
    tuple of scitbx.matrix.col
        ``(h1_xyz, h2_xyz)``, the placed H1 and H2 positions.
    """
    atoms = self.atoms
    # Acceptors (static O/N/F/S/Cl) within range, nearest first. Keep them
    # all so H1 can skip a near acceptor whose own H is in the way.
    acceptors = []
    for i in self.static_tree.query_ball_point(tuple(o_xyz),
                                          _WATER_ACCEPTOR_RADIUS):
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

      Parameters
      ----------
      i : int
          Acceptor atom index.

      Returns
      -------
      scitbx.matrix.col
          The O-H unit direction.
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
    # stay out of the metal's hemisphere (see cation_ok below).
    cation_dirs = []
    for i in self.static_tree.query_ball_point(tuple(o_xyz),
                                          _WATER_CATION_RADIUS):
      if i in own_idx:
        continue
      a = atoms[i]
      if a.element.strip().upper() not in _WATER_CATION_ELEMENTS:
        continue
      v = matrix.col(a.xyz) - o_xyz
      if v.length() < 1e-3:
        continue
      cation_dirs.append(v.normalize())

    def cation_ok(cand):
      """Whether an H candidate is clear of every nearby cation.

      True if ``cand`` is in the hemisphere away from every nearby cation,
      i.e. ``(H - O).(cation - O) <= 0`` for all of them.

      Parameters
      ----------
      cand : scitbx.matrix.col
          Candidate H position.

      Returns
      -------
      bool
          True if the candidate avoids every cation hemisphere.
      """
      if not cation_dirs:
        return True
      dh = cand - o_xyz
      return all(dh.dot(cd) <= 0.0 for cd in cation_dirs)

    # H1: nearest acceptor giving a placement that is away from cations and
    # clash-free; else the best direction over acceptors + a dense fallback
    # sphere, ranked (away-from-cation, clash-free, clearance). Candidate
    # clearances are evaluated in batches (one query set covering all of
    # them) rather than one direction at a time.
    acc_dirs = [accept_dir(i) for i in acceptors]
    acc_pts = [o_xyz + self.oh_length * dv for dv in acc_dirs]
    acc_res = self._clearances(acc_pts, own_idx, own_slots) if acc_pts else []

    # Joint H1/H2 search (opt-in): pick the (H1, H2) orientation maximizing
    # the summed acceptor alignment of both protons, among clash-free,
    # cation-OK placements -- balancing the two H-bonds instead of greedily
    # perfecting H1. H1 candidates are the acceptor directions plus, for each
    # ordered acceptor pair, a "tilt" direction off one acceptor away from
    # the other (so H1 can sit between two acceptors); each viable H1 gets an
    # H2 cone. Falls through to the greedy path if no clash-free pair exists.
    if self.joint and acc_dirs:
      hoh_rad = math.radians(_WATER_HOH_DEG)
      tilts = [_tilt_dir(acc_dirs[a], acc_dirs[b], hoh_rad)
               for a in range(len(acc_dirs))
               for b in range(len(acc_dirs)) if a != b]
      d1_dirs = acc_dirs + tilts
      d1_pts = acc_pts + [o_xyz + self.oh_length * dv for dv in tilts]
      d1_res = acc_res + self._clearances(d1_pts[len(acc_pts):], own_idx, own_slots)

      def _align(dh):
        return max((dh.dot(a) for a in acc_dirs), default=0.0)

      viable = [m for m in range(len(d1_dirs))
                if d1_res[m][1] and cation_ok(d1_pts[m])]
      cone = []                       # (m, d2) for each viable H1 m
      for m in viable:
        fp, fq = _ortho_frame(d1_dirs[m])
        for k in range(_WATER_CONE_SAMPLES):
          th = 2.0 * math.pi * k / _WATER_CONE_SAMPLES
          cone.append((m, self.cos_hoh * d1_dirs[m]
                       + self.sin_hoh * (math.cos(th) * fp + math.sin(th) * fq)))
      cone_pts = [o_xyz + self.oh_length * d2 for _, d2 in cone]
      cone_res = self._clearances(cone_pts, own_idx, own_slots) if cone else []
      a1 = {m: _align(d1_dirs[m]) for m in viable}
      best_s = None
      best_pair = None
      for idx, (m, d2) in enumerate(cone):
        if not cone_res[idx][1] or not cation_ok(cone_pts[idx]):
          continue
        s = a1[m] + _align(d2)
        if best_s is None or s > best_s:
          best_s = s
          best_pair = (d1_pts[m], cone_pts[idx])
      if best_pair is not None:
        return best_pair

    d1 = None
    h1_acc = None
    for k, i in enumerate(acceptors):
      if acc_res[k][1] and cation_ok(acc_pts[k]):  # ok and away from cations
        d1 = acc_dirs[k]
        h1_acc = i
        break
    if d1 is None:
      sph_dirs = [matrix.col(v).normalize() for v in _WATER_FALLBACK_DIRECTIONS]
      sph_pts = [o_xyz + self.oh_length * dv for dv in sph_dirs]
      cand_dirs = acc_dirs + sph_dirs
      cand_pts = acc_pts + sph_pts
      cand_res = acc_res + self._clearances(sph_pts, own_idx, own_slots)
      best_j = max(
        range(len(cand_dirs)),
        key=lambda j: (cation_ok(cand_pts[j]),
                       cand_res[j][1], cand_res[j][0]))
      d1 = cand_dirs[best_j]
    h1_xyz = o_xyz + self.oh_length * d1

    # Orthonormal frame for the H2 cone (d1, p, q mutually perpendicular)
    p, q = _ortho_frame(d1)

    # H2: directed at the nearest acceptor not used by H1. Sample the cone
    # and rank each angle (away-from-cation, clash-free, acceptor
    # alignment, clearance); acceptor alignment only counts once the
    # higher-priority constraints are met, and the clearest angle wins if
    # none satisfies them.
    target = None
    for i in acceptors:
      if i == h1_acc:
        continue
      target = accept_dir(i)
      break

    cone_dirs = []
    for k in range(_WATER_CONE_SAMPLES):
      theta = 2.0 * math.pi * k / _WATER_CONE_SAMPLES
      cone_dirs.append(self.cos_hoh * d1
                       + self.sin_hoh * (math.cos(theta) * p + math.sin(theta) * q))
    cone_pts = [o_xyz + self.oh_length * d2 for d2 in cone_dirs]
    cone_res = self._clearances(cone_pts, own_idx, own_slots)

    best_h2_xyz = None
    best_key = None
    for k, d2 in enumerate(cone_dirs):
      min_dist, ok = cone_res[k]
      cat_ok = cation_ok(cone_pts[k])
      good = cat_ok and ok
      align = d2.dot(target) if target is not None else 0.0
      key = (cat_ok, ok, align if good else 0.0, min_dist)
      if best_key is None or key > best_key:
        best_key = key
        best_h2_xyz = cone_pts[k]
    return h1_xyz, best_h2_xyz

  def _apply_sweep(self):
    """Run one relaxation sweep.

    Re-places every recorded water against the current positions of all the
    others (its own H excluded), updating ``placed_coords`` and the atoms in
    place.
    """
    self.placed_tree = KDTree(self.placed_coords) if self.placed_coords else None
    self.pending_tree = None
    self.n_built = len(self.placed_coords)
    for o_xyz, own_idx, slots in self.records:
      own_slots = {slot for _, slot, _ in slots}
      h_xyz = dict(zip((1, 2), self._place_one(o_xyz, own_idx, own_slots)))
      for atom, slot, di in slots:
        self.placed_coords[slot] = tuple(h_xyz[di])
        atom.set_xyz(tuple(h_xyz[di]))

  def _restore(self, coords):
    """Reset all placed H to a coordinate snapshot.

    Parameters
    ----------
    coords : list of tuple
        Per-slot placed-H coordinates to restore (a snapshot of
        ``placed_coords``).
    """
    for _, _, slots in self.records:
      for atom, slot, di in slots:
        self.placed_coords[slot] = coords[slot]
        atom.set_xyz(coords[slot])

  def _clashing_records(self):
    """Find records whose placed H still clash.

    Returns
    -------
    list of int
        Indices into ``records`` of waters with a placed H that fails the
        clash gate.
    """
    self.placed_tree = KDTree(self.placed_coords) if self.placed_coords else None
    self.pending_tree = None
    self.n_built = len(self.placed_coords)
    bad = []
    for ri, (o_xyz, own_idx, slots) in enumerate(self.records):
      own_slots = {slot for _, slot, _ in slots}
      pts = [matrix.col(self.placed_coords[slot]) for _, slot, _ in slots]
      if any(not ok for _, ok in self._clearances(pts, own_idx, own_slots)):
        bad.append(ri)
    return bad

  def run(self):
    """Place H on every bare water, refine, optionally basin-hop, keep best.

    Modifies the hierarchy in place. Returns the kept-state label (see
    :func:`place_water_hydrogens`).
    """
    hier = self.hier
    if self.reorient_existing:
      _strip_water_hydrogens(hier)

    atoms = list(hier.atoms())
    if not atoms:
      return None
    self.atoms = atoms

    # Static neighbours (protein, ligands, water O, any pre-existing H):
    # fixed, so the KDTree is built once. The water H we place live in a
    # separate, much smaller list/tree that is cheap to rebuild as
    # placement and refinement proceed.
    self.static_coords = [tuple(a.xyz) for a in atoms]
    self.static_tree = KDTree(self.static_coords)
    self.static_is_h = [a.element.strip().upper() in ("H", "D") for a in atoms]
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
    for i, a in enumerate(atoms):
      if not self.static_is_h[i]:
        continue
      for j in self.static_tree.query_ball_point(a.xyz, _WATER_NH_BOND):
        if j != i and atoms[j].element.strip().upper() == "N":
          self.donor_n.add(j)

    # Lone-pair lobe directions per acceptor (opt-in; empty when off).
    self.acc_lobes = _acceptor_lobes(atoms, self.static_tree, self.donor_n) \
        if self.lone_pair_directed else {}

    self.cos_hoh = math.cos(math.radians(_WATER_HOH_DEG))
    self.sin_hoh = math.sin(math.radians(_WATER_HOH_DEG))

    # Gather the waters to protonate. Order them most-crowded first: a water
    # boxed in by many neighbours is the hardest to satisfy, so placing it
    # while few water H are fixed gives it the most freedom; roomy waters have
    # options left over and adapt around it. ``crowd`` is the neighbour count
    # within the clash radius (excluding the water's own atoms).
    waters = []
    for ag in hier.atom_groups():
      if not _is_water(ag.resname):
        continue
      if len(ag.atoms()) >= 3:
        continue  # already protonated
      o = next((a for a in ag.atoms()
                if a.element.strip().upper() == "O"), None)
      if o is None:
        continue
      own_idx = {pos_by_id[a.memory_id()] for a in ag.atoms()
                 if a.memory_id() in pos_by_id}
      crowd = sum(
        1 for j in self.static_tree.query_ball_point(o.xyz, _WATER_CLEARANCE_RADIUS)
        if j not in own_idx)
      waters.append((crowd, ag, o, own_idx))
    waters.sort(key=lambda w: w[0], reverse=True)

    # Initial greedy pass over the ordered waters, each avoiding the H already
    # placed on earlier ones (prefix + pending trees, refreshed below).
    # ``records`` keeps per-water (o_xyz, own_idx, placed-H slots) so the
    # refinement sweeps can re-place each water against the *final*
    # environment, breaking the order-dependence of the greedy pass.
    self.records = []   # list of (o_xyz, own_idx, [(atom, slot, di), ...])
    for crowd, ag, o, own_idx in waters:
      if self.element is not None:
        proton_element = self.element
      else:
        proton_element = "D" if ag.resname.strip().upper() == "DOD" else "H"
      o_xyz = matrix.col(o.xyz)
      existing_names = {a.name.strip() for a in ag.atoms()}

      # Refresh the placed-H trees. Fold the pending H into the (rebuilt)
      # prefix tree once enough have accumulated; otherwise just rebuild the
      # small pending tree so this water sees the previous waters' H. Their
      # union is every earlier water's H, exactly as a single rebuilt tree.
      if len(self.placed_coords) - self.n_built >= _WATER_TREE_REBUILD:
        self.placed_tree = KDTree(self.placed_coords)
        self.n_built = len(self.placed_coords)
        self.pending_tree = None
      else:
        pending = self.placed_coords[self.n_built:]
        self.pending_tree = KDTree(pending) if pending else None
      h_xyz = dict(zip((1, 2), self._place_one(o_xyz, own_idx, set())))

      slots = []
      for di in (1, 2):
        proton_name = f" {proton_element}{di} "
        if proton_name.strip() in existing_names:
          continue
        atom = _new_h_atom(proton_name, proton_element,
                           tuple(h_xyz[di]), o.occ, o.b)
        ag.append_atom(atom)
        slots.append((atom, len(self.placed_coords), di))
        self.placed_coords.append(tuple(h_xyz[di]))
      if slots:
        self.records.append((o_xyz, own_idx, slots))

    # Refinement sweeps: re-place each water against the *final* set of all
    # placed H (its own excluded), so a water no longer ignores neighbours
    # placed after it. This relaxes the water-water clashes the greedy pass
    # leaves in dense clusters. With refine_tol > 0 it stops early once a
    # sweep reduces the close (<2.0 A) contact count by fewer than refine_tol
    # (n_refine is then a cap); refine_tol == 0 runs all n_refine sweeps.
    # Refinement isn't strictly monotonic, so the best state is kept and
    # restored at the end -- the result is never worse than any state seen.
    stats = _water_clash_stats(hier) \
        if (self.on_state or self.n_refine or self.n_basin) else None
    if self.on_state is not None:
      self.on_state("initial", stats)
    prev_n20 = stats[1] if stats is not None else None
    best_n20 = prev_n20
    best_coords = list(self.placed_coords) if (self.n_refine or self.n_basin) else None
    best_label = "initial"
    for i in range(self.n_refine):
      self._apply_sweep()
      stats = _water_clash_stats(hier)
      if self.on_state is not None:
        self.on_state(f"sweep {i + 1}", stats)
      if best_n20 is None or stats[1] < best_n20:
        best_n20 = stats[1]
        best_coords = list(self.placed_coords)
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
          _kick_water(self.records[ri], self.placed_coords, self.oh_length,
                      self.cos_hoh, self.sin_hoh, rng)
        for s in range(_WATER_BASIN_RELAX):
          self._apply_sweep()
          stats = _water_clash_stats(hier)
          label = f"basin {it + 1}.{s + 1}"
          if self.on_state is not None:
            self.on_state(label, stats)
          if best_n20 is None or stats[1] < best_n20:
            best_n20 = stats[1]
            best_coords = list(self.placed_coords)
            best_label = label

    if best_coords is not None:
      self._restore(best_coords)
    return best_label if (self.n_refine or self.n_basin) else None


def place_water_hydrogens(hier, oh_length=_WATER_OH_NEUTRON, element=None,
                          n_refine=_WATER_REFINE_SWEEPS,
                          refine_tol=_WATER_REFINE_TOL, n_basin=0,
                          reorient_existing=False, lone_pair_directed=False,
                          joint=False, on_state=None):
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
  all other atoms *and all other waters' placed H* (a small KDTree of the
  placed protons, separate from the static-atom tree, so neighbouring
  waters don't aim H at each other). Positions stay at least
  ``_WATER_MIN_CLEARANCE`` from everything else where the local packing
  allows it; in genuinely crowded pockets the clearest available position
  is used.

  Waters are placed most-crowded first (the hardest to satisfy get the
  most freedom), then up to ``n_refine`` refinement sweeps re-place each
  water against the final environment to relax water-water clashes; the
  best sweep is always kept, so the result is never worse than any state
  seen. New H inherit the parent O's occupancy and B factor.

  Modifies *hier* in place. By default it is idempotent: water residues
  already carrying H are left untouched.

  This is a thin wrapper around :class:`_WaterHydrogenPlacer`; the
  parameters and return value are documented here.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model hierarchy; modified in place.
  oh_length : float, optional
      O-H bond length in A (default neutron, 0.984; pass
      ``_WATER_OH_XRAY`` = 0.957 for X-ray-style placement).
  element : str or None, optional
      Element of the placed atoms: ``"H"`` or ``"D"`` forces that element
      on every water; None (default) picks per residue -- ``"D"`` for DOD,
      ``"H"`` for HOH. Atom names follow the element (``H1``/``H2`` or
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
  reorient_existing : bool, optional
      If True, strip any H already on waters and re-place every water from
      scratch; otherwise already-protonated waters are left untouched.
  lone_pair_directed : bool, optional
      If True, each O-H aims at an acceptor's lone-pair lobe (estimated
      from the acceptor's bonded-neighbour geometry) rather than its
      nucleus, for better D-H...A angles.
  joint : bool, optional
      If True, optimize both O-H together (the orientation maximizing the
      summed acceptor alignment of the pair) rather than greedily fixing H1
      then placing H2 -- so a water flanked by two acceptors can donate
      well to both. Falls back to the greedy path where no clash-free pair
      exists.
  on_state : callable, optional
      If given, called as ``on_state(label, stats)`` once the placement
      reaches each state -- ``"initial"`` after the greedy pass, then
      ``"sweep N"`` after each refinement sweep (and ``"basin N.M"`` during
      basin-hopping) -- where ``stats`` is the ``_water_clash_stats``
      tuple. Lets a caller stream a progress table row by row.

  Returns
  -------
  str or None
      The label of the kept state (``"initial"``, ``"sweep N"`` or
      ``"basin N.M"``) when refinement or basin-hopping ran, else None.
  """
  return _WaterHydrogenPlacer(
    hier, oh_length=oh_length, element=element, n_refine=n_refine,
    refine_tol=refine_tol, n_basin=n_basin,
    reorient_existing=reorient_existing,
    lone_pair_directed=lone_pair_directed, joint=joint,
    on_state=on_state).run()


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
      if a.element.strip().upper() in ("H", "D"):
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
      if a.element.strip().upper() in ("H", "D"):
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
  if any(a.element.strip().upper() == "D" for a in hier.atoms()):
    return True, "D atoms present (no conclusive experiment metadata)"
  return False, "no neutron signal (assuming X-ray)"

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

Needs ``scipy`` (KDTree); the vectorised work is ``scitbx.array_family.flex``.
"""

from __future__ import absolute_import, division, print_function

import math
import random

import iotbx.pdb
from libtbx import group_args
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
_EMPTY_SLOTS = flex.size_t()


def _as_vec3(tuples):
  """Build a flex.vec3_double, tolerating the empty case."""
  return flex.vec3_double(tuples) if tuples else flex.vec3_double()


# Index arrays for the (candidate x neighbour) block, keyed by its shape.
# The block is flattened row-major -- candidate i occupies ``[i * m, (i+1) * m)``
# -- so the neighbour axis is contiguous and a row is a plain slice. The two
# arrays only depend on the shape, so they are built once and shared by every
# water with that shape.
_TILE_CACHE = {}


def _tile(k, m):
  """``(cand_index, nbr_index)`` for a flattened (k, m) block."""
  t = _TILE_CACHE.get((k, m))
  if t is None:
    t = (flex.size_t([i for i in range(k) for _ in range(m)]),
         flex.size_t(list(range(m)) * k))
    _TILE_CACHE[(k, m)] = t
  return t


def _row_mins(d, k, m, rows):
  """Per-row minima of a flattened (k, m) block, for selected rows.

  flex has no segmented reduction, so the block is sorted once and each
  row's smallest entry read off the permutation -- one C++ sort whatever
  the number of rows wanted, which beats slicing rows out one at a time
  for everything but the single-row case. Exact either way: no arithmetic
  touches the values, and ``min`` is order-independent.

  Parameters
  ----------
  d : scitbx.array_family.flex.double
      The block, flattened row-major.
  k, m : int
      Block shape.
  rows : scitbx.array_family.flex.size_t
      Rows whose minimum is wanted.

  Returns
  -------
  scitbx.array_family.flex.double
      One minimum per entry of ``rows``.
  """
  n = rows.size()
  if n == 0:
    return flex.double()
  if n == 1:
    i = int(rows[0]) * m
    return flex.double(1, flex.min(d[i:i + m]))
  # Descending order, then scatter into a per-row slot: the last write for a
  # row is its smallest entry, so every row ends up holding its own minimum.
  perm = flex.sort_permutation(d).reversed()
  mi = flex.size_t(k, 0)
  mi.set_selected(perm / m, perm)
  return d.select(mi).select(rows)


def _rank_top(cat, ok):
  """Rows of maximal ``2 * cat + ok`` rank, and that rank.

  Reproduces ``rank = cat * 2 + ok; top = rank == rank.max()`` without
  materialising the rank array: the four ranks are just the four
  combinations of the two booleans, so the highest non-empty one is the
  answer.

  Returns
  -------
  tuple
      ``(top, rank)`` -- a flex.bool over the candidates and the rank value.
  """
  both = cat & ok
  if both.count(True):
    return both, 3
  if cat.count(True):
    return cat, 2          # cat & ~ok, and cat & ok is empty
  if ok.count(True):
    return ok, 1           # ~cat & ok, and cat is empty
  return flex.bool(cat.size(), True), 0


def _cat_ok(cat, pts, o):
  """Whether each candidate is clear of every nearby cation.

  True where the candidate is in the hemisphere away from every nearby
  cation, i.e. ``(H - O).(cation - O) <= 0`` for all of them.

  Parameters
  ----------
  cat : list of tuple
      O->cation unit directions.
  pts : scitbx.array_family.flex.vec3_double
      Candidate H positions.
  o : tuple of float
      Water O coordinates.

  Returns
  -------
  scitbx.array_family.flex.bool
      One entry per candidate.
  """
  ok = flex.bool(pts.size(), True)
  if not cat:
    return ok
  dx, dy, dz = (pts - o).parts()
  for c in cat:
    bad = (dx * c[0] + dy * c[1] + dz * c[2] > 0.0).iselection()
    if bad.size():
      ok.set_selected(bad, False)
  return ok


# The H2 cone is sampled at fixed azimuths, so its cosines and sines are
# constants; likewise the normalized fallback sphere.
_CONE_COS = flex.double([math.cos(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_CONE_SIN = flex.double([math.sin(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_FALLBACK_DIRS = [matrix.col(v).normalize().elems
                  for v in _WATER_FALLBACK_DIRECTIONS]
_FALLBACK_FX = _as_vec3(_FALLBACK_DIRS)


def _ortho_frame(d1):
  """Build an orthonormal frame around a direction.

  Mirrors ``scitbx.matrix.col.ortho()`` operation for operation.

  Parameters
  ----------
  d1 : tuple of float
      Unit direction (the O-H1 axis).

  Returns
  -------
  tuple of tuple
      Two unit vectors ``p, q`` with ``(d1, p, q)`` mutually perpendicular
      -- the frame for sampling the H2 cone around ``d1``.
  """
  x, y, z = d1[0], d1[1], d1[2]
  a, b, c = abs(x), abs(y), abs(z)
  if c <= a and c <= b:
    p = (-y, x, 0.0)
  elif b <= a and b <= c:
    p = (-z, 0.0, x)
  else:
    p = (0.0, -z, y)
  n = math.sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2])
  p = (p[0] / n, p[1] / n, p[2] / n)
  return p, (y * p[2] - p[1] * z,
             z * p[0] - p[2] * x,
             x * p[1] - p[0] * y)


def _rand_unit(rng):
  """A random unit vector (rejection-free Gaussian normalization).

  Parameters
  ----------
  rng : random.Random
      Seeded RNG (for reproducibility).

  Returns
  -------
  tuple of float
      A uniformly random unit vector.
  """
  while True:
    v = (rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1))
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n > 1e-6:
      return (v[0] / n, v[1] / n, v[2] / n)


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


class _Clearance(object):
  """The result of one clearance test (see :meth:`_WaterHydrogenPlacer._clear`).

  Holds the pass/fail flag per candidate and the squared-distance blocks it
  came from, so the clearance *distance* -- needed only to break ties among
  the top-ranked candidates -- can be reduced out of them on demand instead
  of for every candidate. A block already reduced to its per-candidate
  minima (the cached static block of an invariant candidate set) is carried
  in ``mins`` instead.
  """

  __slots__ = ("ok", "blocks", "mins", "k")

  def __init__(self, k):
    self.ok = flex.bool(k, True)
    self.blocks = []   # (flattened (k, m) squared-distance block, m)
    self.mins = []     # per-candidate minima of already-reduced blocks
    self.k = k

  def best_at(self, rows):
    """Clearance distance of the given candidates.

    Parameters
    ----------
    rows : scitbx.array_family.flex.size_t
        Candidate indices.

    Returns
    -------
    scitbx.array_family.flex.double
        Distance to the nearest non-own atom within
        ``_WATER_CLEARANCE_RADIUS``, or the search radius if nothing is
        near, one entry per requested candidate.
    """
    if rows.size() == 0:
      return flex.double()
    best = flex.double(rows.size(), _WATER_CLEARANCE_RADIUS ** 2)
    for mins in self.mins:
      v = mins.select(rows)
      s = v < best
      best.set_selected(s, v.select(s))
    for d, m in self.blocks:
      v = _row_mins(d, self.k, m, rows)
      s = v < best
      best.set_selected(s, v.select(s))
    return flex.sqrt(best)


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
    # placed_xyz is the same data as a flex.vec3_double.
    self.placed_coords = []
    self.placed_xyz = None
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
    cands : scitbx.array_family.flex.vec3_double
        Candidate H positions.

    Returns
    -------
    tuple
        ``(d, m, ok)``: the flattened (k, m) squared-distance block, its
        neighbour count, and whether each candidate clears every one of
        them. ``m`` is zero (and ``d`` None) when the water has no static
        neighbour at all.
    """
    k = cands.size()
    P = self.w_sxyz[wi]
    m = P.size()
    ok = flex.bool(k, True)
    if not m:
      return None, 0, ok
    ik, im = _tile(k, m)
    # Squared distances go through parts() and pow2, never .dot() or
    # .norms(): those contract to fused multiply-add in C++, while the engine
    # this reproduces rounds every step, so they differ in the last ulp and
    # flip threshold tests. flex.pow2 is also ~11x faster than dx * dx. The
    # same applies to every other squared distance in this file.
    dx, dy, dz = (cands.select(ik) - P.select(im)).parts()
    d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
    bad = (d < self.w_sthr[wi].select(im)).iselection()
    if bad.size():
      ok.set_selected(bad / m, False)
    return d, m, ok

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
    cands : scitbx.array_family.flex.vec3_double
        Candidate H positions; an invariant set (the acceptor directions or
        the fallback sphere).
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
      d, m, ok = self._static_block(wi, cands)
      k = cands.size()
      st = g[key] = (ok, _row_mins(d, k, m, flex.size_t_range(k)) if m
                     else None)
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
    cands : scitbx.array_family.flex.vec3_double
        Candidate H positions.
    nbr_slots : scitbx.array_family.flex.size_t
        Placed-H slots that may lie near this water (its own excluded).
    static_key : str or None, optional
        Geometry-dict slot holding the cached static half of the test, for
        a candidate set that does not move (the acceptor directions, the
        fallback sphere). None for the cone, which is rebuilt around each
        pass's O-H1 axis and so has nothing to reuse.

    Returns
    -------
    _Clearance
        Carries ``ok`` -- whether each candidate clears every heavy atom by
        ``_WATER_MIN_CLEARANCE`` and every hydrogen (static or placed) by
        ``_WATER_MIN_H_CLEARANCE`` -- and the squared-distance blocks the
        clearance distances are read back from. Those distances (the
        distance to the nearest non-own atom within
        ``_WATER_CLEARANCE_RADIUS``, or the search radius if nothing is
        near) are only ever wanted for the few top-ranked candidates, so
        they are left to :meth:`_Clearance.best_at` rather than reduced for
        every candidate here.
    """
    # Squared distances throughout: the thresholds and the cap are exact
    # squares and sqrt is monotone, so both the comparisons and the
    # per-candidate minimum are unchanged.
    k = cands.size()
    cl = _Clearance(k)
    if static_key is not None:
      ok, mins = self._static_clear(wi, cands, static_key)
      cl.ok = ok.deep_copy()
      if mins is not None:
        cl.mins.append(mins)
    else:
      d, m, cl.ok = self._static_block(wi, cands)
      if m:
        cl.blocks.append((d, m))
    m = nbr_slots.size()
    if m:
      ik, im = _tile(k, m)
      ct = cands.select(ik)
      dx, dy, dz = (ct - self.placed_xyz.select(nbr_slots).select(im)).parts()
      d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
      bad = (d < _WATER_MIN_H_CLEARANCE ** 2).iselection()
      if bad.size():
        cl.ok.set_selected(bad / m, False)
      cl.blocks.append((d, m))
    return cl

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
      X = self.placed_xyz[:ns].concatenate(self.wh_xyz)
      W = self.slot_wid[:ns].concatenate(self.wh_wid)
    else:
      X, W = self.wh_xyz, self.wh_wid
    n = X.size()
    if n < 2:
      return n, 0, 0, 0, None
    pairs = KDTree(X.as_numpy_array()).query_pairs(2.0, output_type="ndarray")
    if not len(pairs):
      return n, 0, 0, 0, None
    # KDTree hands back an (n, 2) ndarray; take it into flex here and stay
    # there. flex flattens it row-major, so the pair members are the even
    # and odd entries.
    p = flex.size_t(pairs)
    i = p.select(flex.size_t_range(0, p.size(), 2))
    j = p.select(flex.size_t_range(1, p.size(), 2))
    keep = (W.select(i) != W.select(j)).iselection()  # drop same-water H
    if not keep.size():
      return n, 0, 0, 0, None
    i = i.select(keep)
    j = j.select(keep)
    dx, dy, dz = (X.select(i) - X.select(j)).parts()
    d = flex.sqrt(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz))
    return (n, d.size(), (d < 1.8).count(True), (d < 1.5).count(True),
            flex.min(d))

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
        ``acc_dirs``/``acc_pts`` (O-H unit directions as tuples and the H
        positions built from them), ``cat`` (O->cation unit directions),
        ``sph`` (the lazily built fallback-sphere H positions) and
        ``acc_static``/``sph_static`` (the cached static half of the
        clearance test for those two candidate sets, see
        :meth:`_static_clear`).
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

    o = o_xyz.elems
    acc_dirs = [accept_dir(i).elems for i in acceptors]
    g = {"o": o,
         "acceptors": acceptors,
         "acc_dirs": acc_dirs,
         "acc_pts": _as_vec3(acc_dirs) * self.oh_length + o,
         "cat": [c.elems for c in cation_dirs],
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
    nbr_slots : scitbx.array_family.flex.size_t
        Placed-H slots that may lie near this water (its own excluded).
    fixed_d1 : tuple of float or None, optional
        Unit O-H1 direction to hold fixed instead of searching for one,
        used when completing a water that already carries one proton. The
        returned ``h1_xyz`` then just restates that direction at the
        canonical bond length and the caller ignores it; only ``h2_xyz``
        is new.

    Returns
    -------
    tuple of tuple
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
      acc_cl = self._clear(wi, acc_pts, nbr_slots, "acc_static")
      acc_ok = acc_cl.ok
      acc_cat = _cat_ok(cat, acc_pts, o)

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
        a = acc_dirs[k]
        align = a[0] * d1[0] + a[1] * d1[1] + a[2] * d1[2]
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
          g["sph"] = _FALLBACK_FX * self.oh_length + o
        sph_pts = g["sph"]
        s_cl = self._clear(wi, sph_pts, nbr_slots, "sph_static")
        s_cat = _cat_ok(cat, sph_pts, o)
        if na:
          cand_ok = acc_ok.concatenate(s_cl.ok)
          cand_cat = acc_cat.concatenate(s_cat)
        else:
          cand_ok, cand_cat = s_cl.ok, s_cat
        # (cation-ok, clash-free) as one rank, then clearance, first wins.
        top = _rank_top(cand_cat, cand_ok)[0]
        rows = top.iselection()
        if na:
          # The candidates are the acceptors then the fallback sphere; each
          # half reads its clearances back from its own test.
          lo = rows.select(rows < na)
          hi = rows.select(~(rows < na)) - na
          vals = acc_cl.best_at(lo).concatenate(s_cl.best_at(hi))
          b = flex.max_index(vals)
          pick = int(lo[b]) if b < lo.size() else na + int(hi[b - lo.size()])
        else:
          pick = int(rows[flex.max_index(s_cl.best_at(rows))])
        d1 = acc_dirs[pick] if pick < na else _FALLBACK_DIRS[pick - na]
    h1_xyz = (o[0] + self.oh_length * d1[0],
              o[1] + self.oh_length * d1[1],
              o[2] + self.oh_length * d1[2])

    # Orthonormal frame for the H2 cone (d1, p, q mutually perpendicular)
    p, q = _ortho_frame(d1)

    # H2: rank the cone angles by (away-from-cation, clash-free, best
    # alignment over the acceptors H1 did not take, clearance).
    ns = _WATER_CONE_SAMPLES
    cone_dirs = ((flex.vec3_double(ns, p) * _CONE_COS
                  + flex.vec3_double(ns, q) * _CONE_SIN) * self.sin_hoh
                 + (self.cos_hoh * d1[0], self.cos_hoh * d1[1],
                    self.cos_hoh * d1[2]))
    cone_pts = cone_dirs * self.oh_length + o
    c_cl = self._clear(wi, cone_pts, nbr_slots)
    c_cat = _cat_ok(cat, cone_pts, o)

    top, rank = _rank_top(c_cat, c_cl.ok)
    # The alignment term only ever separates candidates that are both
    # cation-ok and clash-free (rank 3); below that it is zero for every
    # candidate still in the running and cannot break any tie, so the dot
    # products are not worth computing.
    use = [k for k in range(na) if k != h1_k]
    if rank == 3 and use:
      dx, dy, dz = cone_dirs.parts()
      align = None
      for k in use:
        a = acc_dirs[k]
        dot = dx * a[0] + dy * a[1] + dz * a[2]
        if align is None:
          align = dot
        else:
          s = dot > align
          align.set_selected(s, dot.select(s))
      top = top & (align == flex.max(align.select(top)))
    rows = top.iselection()
    return h1_xyz, cone_pts[int(rows[flex.max_index(c_cl.best_at(rows))])]

  def _store(self, slots, h1, h2):
    """Write one water's new H positions to the model and the arrays."""
    for atom, slot, di in slots:
      xyz = h1 if di == 1 else h2
      self.placed_coords[slot] = xyz
      self.placed_xyz[slot] = xyz
      atom.set_xyz(xyz)

  def _apply_sweep(self):
    """Run one relaxation sweep.

    Re-places every recorded water against the current positions of all the
    others (its own H excluded), updating ``placed_coords`` and the atoms in
    place. A water placed later in the sweep therefore sees the H the
    earlier ones just moved, at the position they now hold.
    """
    for wi, slots, fixed_d1 in self.records:
      self._store(slots, *self._place_one(wi, self.w_nbr_slots[wi], fixed_d1))

  def _restore(self, coords):
    """Reset all placed H to a coordinate snapshot.

    Parameters
    ----------
    coords : list of tuple
        Per-slot placed-H coordinates to restore (a snapshot of
        ``placed_coords``).
    """
    for _wi, slots, _fixed in self.records:
      for atom, slot, di in slots:
        self.placed_coords[slot] = coords[slot]
        self.placed_xyz[slot] = coords[slot]
        atom.set_xyz(coords[slot])

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
      pts = self.placed_xyz.select(
        flex.size_t([slot for _, slot, _ in slots]))
      if not self._clear(wi, pts, self.w_nbr_slots[wi]).ok.all_eq(True):
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
    ct, st = math.cos(theta), math.sin(theta)
    d2 = tuple(self.cos_hoh * d1[c] + self.sin_hoh * (ct * p[c] + st * q[c])
               for c in range(3))
    self._store(slots,
                tuple(o[c] + self.oh_length * d1[c] for c in range(3)),
                tuple(o[c] + self.oh_length * d2[c] for c in range(3)))

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
    # tracked separately, by slot, in placed_coords/placed_xyz.
    self.static_xyz = sel.extract_xyz()
    self.static_tree = KDTree(self.static_xyz.as_numpy_array())
    self.static_is_h = flex.bool([a.element_is_hydrogen() for a in atoms])
    self.static_thr = flex.double(len(atoms), _WATER_MIN_CLEARANCE ** 2)
    self.static_thr.set_selected(self.static_is_h,
                                 _WATER_MIN_H_CLEARANCE ** 2)
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
    h_idx = self.static_is_h.iselection()
    if h_idx.size():
      for i, nbrs in zip(h_idx, self.static_tree.query_ball_point(
          self.static_xyz.select(h_idx).as_numpy_array(), _WATER_NH_BOND,
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
            fixed_d1 = v.normalize().elems
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
    self.wh_xyz = _as_vec3(wh_xyz)
    self.wh_wid = flex.size_t(wh_wid) if wh_wid else flex.size_t()

    # Per-water constants. Neither the static atoms nor the water O move, so
    # every neighbour list the placement needs is built once, here, and
    # reused by the greedy pass, every relaxation sweep and every basin
    # round.
    n = len(waters)
    self.w_geom = [None] * n
    self.placed_coords = []
    self.placed_xyz = flex.vec3_double(2 * n, (0.0, 0.0, 0.0))
    self.slot_wid = flex.size_t(2 * n, 0)
    self.records = []   # (water index, [(atom, slot, di), ...], fixed_d1)
    self.w_wnbr = []
    if n:
      o_pts = flex.vec3_double([w[1].xyz for w in waters])
      # Order most-crowded first: a water boxed in by many neighbours is the
      # hardest to satisfy, so placing it while few water H are fixed gives
      # it the most freedom; roomy waters have options left over and adapt
      # around it. ``crowd`` is the neighbour count within the clash radius
      # (excluding the water's own atoms).
      crowd = [sum(1 for j in nb if j not in waters[k][2]) for k, nb in
               enumerate(self.static_tree.query_ball_point(
                 o_pts.as_numpy_array(), _WATER_CLEARANCE_RADIUS,
                 return_sorted=False))]
      order = sorted(range(n), key=lambda k: crowd[k], reverse=True)
      waters = [waters[k] for k in order]
      o_pts = o_pts.select(flex.size_t(order))
      o_np = o_pts.as_numpy_array()
      self.w_own = [w[2] for w in waters]
      self.w_o_col = [matrix.col(w[1].xyz) for w in waters]
      # query_ball_point leaves a single point's indices unsorted but sorts
      # a batch's; the acceptor order is the tie-break of the nearest-first
      # sort, so the batch must not sort either.
      self.w_acc_raw = self.static_tree.query_ball_point(
        o_np, _WATER_ACCEPTOR_RADIUS, return_sorted=False)
      self.w_cat_raw = self.static_tree.query_ball_point(
        o_np, _WATER_CATION_RADIUS, return_sorted=False)
      # One static-neighbour block per water: every candidate H lies on the
      # O-H sphere about the O, so a single ball of oh_length + clearance
      # covers every candidate's own clearance ball.
      self.w_sxyz = []
      self.w_sthr = []
      for wi, nb in enumerate(self.static_tree.query_ball_point(
          o_np, self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
          return_sorted=False)):
        own = self.w_own[wi]
        idx = flex.size_t([int(j) for j in nb if j not in own])
        self.w_sxyz.append(self.static_xyz.select(idx))
        self.w_sthr.append(self.static_thr.select(idx))
      # A placed H sits within oh_length of its own O, so only waters whose
      # O lie within clearance + 2 oh_length can ever hold one near this
      # water's candidates. A neighbour map over the (fixed) water O
      # therefore replaces the growing KDTree over the placed H outright.
      self.w_wnbr = [[j for j in nb if j != wi] for wi, nb in
                     enumerate(KDTree(o_np).query_ball_point(
                       o_np,
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
      return flex.size_t(nbr) if nbr else _EMPTY_SLOTS

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
        xyz = h1 if di == 1 else h2
        atom = _new_h_atom(proton_name, proton_element, xyz, o.occ, o.b,
                           o.hetero)
        ag.append_atom(atom)
        slot = len(self.placed_coords)
        slots.append((atom, slot, di))
        self.placed_coords.append(xyz)
        self.placed_xyz[slot] = xyz
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
    best_coords = list(self.placed_coords) if (self.n_refine or self.n_basin) else None
    best_label = "initial"
    for i in range(self.n_refine):
      self._apply_sweep()
      stats = self._stats()
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
          self._kick(ri, rng)
        for s in range(_WATER_BASIN_RELAX):
          self._apply_sweep()
          stats = self._stats()
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

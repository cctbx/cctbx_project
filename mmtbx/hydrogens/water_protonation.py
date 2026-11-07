"""H-bond-aware placement of the two hydrogens on bare water oxygens.

For every bare water oxygen (any common water residue: HOH, DOD, H2O, WAT,
OH2, ...) the two H are placed pointing at H-bond acceptors, clear of the
whole structure (including H placed on other waters) and out of the
hemisphere of nearby metal cations. Geometry only: no map, no monomer
library. O-H and H-O-H are cctbx's restraint targets for water (0.850 A for
X-ray, 0.980 A for neutron data, 103.91 deg), or on request the gas-phase
molecule (0.957 A, 104.5 deg). Given a crystal symmetry, waters at a lattice
contact also see the neighbouring asymmetric units.

The public entry point :func:`place_water_hydrogens` modifies a hierarchy
in place; :class:`mmtbx.programs.water_protonation.Program` wraps it as the
``mmtbx.development.virvel`` command line.

Needs ``scipy`` (KDTree); the vectorised work is
``scitbx.array_family.flex``.
"""

import math

import iotbx.pdb
from cctbx import crystal
from cctbx.crystal import super_cell
from iotbx.pdb.utils import check_for_missing_elements
from libtbx import group_args
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import fibonacci_sphere_points


# ---------------------------------------------------------------------------
# Constants for the H-bond-aware water-H placer (``place_water_hydrogens``).
# ---------------------------------------------------------------------------

# Water geometry, (O-H in A, H-O-H in deg): cctbx's restraint targets for HOH
# (chem_data/geostd/h/data_HOH.cif), O-H to the H electron density for X-ray
# and to the nucleus for neutron data; and the isolated molecule.
_WATER_GEOMETRY = {
  "xray":      (0.850, 103.91),
  "neutron":   (0.980, 103.91),
  "gas_phase": (0.957, 104.5),
}
_WATER_ACCEPTOR_RADIUS = 3.5  # max distance to search for H-bond acceptors (A)
_WATER_ACCEPTOR_ELEMENTS = frozenset({"O", "N", "F", "S", "CL"})
_WATER_NH_BOND = 1.3          # max N-H distance for "N carries an H" (A)
_WATER_N_BOND = 1.8           # max bond length from an N to a heavy atom (A)
_WATER_PLANAR_N = 350.0       # min bond-angle sum (deg) of a planar N
_WATER_CONE_SAMPLES = 36      # angular samples around the O-H1 cone
# Element-aware "clash-free" thresholds: a candidate H must clear every
# heavy atom by _WATER_MIN_CLEARANCE and every hydrogen by the larger
# _WATER_MIN_H_CLEARANCE.
_WATER_MIN_CLEARANCE = 1.5
_WATER_MIN_H_CLEARANCE = 2.0
_WATER_CLEARANCE_RADIUS = 3.0  # neighbour search radius for clearance (A)
# Distance under which a symmetry equivalent is the same atom, not a copy.
_WATER_SYM_EQUIV_TOL = 0.5

# Basin-hopping (opt-in): each round restarts from the best state, randomly
# re-orients the still-clashing waters, and relaxes. Seeded.
_WATER_BASIN_SEED = 0
_WATER_BASIN_RELAX = 2        # relaxation sweeps after each random kick
# Cations within _WATER_CATION_RADIUS of the water O are repellers, keeping
# both H in the hemisphere away from them, rather than acceptors.
_WATER_CATION_ELEMENTS = frozenset({
  "LI", "NA", "K", "RB", "CS",                                # alkali
  "BE", "MG", "CA", "SR", "BA",                               # alkaline earth
  "SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN",  # 3d
  "Y", "ZR", "NB", "MO", "TC", "RU", "RH", "PD", "AG", "CD",  # 4d
  "HF", "TA", "W", "RE", "OS", "IR", "PT", "AU", "HG",        # 5d
  "AL", "GA", "IN", "SN", "SB", "TL", "PB", "BI",             # p-block
  "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",             # lanthanides
  "TB", "DY", "HO", "ER", "TM", "YB", "LU",
  "TH", "U",                                                  # actinides
})
_WATER_CATION_RADIUS = 3.0
_WATER_METAL_COORD_RADIUS = 2.6  # max first-shell metal bond (A)
_WATER_H1_ACC_ALIGN = 0.7        # min cos(O-H, O->acceptor) to donate to it
# Cluster search (opt-in): a proton scores its distance to the nearest
# acceptor, at most _WATER_CLUSTER_CAP (A); a cluster with up to
# _WATER_CLUSTER_EXHAUSTIVE combinations of candidates is searched in full.
_WATER_CLUSTER_CAP = 3.0
_WATER_CLUSTER_EXHAUSTIVE = 100000


def _has_deuterium(hier):
  """True if the model carries any D atom or DOD water."""
  if (hier.atoms().extract_element(strip=True) == "D").count(True):
    return True
  return any(ag.resname.strip().upper() == "DOD" for ag in hier.atom_groups())


def count_environment_hydrogens(hier):
  """Count the H and D atoms outside the water residues.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Model to count over.

  Returns
  -------
  int
      Number of H and D atoms in non-water atom groups.
  """
  return sum(_n_hd(ag) for ag in hier.atom_groups()
             if not _is_water(ag.resname))


def _hd_flags(ag):
  """``(atoms array, per-atom H/D flags)`` for one atom_group."""
  # Callers index this array: an af_shared_atom has no __iter__.
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


def _tile(k, m):
  """``(cand_index, nbr_index)`` for a flattened (k, m) block.

  The block is flattened row-major: candidate i occupies
  ``[i * m, (i+1) * m)``, so a row is a plain slice.
  """
  r = flex.size_t_range(k * m)
  return r / m, r % m


def _row_mins(d, k, m, rows):
  """Per-row minima of a flattened (k, m) block, for the rows in ``rows``;
  an empty row's minimum is +inf."""
  n = rows.size()
  if n == 0:
    return flex.double()
  if m == 0:
    return flex.double(n, math.inf)
  if n == 1:
    i = int(rows[0]) * m
    return flex.double(1, flex.min(d[i:i + m]))
  # Descending order, then scatter into a per-row slot: the last write for a
  # row is its smallest entry.
  perm = flex.sort_permutation(d).reversed()
  mi = flex.size_t(k, 0)
  mi.set_selected(perm / m, perm)
  return d.select(mi).select(rows)


def _rank_top(cat, ok):
  """Candidates of maximal ``2 * cat + ok`` rank, as ``(top, rank)``."""
  both = cat & ok
  if both.count(True):
    return both, 3
  if cat.count(True):
    return cat, 2          # cat & ~ok, and cat & ok is empty
  if ok.count(True):
    return ok, 1           # ~cat & ok, and cat is empty
  return flex.bool(cat.size(), True), 0


def _cat_ok(cat, pts, o):
  """Whether each candidate lies in the hemisphere away from each cation.

  True where ``(H - O).(cation - O) <= 0`` for every O->cation direction in
  ``cat``, with ``pts`` the candidate H positions and ``o`` the water O.
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


# Fixed H2-cone azimuth cosines and sines, and the normalized fallback sphere.
_CONE_COS = flex.double([math.cos(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_CONE_SIN = flex.double([math.sin(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_FALLBACK_DIRS = [matrix.col(v).normalize().elems
                  for v in fibonacci_sphere_points(64)]
_FALLBACK_FX = flex.vec3_double(_FALLBACK_DIRS)


def _split_dirs(a1, a2, cos_half, sin_half):
  """O-H directions equally far off two acceptors.

  Both lie in the plane of the unit acceptor directions ``a1`` and ``a2``,
  symmetric about their bisector and H-O-H apart (``cos_half``/``sin_half``
  of half that angle), so each proton is equally far off its acceptor.
  None when the acceptors are (anti)collinear and the plane is undefined.
  """
  b = (a1[0] + a2[0], a1[1] + a2[1], a1[2] + a2[2])
  bn = math.sqrt(b[0] * b[0] + b[1] * b[1] + b[2] * b[2])
  if bn < 1e-6:
    return None
  b = (b[0] / bn, b[1] / bn, b[2] / bn)
  t = a1[0] * b[0] + a1[1] * b[1] + a1[2] * b[2]
  w = (a1[0] - t * b[0], a1[1] - t * b[1], a1[2] - t * b[2])
  wn = math.sqrt(w[0] * w[0] + w[1] * w[1] + w[2] * w[2])
  if wn < 1e-6:
    return None
  w = (w[0] / wn, w[1] / wn, w[2] / wn)
  return (tuple(b[c] * cos_half + w[c] * sin_half for c in range(3)),
          tuple(b[c] * cos_half - w[c] * sin_half for c in range(3)))


def _water_residue_groups(hier):
  """Each residue group holding a water, as ``(residue_group, water atom
  groups)``."""
  for rg in hier.residue_groups():
    ags = [ag for ag in rg.atom_groups() if _is_water(ag.resname)]
    if ags:
      yield rg, ags


def _water_conformers(rg, ags):
  """The water conformers of residue group ``rg``, whose water atom groups
  are ``ags``.

  Returns ``(altloc, atoms, target)`` per conformer: ``atoms`` are its blank
  atoms plus its altloc's, and ``target`` the atom group its new H join. A
  residue without altlocs is one conformer per atom group, as it stands; one
  with altlocs has one per altloc, from ``residue_group.conformers()``,
  leaving out another residue's (HOH in A, SO4 in B).
  """
  if not rg.have_conformers():
    return [(ag.altloc, ag.atoms(), ag) for ag in ags]
  by_altloc = {ag.altloc: ag for ag in ags}
  blank = by_altloc.get("")
  out = []
  seen = set()
  for cf in rg.conformers():
    res = cf.residues()[0]
    if not _is_water(res.resname):
      continue
    target = by_altloc.get(cf.altloc, blank)
    if target is None or target.memory_id() in seen:
      continue
    seen.add(target.memory_id())
    out.append((target.altloc, res.atoms(), target))
  return out


def _water_o_sites(ags):
  """``{altloc: O site}`` over the water atom groups ``ags`` holding an O."""
  sites = {}
  for ag in ags:
    ats = ag.atoms()
    els = ats.extract_element(strip=True)
    for k in range(ats.size()):
      if els[k].upper() == "O":
        sites[ag.altloc] = ats[k].xyz
        break
  return sites


def _h_anchor(xyz, altloc, o_sites):
  """The O site a water H at ``xyz`` in ``altloc`` is bound to: its own atom
  group's (``o_sites`` from :func:`_water_o_sites`), else the nearest of its
  residue's, as for a blank H on a water split between altlocs; None
  without one."""
  if altloc in o_sites:
    return o_sites[altloc]
  if not o_sites:
    return None
  h = matrix.col(xyz)
  return min(o_sites.values(), key=lambda o: (matrix.col(o) - h).length())


def _strip_water_hydrogens(hier):
  """Remove every H/D from water residues.

  Returns what was removed, as ``{atom_group memory_id: [(element,
  occupancy), ...]}`` over the water atom groups that carried any.
  """
  stripped = {}
  for ag in hier.atom_groups():
    if not _is_water(ag.resname):
      continue
    ats, hd = _hd_flags(ag)
    removed = [ats[k] for k in range(ats.size()) if hd[k]]
    if not removed:
      continue
    stripped[ag.memory_id()] = [(a.element.strip().upper(), a.occ)
                                for a in removed]
    for a in removed:
      ag.remove_atom(a)
  return stripped


def _free_proton_name(existing_names, element):
  """First of ``<element>1`` / ``<element>2`` not in ``existing_names``.

  Returns a 4-character PDB atom name, or None if both are taken.
  """
  for di in (1, 2):
    name = f"{element}{di}"
    if name not in existing_names:
      return f" {name} "
  return None


def _compatible(ca, cb):
  """Whether atoms of conformer indices ``ca`` and ``cb`` coexist.

  cctbx's rule: atoms in two different altlocs never meet, and a blank atom
  (index 0) meets every altloc.
  """
  return not ca or not cb or ca == cb


def _symmetry_environment(hier, crystal_symmetry, radius):
  """Copies of the atoms crystal symmetry places within ``radius`` of a water.

  Each copy carries its whole residue group (see
  :meth:`iotbx.pdb.hierarchy.root.sym_equiv_residue_groups`). Copies are
  environment only; the placed protons of symmetry-equivalent waters are
  tracked separately (see :func:`_water_sym_equiv_neighbours`). A copied
  atom closer than ``_WATER_SYM_EQUIV_TOL`` to an atom of the same element,
  in the model or an earlier copy, counts as that atom and is left out.

  Returns ``(hierarchies, atoms, xyz, source)``: the sub-hierarchies, which
  own the atom objects and must be kept alive, their atoms in coordinate
  order, their sites, and the index of the model atom each was copied from.
  All four are empty when symmetry places nothing in range.
  """
  o_sel = flex.size_t([a.i_seq for a in hier.atoms()
                       if _is_water(a.parent().resname)
                       and a.element.strip().upper() == "O"])
  hiers = []
  atoms = []
  xyz = flex.vec3_double()
  source = flex.size_t()
  for equiv in hier.sym_equiv_residue_groups(
      crystal_symmetry, radius, selection=o_sel,
      min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
    sub_atoms = equiv.hierarchy.atoms()
    hiers.append(equiv.hierarchy)
    atoms.extend(list(sub_atoms))
    xyz.extend(sub_atoms.extract_xyz())
    source.extend(equiv.i_seqs)
  return hiers, atoms, xyz, source


def _sym_equiv_pairs(sites, crystal_symmetry, radius):
  """Every ``(i, j, op)`` whose equivalent of site ``j`` lies within
  ``radius`` of site ``i``.

  ``op`` acts on fractional coordinates. A contact is listed from both ends,
  ``(j, i)`` under the inverse operator, except that a site on a special
  position can see more equivalents of its partner than the partner sees of
  it. An equivalent closer than ``_WATER_SYM_EQUIV_TOL`` to its own site is
  that same site, and is dropped.
  """
  if not sites.size():
    return []
  tables = super_cell.get_sym_equiv_tables(
    sites_cart=sites, crystal_symmetry=crystal_symmetry, radius=radius,
    min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL)
  found = []
  for i in range(sites.size()):
    siiu, _ = super_cell.get_siiu(
      crystal_symmetry=crystal_symmetry, selection=[i],
      symmetry_tables=tables)
    for j, ops in siiu.items():
      for op in ops:
        found.append((i, j, op))
  return found


def _water_sym_equiv_neighbours(sites, crystal_symmetry, radius, conf=None):
  """Symmetry equivalents of ``sites`` that fall within ``radius`` of a site.

  ``conf`` holds per-site conformer indices; a pair whose sites never
  coexist (see :func:`_compatible`) is dropped. None keeps every pair.

  Returns ``(by_site, transforms, own)``: per site, the sorted ``(block, other
  site)`` pairs whose equivalent is in range; per block the cartesian
  ``(rotation, translation)`` that produces it, block 0, the untransformed
  one, holding None; and per site the blocks whose operator brings the
  site's own equivalent into range. An equivalent closer than
  ``_WATER_SYM_EQUIV_TOL`` to its own site is that same site, and is
  dropped.
  """
  n = sites.size()
  by_site = [[] for _ in range(n)]
  own = [[] for _ in range(n)]
  transforms = [None]
  found = _sym_equiv_pairs(sites, crystal_symmetry, radius)
  if conf is not None:
    found = [(wi, wj, op) for wi, wj, op in found
             if _compatible(conf[wi], conf[wj])]
  if not found:
    return by_site, transforms, own
  # Blocks numbered by operator, independent of the order the tables list
  # pairs in.
  unit_cell = crystal_symmetry.unit_cell()
  blocks = {}
  for key, op in sorted({str(op): op for _, _, op in found}.items()):
    blocks[key] = len(blocks) + 1
    transforms.append((
      matrix.sqr(unit_cell.matrix_cart(op.r())),
      matrix.col(unit_cell.orthogonalize(op.t().as_double()))))
  for wi, wj, op in found:
    if wi == wj:
      own[wj].append(blocks[str(op)])
    else:
      by_site[wi].append((blocks[str(op)], wj))
  for pairs in by_site:
    pairs.sort()
  for pairs in own:
    pairs.sort()
  return by_site, transforms, own


def _donating(pts, o, acceptor, oh_length):
  """Which candidates aim a proton at an acceptor within the donation limit.

  ``pts`` holds the candidates' protons as consecutive pairs about the O
  ``o``; ``o`` and ``acceptor`` are ``matrix.col``. Returns one flag per
  candidate.
  """
  u = (acceptor - o).normalize()
  dx, dy, dz = (pts - o.elems).parts()
  near = (dx / oh_length * u[0] + dy / oh_length * u[1]
          + dz / oh_length * u[2]) >= _WATER_H1_ACC_ALIGN
  k = pts.size()
  return (near.select(flex.size_t_range(0, k, 2))
          | near.select(flex.size_t_range(1, k, 2)))


def _unrepeated_pairs(pts):
  """Proton pairs without those repeating an earlier one.

  ``pts`` holds the pairs as consecutive protons. A pair repeats an earlier
  one when its two protons lie within 1e-6 A of the earlier pair's in all,
  matched either way round. Returns the pairs kept, in order.
  """
  n = pts.size() // 2
  ia, ib = _tile(n, n)

  def apart(p, q):
    dx, dy, dz = (pts.select(p) - pts.select(q)).parts()
    return flex.sqrt(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz))

  # Pairs a and b at a * n + b: their protons' summed distance, matched the
  # closer way round.
  a1, a2, b1, b2 = ia * 2, ia * 2 + 1, ib * 2, ib * 2 + 1
  gap = apart(a1, b1) + apart(a2, b2)
  swap = apart(a1, b2) + apart(a2, b1)
  closer = swap < gap
  gap.set_selected(closer, swap.select(closer))
  kept = [0]
  for a in range(1, n):
    repeat = False
    for b in kept:
      if gap[a * n + b] < 1e-6:
        repeat = True
        break
    if not repeat:
      kept.append(a)
  sel = flex.size_t()
  for a in kept:
    sel.append(2 * a)
    sel.append(2 * a + 1)
  return pts.select(sel)


def _search_exhaustive(scores, allowed):
  """The allowed combination of candidates with the lowest total score.

  ``scores`` holds per water the score of each candidate, candidate 0 being
  its current orientation. ``allowed`` maps a pair of waters ``(i, j)``,
  ``i < j``, to whether each combination of their candidates may coexist,
  flattened row-major: candidate ``a`` of ``i`` with candidate ``b`` of ``j``
  at ``a * n_j + b``. Every combination is enumerated, the last water
  varying fastest, so combination 0 keeps every current orientation and
  wins a tie. Returns the chosen candidate per water.
  """
  sizes = [s.size() for s in scores]
  n = 1
  for size in sizes:
    n *= size
  combo = flex.size_t_range(n)
  pick = [None] * len(sizes)
  stride = 1
  for i in reversed(range(len(sizes))):
    pick[i] = (combo / stride) % sizes[i]
    stride *= sizes[i]
  total = flex.double(n, 0.0)
  for i, s in enumerate(scores):
    total += s.select(pick[i])
  ok = flex.bool(n, True)
  for (i, j), table in allowed.items():
    ok &= table.select(pick[i] * sizes[j] + pick[j])
  rows = ok.iselection()
  k = rows[flex.min_index(total.select(rows))]
  return [int(p[k]) for p in pick]


def _search_icm(scores, allowed):
  """A locally best allowed combination, by iterated conditional modes.

  Starting from every water's current orientation, each water in turn takes
  its lowest-scoring candidate allowed with the others' current choices,
  until a round changes nothing. Every change lowers the total score, so
  the rounds end. ``scores`` and ``allowed`` as for
  :func:`_search_exhaustive`. Returns the chosen candidate per water.
  """
  sizes = [s.size() for s in scores]
  # Per water, its pair tables as (other water, table, stride of this
  # water's candidate, stride of the other's).
  pairs = [[] for _ in sizes]
  for (i, j), table in allowed.items():
    pairs[i].append((j, table, sizes[j], 1))
    pairs[j].append((i, table, 1, sizes[j]))
  pick = [0] * len(sizes)
  changed = True
  while changed:
    changed = False
    for i, s in enumerate(scores):
      best = pick[i]
      for a in range(sizes[i]):
        if not s[a] < s[best]:
          continue
        fits = True
        for j, table, own_stride, other_stride in pairs[i]:
          if not table[a * own_stride + pick[j] * other_stride]:
            fits = False
            break
        if fits:
          best = a
      if best != pick[i]:
        pick[i] = best
        changed = True
  return pick


class _Clearance:
  """Result of one clearance test (see :meth:`_WaterHydrogenPlacer._clear`).

  ``ok`` is the per-candidate pass/fail flag; ``blocks`` and ``mins`` are the
  squared-distance blocks it came from (raw, and already reduced to
  per-candidate minima) that :meth:`best_at` reduces the clearance distance
  out of on demand.
  """

  __slots__ = ("ok", "blocks", "mins", "k")

  def __init__(self, k):
    self.ok = flex.bool(k, True)
    self.blocks = []   # (flattened (k, m) squared-distance block, m)
    self.mins = []     # per-candidate minima of already-reduced blocks
    self.k = k

  def best_at(self, rows):
    """Clearance distance of the candidates ``rows``.

    Distance to the nearest non-own atom within ``_WATER_CLEARANCE_RADIUS``, or
    the search radius if nothing is near.
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


class _WaterHydrogenPlacer:
  """Engine behind :func:`place_water_hydrogens`.

  Holds the placement state (static-atom KDTree, donor/acceptor bookkeeping,
  per-water neighbour blocks, placed protons) and the placement, refinement
  and basin-hopping logic.

  Only the placed water H move; each water's static neighbour block and
  water-neighbour list are built once in :meth:`_build_neighbourhoods`.

  The constructor records the parameters; :meth:`run` does the work, in
  place. Parameter semantics: :func:`place_water_hydrogens`.
  """

  def __init__(self, hier, oh_length, element, n_refine, refine_tol,
               n_basin, cluster_search, existing_h, on_state,
               crystal_symmetry, geometry):
    self.hier = hier
    self.oh_length = oh_length
    self.geometry = geometry
    self.basin_rounds_run = 0
    self.element = element
    self.n_refine = n_refine
    self.refine_tol = refine_tol
    self.n_basin = n_basin
    self.cluster_search = cluster_search
    self.existing_h = existing_h
    self.on_state = on_state
    self.crystal_symmetry = crystal_symmetry

    # Placed-H coordinates, one slot per proton, filled in _greedy_pass();
    # placed_xyz is the same data as a flex.vec3_double.
    self.placed_coords = []
    self.placed_xyz = None
    # Sweep active-set stamps (see _dirty), sized in _build_neighbourhoods():
    # off one monotonic counter, the tick at which each water was last placed
    # and the tick at which its H last moved.
    self.w_placed_at = []
    self.w_moved_at = []
    self.tick = 0
    # (water index, [(atom, slot, di), ...], fixed_d1)
    self.records = []
    # (residue id, (metal element, distance) or None, action)
    self.partial_waters = []

  def _static_block(self, wi, cands):
    """Candidates against water ``wi``'s static neighbours.

    Returns ``(d, m, ok)``: the flattened (k, m) squared-distance block, its
    neighbour count, and whether each candidate clears every one of them. ``m``
    is zero (and ``d`` None) when the water has no static neighbour at all.
    """
    k = cands.size()
    P = self.w_sxyz[wi]
    m = P.size()
    ok = flex.bool(k, True)
    if not m:
      return None, 0, ok
    ik, im = _tile(k, m)
    # Squared distances via parts() and pow2, never .dot() or .norms() (fused
    # multiply-add changes the last bits); likewise throughout this file.
    dx, dy, dz = (cands.select(ik) - P.select(im)).parts()
    d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
    bad = (d < self.w_sthr[wi].select(im)).iselection()
    if bad.size():
      ok.set_selected(bad / m, False)
    return d, m, ok

  def _static_clear(self, wi, cands, key):
    """Static half of the clearance test, reduced and cached per water.

    ``cands`` must be an invariant candidate set (the acceptor directions,
    the fallback sphere, or one acceptor pair's even-split protons) and
    ``key`` the geometry-dict slot to cache under.
    Returns ``(ok, mins)``: the per-candidate flag and the per-candidate
    minimum squared distance, the latter None when the water has no static
    neighbour.
    """
    g = self.w_geom[wi]
    st = g.get(key)
    if st is None:
      d, m, ok = self._static_block(wi, cands)
      k = cands.size()
      st = g[key] = (ok, _row_mins(d, k, m, flex.size_t_range(k)) if m
                     else None)
    return st

  def _own_sym_equiv_clear(self, wi, cands, own_fixed=()):
    """Candidates against the equivalents of water ``wi``'s own protons.

    A water whose own equivalent is in reach cannot be scored against a
    standing point set: moving a candidate moves its equivalent with it. Each
    operator is applied to the whole candidate array instead, which gives
    every candidate's distance to its own equivalent, and to the equivalent of
    each proton in ``own_fixed``, this water's protons already settled this
    pass.

    Returns ``(ok, mins)``, both None when no operator brings this water's own
    equivalent into range.
    """
    ops = self.w_self_ops[wi]
    if not ops:
      return None, None
    k = cands.size()
    ok = flex.bool(k, True)
    best = None
    thr = _WATER_MIN_H_CLEARANCE ** 2
    for rot, trn in ops:
      equiv = rot * cands + trn
      pairs = [(cands, equiv)]
      for p in own_fixed:
        # The candidate against a fixed proton's equivalent, which stands
        # still, and the fixed proton against the candidate's, which does not.
        pairs.append((cands, (rot * flex.vec3_double([p]) + trn)[0]))
        pairs.append((equiv, p))
      for left, right in pairs:
        dx, dy, dz = (left - right).parts()
        d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
        bad = (d < thr).iselection()
        if bad.size():
          ok.set_selected(bad, False)
        if best is None:
          best = d
        else:
          closer = d < best
          best.set_selected(closer, d.select(closer))
    return ok, best

  def _clear(self, wi, cands, nbr_slots, static_key=None, own_fixed=()):
    """Clearance test over the candidate H positions of one water.

    One dense (candidate x neighbour) block against water ``wi``'s static
    block (a ball of ``oh_length + _WATER_CLEARANCE_RADIUS`` about the O, own
    atoms dropped) and its neighbouring protons.

    ``nbr_slots`` are the pool slots that may lie near this water, its own
    excluded; they index symmetry equivalents as readily as placed protons.
    ``static_key`` names the geometry-dict slot holding the cached
    static half, for a candidate set that does not move; None for the cone,
    which is rebuilt around each pass's O-H1 axis. ``own_fixed`` are this
    water's protons already settled this pass, which only matter to a water
    whose own equivalent is in reach (see :meth:`_own_sym_equiv_clear`).

    Returns a :class:`_Clearance` whose ``ok`` is True where the candidate
    clears every heavy atom by ``_WATER_MIN_CLEARANCE`` and every hydrogen
    (static or placed) by ``_WATER_MIN_H_CLEARANCE``; the clearance distances
    are left to :meth:`_Clearance.best_at`.
    """
    # Squared distances throughout.
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
      dx, dy, dz = (ct - self.pool_xyz.select(nbr_slots).select(im)).parts()
      d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
      bad = (d < _WATER_MIN_H_CLEARANCE ** 2).iselection()
      if bad.size():
        cl.ok.set_selected(bad / m, False)
      cl.blocks.append((d, m))
    s_ok, s_mins = self._own_sym_equiv_clear(wi, cands, own_fixed)
    if s_ok is not None:
      cl.ok &= s_ok
      cl.mins.append(s_mins)
    return cl

  def _stats(self):
    """:func:`_water_clash_stats` from the engine's own arrays."""
    ns = len(self.placed_coords)
    if ns:
      X = self.placed_xyz[:ns].concatenate(self.wh_xyz)
      W = self.slot_wid[:ns].concatenate(self.wh_wid)
      C = self.slot_conf[:ns].concatenate(self.wh_conf)
    else:
      X, W, C = self.wh_xyz, self.wh_wid, self.wh_conf
    anchor = self.slot_anchor[:ns] + self.wh_anchor
    # Without water H in an altloc every pair coexists.
    if (C == 0).all_eq(True):
      C = None
    d = _water_h_contacts(X, W, C)[2]
    if self.crystal_symmetry is not None:
      # The slots are final once the greedy pass is done, so one pairing
      # serves every state.
      if self.h_sym_pairs is None:
        self.h_sym_pairs = _water_h_sym_pairs(
          X, W, anchor, self.crystal_symmetry, C)
      for _op, _i, _j, d_sym in _water_h_sym_contacts(
          X, self.h_sym_pairs, self.crystal_symmetry.unit_cell()):
        d = d.concatenate(d_sym)
    return _contact_stats(X.size(), d)

  def _nearest_cation(self, o_xyz, own_idx, c):
    """Closest metal cation coordinating a water O, if any.

    Returns ``(element, distance)`` or None, within
    ``_WATER_METAL_COORD_RADIUS`` (a first-shell bond), excluding the
    water's own atoms ``own_idx`` and cations of an altloc other than the
    water's conformer index ``c``. Reporting only; placement ignores it.
    """
    o = matrix.col(o_xyz)
    best = None
    for i in self.static_tree.query_ball_point(tuple(o_xyz),
                                               _WATER_METAL_COORD_RADIUS):
      if i in own_idx:
        continue
      if c and not _compatible(c, self.conf[i]):
        continue
      el = self.atoms[i].element.strip().upper()
      if el not in _WATER_CATION_ELEMENTS:
        continue
      d = (matrix.col(self.atoms[i].xyz) - o).length()
      if best is None or d < best[1]:
        best = (el, d)
    return best

  def _site_symmetry(self, o_pts):
    """Site symmetry of the water O ``o_pts``, for the waters on a special
    position.

    Fills ``w_site_ops`` with the cartesian (rotation, translation) of each
    such water's non-identity site-symmetry operators and, on a two-fold,
    ``w_axis`` with the unit axis direction. Under these operators the
    water's own atoms map onto the water itself (a symmetric water's other
    H) or onto its other orientation, never onto a neighbour.
    """
    sps = crystal.special_position_settings(
      self.crystal_symmetry, min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL)
    uc = self.crystal_symmetry.unit_cell()
    table = sps.site_symmetry_table(sites_cart=o_pts)
    for wi in table.special_position_indices():
      ops = [m for m in table.get(wi).matrices() if not m.is_unit_mx()]
      self.w_site_ops[wi] = tuple(
        (matrix.sqr(uc.matrix_cart(m.r())),
         matrix.col(uc.orthogonalize(m.t().as_double()))) for m in ops)
      info = ops[0].r().info()
      if len(ops) == 1 and info.type() == 2:
        axis = matrix.col(uc.orthogonalize(info.ev())).normalize()
        self.w_axis[wi] = axis.elems

  def _own_equiv(self, wi, j):
    """Whether static atom ``j`` is the equivalent of one of water ``wi``'s
    own atoms under the water's site symmetry."""
    if j < self.n_model_atoms:
      return False
    src = self.sym_src[j - self.n_model_atoms]
    if src not in self.w_own[wi]:
      return False
    x = matrix.col(self.static_xyz[src])
    y = matrix.col(self.static_xyz[j])
    for rot, trn in self.w_site_ops[wi]:
      if (rot * x + trn - y).length() < 0.01:
        return True
    return False

  def _planar_n(self, i):
    """Whether static atom ``i`` is an N bonded to three heavy atoms in a
    plane, as a proline, amide or nucleoside N or one bound to a metal,
    whose lone pair is not free to accept an H-bond.

    A bond is a heavy atom within ``_WATER_N_BOND``, or a metal within
    ``_WATER_METAL_COORD_RADIUS``, in the N's conformer or blank, and the
    three bonds of one conformer are planar when their angles sum to more
    than ``_WATER_PLANAR_N``. Each N is decided once.
    """
    planar = self.planar_n.get(i)
    if planar is not None:
      return planar
    planar = False
    if self.atoms[i].element.strip().upper() == "N":
      x = matrix.col(self.static_xyz[i])
      ci = self.conf[i]
      bonds = []
      for j in self.static_tree.query_ball_point(self.static_xyz[i],
                                                 _WATER_METAL_COORD_RADIUS):
        if j == i or self.static_is_h[j] or not _compatible(ci, self.conf[j]):
          continue
        v = matrix.col(self.static_xyz[j]) - x
        d = v.length()
        if self.atoms[j].element.strip().upper() in _WATER_CATION_ELEMENTS:
          limit = _WATER_METAL_COORD_RADIUS
        else:
          limit = _WATER_N_BOND
        if 1e-3 < d < limit:
          bonds.append((self.conf[j], v / d))
      confs = {c for c, _ in bonds if c}
      if not confs:
        confs = {0}
      for c in confs:
        u = [v for cj, v in bonds if not cj or cj == c]
        if len(u) != 3:
          continue
        angles = (u[0].angle(u[1], deg=True) + u[0].angle(u[2], deg=True)
                  + u[1].angle(u[2], deg=True))
        if angles > _WATER_PLANAR_N:
          planar = True
    self.planar_n[i] = planar
    return planar

  def _water_geom(self, wi):
    """Fixed per-water geometry, computed once.

    Returns the dict of ``o`` (O coordinates), ``acceptors`` (atom indices,
    nearest first), ``acc_water`` (whether each is a water atom),
    ``acc_dirs``/``acc_pts`` (unit directions from the O toward each
    acceptor, as tuples, and the H positions built from them), ``cat``
    (O->cation unit directions), ``sph`` (the lazily built fallback-sphere H
    positions) and ``acc_static``/``sph_static`` (the cached static half of
    the clearance test for those two candidate sets, see
    :meth:`_static_clear`). :meth:`_place_one` adds, per acceptor pair
    ``(k1, k2)``, the even split's ``("split", k1, k2)`` (its two H positions
    and whether both are away from cations, None when the acceptors are
    collinear) and ``("split_static", k1, k2)`` (the cached static half of
    its clearance test).
    """
    g = self.w_geom[wi]
    if g is not None:
      return g
    atoms = self.atoms
    o_xyz = self.w_o_col[wi]
    own_idx = self.w_own[wi]
    c = self.w_conf[wi]
    # Acceptors (static O/N/F/S/Cl) within range, nearest first.
    acceptors = []
    for i in self.w_acc_raw[wi]:
      if i in own_idx:
        continue
      if any(_compatible(c, ch) for ch in self.donor_n.get(i, ())):
        continue  # protonated N: a donor, not a usable acceptor
      if self._planar_n(i):
        continue  # no free lone pair
      a = atoms[i]
      if a.element.strip().upper() not in _WATER_ACCEPTOR_ELEMENTS:
        continue
      d = (matrix.col(a.xyz) - o_xyz).length()
      if d < 1e-3:
        continue  # atom coincident with O (alt-conf / overlap): no direction
      acceptors.append((d, i))
    acceptors.sort(key=lambda t: t[0])
    acceptors = [i for _, i in acceptors]

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

    o = o_xyz.elems
    acc_dirs = [(matrix.col(atoms[i].xyz) - o_xyz).normalize().elems
                for i in acceptors]
    g = {"o": o,
         "acceptors": acceptors,
         "acc_water": [_is_water(atoms[i].parent().resname) for i in acceptors],
         "acc_dirs": acc_dirs,
         "acc_pts": flex.vec3_double(acc_dirs) * self.oh_length + o,
         "cat": [c.elems for c in cation_dirs],
         "sph": None,
         "acc_static": None,
         "sph_static": None}
    self.w_geom[wi] = g
    return g

  def _place_symmetric(self, wi, nbr_slots, fixed_d1):
    """The two H of a water on a two-fold, each other's equivalents.

    H2 is H1 turned 180 deg about the axis through the O, so both lie on one
    cone of half the H-O-H angle about either end of the axis. A given H1
    fixes H2; otherwise the pair is chosen over sampled azimuths, plus the
    azimuth of each acceptor, ranked (both away from cations, both
    clash-free, summed alignment of each H with its best acceptor, the
    smaller of the two clearances).

    Returns ``(h1_xyz, h2_xyz)``.
    """
    g = self._water_geom(wi)
    o = g["o"]
    a = self.w_axis[wi]
    axis = matrix.col(a)
    turn = axis.axis_and_angle_as_r3_rotation_matrix(180, deg=True)
    oh = self.oh_length

    def at(d):
      return (o[0] + oh * d[0], o[1] + oh * d[1], o[2] + oh * d[2])

    if fixed_d1 is not None:
      return at(fixed_d1), at((turn * matrix.col(fixed_d1)).elems)
    p = axis.ortho().normalize()
    p, q = p.elems, axis.cross(p).elems
    ch, sh = self.cos_half_hoh, self.sin_half_hoh
    acc_dirs = g["acc_dirs"]
    psis = [math.pi * k / _WATER_CONE_SAMPLES
            for k in range(2 * _WATER_CONE_SAMPLES)]
    for u in acc_dirs:
      psis.append(math.atan2(u[0] * q[0] + u[1] * q[1] + u[2] * q[2],
                             u[0] * p[0] + u[1] * p[1] + u[2] * p[2]))
    d1s = []
    for s in (1.0, -1.0):
      for psi in psis:
        cp, sp = math.cos(psi), math.sin(psi)
        d1s.append(tuple(s * ch * a[c] + sh * (cp * p[c] + sp * q[c])
                         for c in range(3)))
    d2s = [(turn * matrix.col(d)).elems for d in d1s]
    pts1 = flex.vec3_double([at(d) for d in d1s])
    pts2 = flex.vec3_double([at(d) for d in d2s])
    cl1 = self._clear(wi, pts1, nbr_slots)
    cl2 = self._clear(wi, pts2, nbr_slots)
    cat = g["cat"]
    top, rank = _rank_top(_cat_ok(cat, pts1, o) & _cat_ok(cat, pts2, o),
                          cl1.ok & cl2.ok)
    if rank == 3 and acc_dirs:
      align = flex.double()
      for d1, d2 in zip(d1s, d2s):
        align.append(
          max(u[0] * d1[0] + u[1] * d1[1] + u[2] * d1[2] for u in acc_dirs)
          + max(u[0] * d2[0] + u[1] * d2[1] + u[2] * d2[2] for u in acc_dirs))
      top = top & (align == flex.max(align.select(top)))
    rows = top.iselection()
    clear = cl1.best_at(rows)
    c2 = cl2.best_at(rows)
    lower = c2 < clear
    clear.set_selected(lower, c2.select(lower))
    j = int(rows[flex.max_index(clear)])
    return pts1[j], pts2[j]

  def _place_one(self, wi, nbr_slots, fixed_d1=None):
    """The two H positions for one water O, clash-aware.

    ``nbr_slots`` are the pool slots that may lie near this water, its own
    excluded. ``fixed_d1`` holds the unit O-H1 direction fixed instead of
    searching for one, for a water that already carries a proton: the returned
    ``h1_xyz`` then only restates it and just ``h2_xyz`` is new. A water on a
    two-fold goes to :meth:`_place_symmetric`.

    Returns ``(h1_xyz, h2_xyz)``.
    """
    if self.w_axis[wi] is not None:
      return self._place_symmetric(wi, nbr_slots, fixed_d1)
    g = self._water_geom(wi)
    o = g["o"]
    acceptors = g["acceptors"]
    acc_dirs = g["acc_dirs"]
    acc_pts = g["acc_pts"]
    cat = g["cat"]
    na = len(acceptors)

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
        a = acc_dirs[k]
        align = a[0] * d1[0] + a[1] * d1[1] + a[2] * d1[2]
        if align > best_align:
          best_align = align
          h1_k = k
    else:
      if na:
        acc_cl = self._clear(wi, acc_pts, nbr_slots, "acc_static")
        acc_ok = acc_cl.ok
        acc_cat = _cat_ok(cat, acc_pts, o)
      d1 = None
      h1_k = -1
      for k in range(na):
        if acc_ok[k] and acc_cat[k]:  # clash-free and away from cations
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
          # The candidates are the acceptors, then the fallback sphere; each
          # half reads its clearances back from its own test.
          lo = rows.select(rows < na)
          hi = rows.select(~(rows < na)) - na
          vals = acc_cl.best_at(lo).concatenate(s_cl.best_at(hi))
          b = flex.max_index(vals)
          pick = int(lo[b]) if b < lo.size() else na + int(hi[b - lo.size()])
        else:
          pick = int(rows[flex.max_index(s_cl.best_at(rows))])
        if pick < na:
          d1 = acc_dirs[pick]
          h1_k = pick
        else:
          d1 = _FALLBACK_DIRS[pick - na]
    h1_xyz = (o[0] + self.oh_length * d1[0],
              o[1] + self.oh_length * d1[1],
              o[2] + self.oh_length * d1[2])

    v = matrix.col(d1)
    p = v.ortho().normalize()
    p, q = p.elems, v.cross(p).elems

    # H2: rank the cone angles by (away-from-cation, clash-free, best
    # alignment over the acceptors H1 did not take, clearance).
    ns = _WATER_CONE_SAMPLES
    cone_dirs = ((flex.vec3_double(ns, p) * _CONE_COS
                  + flex.vec3_double(ns, q) * _CONE_SIN) * self.sin_hoh
                 + (self.cos_hoh * d1[0], self.cos_hoh * d1[1],
                    self.cos_hoh * d1[2]))
    cone_pts = cone_dirs * self.oh_length + o
    c_cl = self._clear(wi, cone_pts, nbr_slots, own_fixed=(h1_xyz,))
    c_cat = _cat_ok(cat, cone_pts, o)

    top, rank = _rank_top(c_cat, c_cl.ok)
    # Alignment ranks only cation-ok, clash-free candidates (rank 3).
    use = [k for k in range(na) if k != h1_k]
    if rank == 3 and use:
      dx, dy, dz = cone_dirs.parts()
      dots = []
      for k in use:
        a = acc_dirs[k]
        dots.append(dx * a[0] + dy * a[1] + dz * a[2])
      # Of the acceptors H2 can reach (within the donation limit at a
      # cation-ok, clash-free cone angle), only those outside the waters
      # count when it can reach both kinds.
      reach = [i for i in range(len(use))
               if flex.max(dots[i].select(top)) >= _WATER_H1_ACC_ALIGN]
      keep = [i for i in reach if not g["acc_water"][use[i]]]
      if keep and len(keep) < len(reach):
        use = [use[i] for i in keep]
        dots = [dots[i] for i in keep]
      align = dots[0]
      for dot in dots[1:]:
        s = dot > align
        align.set_selected(s, dot.select(s))
      top = top & (align == flex.max(align.select(top)))
    rows = top.iselection()
    j = int(rows[flex.max_index(c_cl.best_at(rows))])
    h2_xyz = cone_pts[j]

    # Both protons on acceptors: turn the water so each is equally far off its
    # acceptor, if that stays clash-free and away from cations. Not for a
    # given H1.
    if fixed_d1 is not None or h1_k < 0 or rank != 3 or not use:
      return h1_xyz, h2_xyz
    d2 = cone_dirs[j]
    k2, best = -1, _WATER_H1_ACC_ALIGN
    for k in use:
      a = acc_dirs[k]
      dot = d2[0] * a[0] + d2[1] * a[1] + d2[2] * a[2]
      if dot >= best:
        k2, best = k, dot
    if k2 < 0:
      return h1_xyz, h2_xyz
    # Positions, cation test and static clearance, cached per water and
    # acceptor pair.
    key = (h1_k, k2)
    split = g.get(("split",) + key, False)
    if split is False:
      split = None
      dirs = _split_dirs(acc_dirs[h1_k], acc_dirs[k2],
                         self.cos_half_hoh, self.sin_half_hoh)
      if dirs is not None:
        pts = flex.vec3_double()
        for d in dirs:
          pts.append(tuple(o[c] + self.oh_length * d[c] for c in range(3)))
        split = (pts, _cat_ok(cat, pts, o).all_eq(True))
      g[("split",) + key] = split
    if split is None or not split[1]:
      return h1_xyz, h2_xyz
    pts = split[0]
    s1, s2 = pts[0], pts[1]
    if self._clear(wi, pts, nbr_slots, ("split_static",) + key,
                   own_fixed=(s1,)).ok.all_eq(True):
      return s1, s2
    return h1_xyz, h2_xyz

  def _pool_write(self, wi, slot, xyz):
    """Put one proton of water ``wi`` at ``xyz``, its equivalents with it."""
    self.placed_xyz[slot] = xyz
    if self.pool_xyz is self.placed_xyz:
      return
    self.pool_xyz[slot] = xyz
    v = matrix.col(xyz)
    for b, rot, trn in self.w_blocks[wi]:
      self.pool_xyz[b * self.n_slots + slot] = (rot * v + trn).elems

  def _store(self, wi, slots, h1, h2):
    """Write one water's new H positions to the model and the arrays.

    True if any of them moved (see :meth:`_dirty`).
    """
    moved = False
    for atom, slot, di in slots:
      xyz = h1 if di == 1 else h2
      if xyz != self.placed_coords[slot]:
        moved = True
      self.placed_coords[slot] = xyz
      self._pool_write(wi, slot, xyz)
      atom.set_xyz(xyz)
    return moved

  def _dirty(self, wi):
    """Whether re-placing water ``wi`` could move its H.

    True if the water was kicked, or if a placed H of the waters in
    ``w_wnbr[wi]`` or their symmetry equivalents in ``w_sym_nbr[wi]`` moved
    since it was last placed. The tick is bumped once per water per sweep.
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
    for _b, wj in self.w_sym_nbr[wi]:
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
    """Capture the placed-H state as ``(coords, placed_at, moved_at)``."""
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
      own = [slot for _, slot, _ in slots]
      pts = self.placed_xyz.select(flex.size_t(own))
      # tuple() of a vec3_double would iterate it through IndexError.
      own_fixed = tuple(self.placed_coords[slot] for slot in own)
      if not self._clear(wi, pts, self.w_nbr_slots[wi],
                         own_fixed=own_fixed).ok.all_eq(True):
        bad.append(ri)
    return bad

  def _kick(self, ri, rng):
    """Re-orient the water ``records[ri]`` randomly, from seeded ``rng``.

    H1 goes along a random axis and H2 to a random azimuth on the H-O-H cone.
    A water being completed has a fixed O-H1, so only the azimuth is random.
    On a two-fold, H1 goes to a random azimuth on the cone about either end
    of the axis (unless given) and H2 to its equivalent, so a completed water
    there stays put.

    Returns True if any of the water's H moved.
    """
    wi, slots, fixed_d1 = self.records[ri]
    o = self._water_geom(wi)["o"]
    a = self.w_axis[wi]
    if a is None:
      d1 = fixed_d1
      if d1 is None:
        d1 = rng.random_double_point_on_sphere()
      v = matrix.col(d1)
      p = v.ortho().normalize()
      p, q = p.elems, v.cross(p).elems
      theta = 2.0 * math.pi * rng.random_double()
      ct, st = math.cos(theta), math.sin(theta)
      d2 = tuple(self.cos_hoh * d1[c] + self.sin_hoh * (ct * p[c] + st * q[c])
                 for c in range(3))
    else:
      axis = matrix.col(a)
      d1 = fixed_d1
      if d1 is None:
        p = axis.ortho().normalize()
        p, q = p.elems, axis.cross(p).elems
        s = 1.0 if rng.random_double() < 0.5 else -1.0
        theta = 2.0 * math.pi * rng.random_double()
        ct, st = math.cos(theta), math.sin(theta)
        d1 = tuple(s * self.cos_half_hoh * a[c]
                   + self.sin_half_hoh * (ct * p[c] + st * q[c])
                   for c in range(3))
      turn = axis.axis_and_angle_as_r3_rotation_matrix(180, deg=True)
      d2 = (turn * matrix.col(d1)).elems
    self.tick += 1
    moved = self._store(wi, slots,
                        tuple(o[c] + self.oh_length * d1[c] for c in range(3)),
                        tuple(o[c] + self.oh_length * d2[c] for c in range(3)))
    if moved:
      self.w_moved_at[wi] = self.tick
    return moved

  def run(self):
    """Place H on every bare water, refine, optionally basin-hop, keep best.

    Modifies the hierarchy in place. Returns the kept-state label (see
    :func:`place_water_hydrogens`).
    """
    hier = self.hier
    assert hier.models_size() <= 1, (
      f"place_water_hydrogens takes one model, not {hier.models_size()}")
    check_for_missing_elements(hier)
    # Resolve before stripping, which would remove the D this keys on.
    if self.geometry is None:
      self.geometry = "neutron" if _has_deuterium(hier) else "xray"
    oh, self.hoh_deg = _WATER_GEOMETRY[self.geometry]
    if self.oh_length is None:
      self.oh_length = oh
    # The strip records what each water carried.
    stripped = {}
    if self.existing_h == "reorient":
      stripped = _strip_water_hydrogens(hier)

    sel = hier.atoms()
    atoms = list(sel)
    if not atoms:
      return None
    sel.reset_i_seq()
    # Conformer index per atom, 0 for blank (see _compatible).
    ci = hier.get_conformer_indices()
    if " " in ci.index_altloc_mapping:
      hier.overall_counts().raise_improper_alt_conf_if_necessary()
    conf = ci.conformer_indices
    altloc_index = dict(ci.index_altloc_mapping)
    waters, single = self._collect_waters(conf, altloc_index, stripped)
    self._build_environment(sel, atoms, conf, waters, single)
    waters = self._build_neighbourhoods(waters)
    self._greedy_pass(waters, stripped)
    kept_label = self._refine_and_hop()
    if self.cluster_search:
      self._search_clusters()
    return kept_label

  def _collect_waters(self, conf, altloc_index, stripped):
    """Collect the waters to place and the single-H waters to report."""
    hier = self.hier
    # The waters to protonate, one per water conformer (see
    # _water_conformers), as (target atom group, O, own atom indices,
    # fixed_d1, existing H, atom names, water id, conformer index, occupancy,
    # elements of the H a reorient stripped); every water H is also listed
    # once for the clash statistics.
    waters = []
    # Single-H waters to report, once per such H, annotated once the tree
    # exists.
    single = {}
    wh_xyz = []
    wh_wid = []
    wh_conf = []
    wh_anchor = []
    # cctbx counts an O at occupancy 1 on a site of order n once (a lower
    # occupancy already carries the multiplicity): each of its H accounts
    # for n times its occupancy, and new H take 1/n.
    sps = None
    if self.crystal_symmetry is not None:
      sps = crystal.special_position_settings(
        self.crystal_symmetry, min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL)
    for wid, (rg, ags) in enumerate(_water_residue_groups(hier)):
      o_sites = _water_o_sites(ags)
      blank = None
      for ag in ags:
        if not ag.altloc:
          blank = ag
        ats, hd = _hd_flags(ag)
        for k in range(ats.size()):
          if hd[k]:
            a = ats[k]
            wh_xyz.append(a.xyz)
            wh_wid.append(wid)
            wh_conf.append(conf[a.i_seq])
            wh_anchor.append(_h_anchor(a.xyz, ag.altloc, o_sites))
      for altloc, ats, target in _water_conformers(rg, ags):
        o = None
        existing = []
        names = set()
        own_idx = set()
        els = ats.extract_element(strip=True)
        for k in range(ats.size()):
          a = ats[k]
          names.add(a.name.strip())
          own_idx.add(a.i_seq)
          if els[k] in ("H", "D"):
            existing.append(a)
          elif o is None and els[k].upper() == "O":
            o = a
        if o is None:
          continue
        # What a reorient stripped from this conformer: its blank H and its
        # own altloc's, each with the atom group that held it.
        carried = []
        for g in ([target] if target is blank else [blank, target]):
          if g is None:
            continue
          for el, occ in stripped.get(g.memory_id(), ()):
            carried.append((el, occ, g))
        n_site = 1
        if sps is not None and o.occ == 1:
          n_site = sps.site_symmetry(site_cart=o.xyz).n_matrices()
        if n_site > 1:
          n_h = sum(int(n_site * a.occ + 0.5) for a in existing)
          n_carried = sum(int(n_site * q + 0.5) for _, q, _ in carried)
        else:
          n_h = len(existing)
          n_carried = len(carried)
        fixed_d1 = None
        skip = n_h >= 2                    # already protonated
        if existing and not skip:
          if self.existing_h != "complete":
            skip = True
          else:
            v = matrix.col(existing[0].xyz) - matrix.col(o.xyz)
            if v.length() < 1e-3:
              skip = True  # H coincident with O: no direction for a cone
            else:
              fixed_d1 = v.normalize().elems
        is_single = (n_carried == 1 if self.existing_h == "reorient"
                     else n_h == 1)
        if skip and not is_single:
          continue                         # nothing to place and nothing to say
        c = 0
        if altloc:
          # An altloc the strip emptied gets an index of its own.
          if altloc not in altloc_index:
            altloc_index[altloc] = max(altloc_index.values()) + 1
          c = altloc_index[altloc]
        if is_single:
          # A blank H shared by two conformers is one water to report.
          holder = (carried[0][2] if self.existing_h == "reorient"
                    else existing[0].parent())
          if self.existing_h == "reorient":
            action = "stripped"
          elif fixed_d1 is not None:
            action = "completed"
          else:
            action = "kept"
          key = (rg.memory_id(), holder.memory_id())
          if key not in single:
            single[key] = (_water_id(holder), o.xyz, own_idx, c, action)
        if skip:
          continue
        # New H take the conformer's occupancy: its own altloc's atoms', or
        # what the strip took from them, else the O's; divided by the site
        # order of an O at occupancy 1.
        occ = o.occ
        if altloc:
          tats = target.atoms()
          if tats.size():
            occ = tats[0].occ
          else:
            for _, q, g in carried:
              if g is target:
                occ = q
                break
        occ /= n_site
        waters.append((target, o, own_idx, fixed_d1, existing, names, wid, c,
                       occ, [el for el, _, _ in carried]))
    single = list(single.values())
    self.wh_xyz = flex.vec3_double(wh_xyz)
    self.wh_wid = flex.size_t(wh_wid)
    self.wh_conf = flex.size_t(wh_conf)
    self.wh_anchor = wh_anchor
    # Water H pairs across symmetry for the clash counts, built by _stats.
    self.h_sym_pairs = None
    return waters, single

  def _build_environment(self, sel, atoms, conf, waters, single):
    """Build the static environment and annotate the single-H waters."""
    from scipy.spatial import KDTree
    hier = self.hier
    # Symmetry-equivalent environment atoms, appended after the model's own so
    # the indices 0..n-1 stay valid.
    self.sym_hier = []
    sym_atoms = []
    sym_xyz = flex.vec3_double()
    sym_src = flex.size_t()
    if self.crystal_symmetry is not None and (waters or single):
      self.sym_hier, sym_atoms, sym_xyz, sym_src = _symmetry_environment(
        hier, self.crystal_symmetry,
        max(self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
            _WATER_ACCEPTOR_RADIUS))
    self.n_model_atoms = len(atoms)
    self.sym_src = sym_src
    atoms = atoms + sym_atoms
    self.atoms = atoms
    # A copy keeps its source's altloc.
    self.conf = conf.concatenate(conf.select(sym_src))

    # Static neighbours (protein, ligands, water O, pre-existing H) never
    # move; the placed water H are tracked by slot in placed_coords/placed_xyz.
    self.static_xyz = sel.extract_xyz()
    self.static_xyz.extend(sym_xyz)
    self.static_tree = KDTree(self.static_xyz.as_numpy_array())
    _el = sel.extract_element(strip=True)
    _el.extend(flex.std_string([a.element.strip() for a in sym_atoms]))
    self.static_is_h = (_el == "H") | (_el == "D")
    self.static_thr = flex.double(len(atoms), _WATER_MIN_CLEARANCE ** 2)
    self.static_thr.set_selected(self.static_is_h,
                                 _WATER_MIN_H_CLEARANCE ** 2)

    # N atoms carrying an H are donors, not acceptors, to the waters that
    # meet the H: per N, the conformer indices of its H (see _compatible).
    self.donor_n = {}
    h_idx = self.static_is_h.iselection()
    if h_idx.size():
      for i, nbrs in zip(h_idx, self.static_tree.query_ball_point(
          self.static_xyz.select(h_idx).as_numpy_array(), _WATER_NH_BOND,
          return_sorted=False)):
        for j in nbrs:
          if (j != i and atoms[j].element.strip().upper() == "N"
              and _compatible(self.conf[i], self.conf[j])):
            if j not in self.donor_n:
              self.donor_n[j] = set()
            self.donor_n[j].add(self.conf[i])
    # Per static N looked at so far, whether it is planar (see _planar_n).
    self.planar_n = {}

    self.cos_hoh = math.cos(math.radians(self.hoh_deg))
    self.sin_hoh = math.sin(math.radians(self.hoh_deg))
    self.cos_half_hoh = math.cos(math.radians(self.hoh_deg) / 2.0)
    self.sin_half_hoh = math.sin(math.radians(self.hoh_deg) / 2.0)

    # The cation coordinating each single-H water, symmetry copies included.
    self.partial_waters = [
      (rid, self._nearest_cation(o_xyz, own_idx, c), action)
      for rid, o_xyz, own_idx, c, action in single]

  def _build_neighbourhoods(self, waters):
    """Order the waters most-crowded first and build their neighbour lists."""
    from scipy.spatial import KDTree
    # Per-water neighbour lists, built once for every pass.
    n = len(waters)
    self.w_geom = [None] * n
    self.placed_coords = []
    self.placed_xyz = flex.vec3_double(2 * n, (0.0, 0.0, 0.0))
    # Protons the clearance test draws on: the placed ones, plus one
    # transformed block per symmetry operator that brings another water's
    # protons into reach, block b slot s at b * n_slots + s. Without symmetry
    # the pool is the placed array itself.
    self.n_slots = 2 * n
    self.pool_xyz = self.placed_xyz
    self.sym_tf = [None]
    self.w_sym_nbr = [[] for _ in range(n)]
    self.w_blocks = [()] * n
    # Cartesian (rotation, translation) per operator that brings the water's
    # own equivalent into range.
    self.w_self_ops = [()] * n
    # Per water on a special position (see _site_symmetry): its site
    # symmetry's operators and, on a two-fold, the axis.
    self.w_site_ops = [()] * n
    self.w_axis = [None] * n
    self.slot_wid = flex.size_t(2 * n, 0)
    self.slot_conf = flex.size_t(2 * n, 0)
    # The O site each placed H is bound to, for the clash statistics.
    self.slot_anchor = [None] * (2 * n)
    self.w_conf = []
    self.records = []   # (water index, [(atom, slot, di), ...], fixed_d1)
    self.w_wnbr = []
    # Every water is dirty for the first sweep.
    self.w_placed_at = [-1] * n
    self.w_moved_at = [0] * n
    self.tick = 0
    if n:
      o_pts = flex.vec3_double([w[1].xyz for w in waters])
      w_conf = [w[7] for w in waters]
      # Order most-crowded first. ``crowd`` is the neighbour count within
      # the clash radius, the water's own atoms and atoms of other altlocs
      # excluded.
      crowd = []
      for k, nb in enumerate(self.static_tree.query_ball_point(
          o_pts.as_numpy_array(), _WATER_CLEARANCE_RADIUS,
          return_sorted=False)):
        count = 0
        for j in nb:
          if j in waters[k][2]:
            continue
          if w_conf[k] and not _compatible(w_conf[k], self.conf[j]):
            continue
          count += 1
        crowd.append(count)
      order = sorted(range(n), key=lambda k: crowd[k], reverse=True)
      waters = [waters[k] for k in order]
      # Conformer index per water; the filters below run only for the
      # waters in an altloc, a blank water meeting every atom.
      self.w_conf = [w_conf[k] for k in order]
      o_pts = o_pts.select(flex.size_t(order))
      o_np = o_pts.as_numpy_array()
      self.w_own = [w[2] for w in waters]
      self.w_o_col = [matrix.col(w[1].xyz) for w in waters]
      if self.crystal_symmetry is not None:
        self._site_symmetry(o_pts)
      # query_ball_point sorts a batch's indices but not a single point's, and
      # the acceptor order is that sort's tie-break: the batch must not sort.
      self.w_acc_raw = self.static_tree.query_ball_point(
        o_np, _WATER_ACCEPTOR_RADIUS, return_sorted=False)
      self.w_cat_raw = self.static_tree.query_ball_point(
        o_np, _WATER_CATION_RADIUS, return_sorted=False)
      for wi, c in enumerate(self.w_conf):
        if c:
          self.w_acc_raw[wi] = [j for j in self.w_acc_raw[wi]
                                if _compatible(c, self.conf[j])]
          self.w_cat_raw[wi] = [j for j in self.w_cat_raw[wi]
                                if _compatible(c, self.conf[j])]
      # One static-neighbour block per water: every candidate H lies on the
      # O-H sphere about the O, so a single ball of oh_length + clearance
      # covers every candidate's own clearance ball.
      self.w_sxyz = []
      self.w_sthr = []
      for wi, nb in enumerate(self.static_tree.query_ball_point(
          o_np, self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
          return_sorted=False)):
        own = self.w_own[wi]
        if self.w_site_ops[wi]:
          # The water's own equivalents are the water itself.
          idx = flex.size_t([int(j) for j in nb if j not in own
                             and not self._own_equiv(wi, j)])
        else:
          idx = flex.size_t([int(j) for j in nb if j not in own])
        c = self.w_conf[wi]
        if c:
          cj = self.conf.select(idx)
          idx = idx.select(((cj == 0) | (cj == c)).iselection())
        self.w_sxyz.append(self.static_xyz.select(idx))
        self.w_sthr.append(self.static_thr.select(idx))
      # A placed H sits within oh_length of its own O, so only waters whose
      # O lie within clearance + 2 oh_length can hold one near this water's
      # candidates.
      r_wh = _WATER_CLEARANCE_RADIUS + 2.0 * self.oh_length + 0.01
      split = any(self.w_conf)
      self.w_wnbr = []
      for wi, nb in enumerate(KDTree(o_np).query_ball_point(
          o_np, r_wh, return_sorted=False)):
        self.w_wnbr.append([j for j in nb if j != wi and (
                              not split or _compatible(self.w_conf[wi],
                                                       self.w_conf[j]))])
      # The same test against the symmetry equivalents of these waters.
      if self.crystal_symmetry is not None:
        self.w_sym_nbr, self.sym_tf, own_blocks = _water_sym_equiv_neighbours(
          o_pts, self.crystal_symmetry, r_wh,
          self.w_conf if split else None)
        # matrix_cart's raw tuple is what multiplies a vec3_double array;
        # the matrix.sqr wrapper _pool_write uses does not.
        self.w_self_ops = []
        for wi in range(n):
          self.w_self_ops.append(tuple(
            (self.sym_tf[b][0].elems, self.sym_tf[b][1].elems)
            for b in own_blocks[wi]))
        if len(self.sym_tf) > 1:
          in_block = [set() for _ in range(n)]
          for wi in range(n):
            for b, wj in self.w_sym_nbr[wi]:
              in_block[wj].add(b)
          self.w_blocks = []
          for bs in in_block:
            self.w_blocks.append(tuple((b,) + self.sym_tf[b]
                                       for b in sorted(bs)))
          self.pool_xyz = flex.vec3_double(
            len(self.sym_tf) * self.n_slots, (0.0, 0.0, 0.0))
    return waters

  def _greedy_pass(self, waters, stripped):
    """Place every water once, in order, each against the H placed so far."""
    from mmtbx.ligands.hierarchy_utils import _new_atom
    hier = self.hier
    n = len(waters)
    # Initial greedy pass over the ordered waters, each avoiding the H
    # already placed on earlier ones. ``records`` keeps per-water (index,
    # placed-H slots, fixed_d1) for the refinement sweeps.
    self.slots_of = [()] * n

    # A bare water takes D when the model has D anywhere, the stripped
    # protons included; resolved at the first bare water.
    model_has_d = None
    for wi in range(n):
      (ag, o, own_idx, fixed_d1, existing, existing_names, wid, c, occ,
       carried) = waters[wi]
      if self.element is not None:
        proton_element = self.element
      elif existing:
        proton_element = existing[0].element.strip().upper()
      elif carried:
        proton_element = carried[0]
      else:
        if model_has_d is None:
          model_has_d = _has_deuterium(hier)
          for removed in stripped.values():
            for e, _ in removed:
              if e == "D":
                model_has_d = True
        proton_element = "D" if model_has_d else "H"
      h1, h2 = self._place_one(wi, self._nbr_slots(wi), fixed_d1)

      # Completing builds only the cone proton; the name takes the free slot.
      slots = []
      for di in ((2,) if fixed_d1 is not None else (1, 2)):
        proton_name = _free_proton_name(existing_names, proton_element)
        if proton_name is None:
          continue
        existing_names.add(proton_name.strip())
        xyz = h1 if di == 1 else h2
        atom = _new_atom(proton_name, proton_element, xyz, occ, o.b,
                         o.hetero, segid=o.segid)
        ag.append_atom(atom)
        slot = len(self.placed_coords)
        slots.append((atom, slot, di))
        self.placed_coords.append(xyz)
        self._pool_write(wi, slot, xyz)
        self.slot_wid[slot] = wid
        self.slot_conf[slot] = c
        self.slot_anchor[slot] = o.xyz
      if slots:
        self.records.append((wi, slots, fixed_d1))
        self.slots_of[wi] = tuple(s for _, s, _ in slots)

    # The placed-H neighbourhood of every water, now that all of them exist.
    self.w_nbr_slots = [self._nbr_slots(wi) for wi in range(n)]

  def _nbr_slots(self, wi, skip=()):
    """Pool slots of the waters and equivalents neighbouring water ``wi``,
    leaving out those of the waters in ``skip``."""
    nbr = []
    for wj in self.w_wnbr[wi]:
      if wj not in skip:
        nbr.extend(self.slots_of[wj])
    for b, wj in self.w_sym_nbr[wi]:
      if wj not in skip:
        off = b * self.n_slots
        nbr.extend(off + s for s in self.slots_of[wj])
    return flex.size_t(nbr)

  def _refine_and_hop(self):
    """Refine and basin-hop, keep the best state and return its label."""
    # Refinement sweeps: re-place each water against all other placed H,
    # stopping once a sweep removes fewer than refine_tol close (< 2.0 A)
    # contacts (0: all n_refine sweeps). The best state is kept.
    stats = None
    if self.on_state or self.n_refine or self.n_basin:
      stats = self._stats()
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
      if (self.refine_tol and prev_n20 is not None
          and prev_n20 - stats[1] < self.refine_tol):
        break  # gain below tolerance: converged
      prev_n20 = stats[1]

    # Basin-hopping (optional): each round restarts from the best state,
    # kicks the still-clashing waters to a random orientation, relaxes, and
    # keeps the result if it improved; it stops once no kick moves a water.
    # Deterministic (seeded).
    if self.n_basin and best_coords is not None:
      rng = flex.mersenne_twister(seed=_WATER_BASIN_SEED)
      for it in range(self.n_basin):
        self._restore(best_coords)
        moved = False
        for ri in self._clashing_records():
          if self._kick(ri, rng):
            moved = True
        if not moved:
          break  # no clashing water can move
        self.basin_rounds_run += 1
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

  def _search_clusters(self):
    """Re-choose the orientations of each cluster of waters jointly.

    Each water of a cluster (see :meth:`_water_clusters`) has its candidates
    from :meth:`_cluster_candidates`. The search keeps the combination of
    lowest total score that :meth:`_cluster_allowed` allows: exhaustive up to
    ``_WATER_CLUSTER_EXHAUSTIVE`` combinations, else by iterated conditional
    modes. Clusters go in record order, each against the waters outside it
    as they stand.
    """
    clusters, near = self._water_clusters()
    slots_of = {wi: slots for wi, slots, _fixed in self.records}
    for members in clusters:
      skip = set(members)
      cands = []
      scores = []
      for wi in members:
        pts, score = self._cluster_candidates(wi, slots_of[wi],
                                              self._nbr_slots(wi, skip))
        cands.append(pts)
        scores.append(score)
      allowed = self._cluster_allowed(members, cands, near)
      n = 1
      for s in scores:
        n *= s.size()
      if n <= _WATER_CLUSTER_EXHAUSTIVE:
        pick = _search_exhaustive(scores, allowed)
      else:
        pick = _search_icm(scores, allowed)
      for wi, pts, a in zip(members, cands, pick):
        if a:
          self._store(wi, slots_of[wi], pts[2 * a], pts[2 * a + 1])
    if self.on_state is not None:
      self.on_state("cluster search", self._stats())

  def _water_clusters(self):
    """The clusters of waters :meth:`_search_clusters` re-chooses.

    A cluster is a connected set of two or more waters whose O lie within
    ``_WATER_ACCEPTOR_RADIUS`` of each other, symmetry equivalents included.
    Only bare waters off a two-fold join one.

    Returns ``(clusters, near)``: the clusters, each a list of water indices
    in record order, and per clustered water its ``(block, water)`` pairs in
    reach, block 0 for the untransformed water, a water whose O coincides
    with its own left out (no direction to donate along, as in
    :meth:`_water_geom`).
    """
    free = [wi for wi, slots, fixed_d1 in self.records
            if fixed_d1 is None and self.w_axis[wi] is None
            and len(slots) == 2]
    near = {wi: [] for wi in free}
    linked = {wi: set() for wi in free}
    for wi in free:
      oi = self.w_o_col[wi]
      pairs = [(0, wj, self.w_o_col[wj]) for wj in self.w_wnbr[wi]
               if wj in near]
      for b, wj in self.w_sym_nbr[wi]:
        if wj in near:
          rot, trn = self.sym_tf[b]
          pairs.append((b, wj, rot * self.w_o_col[wj] + trn))
      for b, wj, oj in pairs:
        d = (oj - oi).length()
        if d > _WATER_ACCEPTOR_RADIUS:
          continue
        linked[wi].add(wj)
        linked[wj].add(wi)
        if d >= 1e-3:
          near[wi].append((b, wj))
    clusters = []
    seen = set()
    for wi in free:
      if wi in seen:
        continue
      seen.add(wi)
      members = [wi]
      todo = [wi]
      while todo:
        for wj in linked[todo.pop()]:
          if wj not in seen:
            seen.add(wj)
            members.append(wj)
            todo.append(wj)
      if len(members) > 1:
        clusters.append(sorted(members))
    return clusters, near

  def _cluster_candidates(self, wi, slots, fixed):
    """Candidate orientations of water ``wi`` for the cluster search.

    The water's current orientation first, then the even split over each
    acceptor pair whose protons both lie within the donation limit of their
    acceptors, then for each acceptor H1 on it with H2 where
    :meth:`_place_one` puts it. A new candidate must be away from cations and
    clear of the static atoms and of the protons on the pool slots
    ``fixed``, those of the waters outside the cluster; one repeating an
    earlier candidate is dropped. ``slots`` are the water's own record
    slots.

    Returns ``(pts, score)``: the candidates' protons as consecutive pairs,
    and per candidate the summed distance of its two protons to their
    nearest acceptor, each at most ``_WATER_CLUSTER_CAP``.
    """
    g = self._water_geom(wi)
    o = g["o"]
    acc_dirs = g["acc_dirs"]
    oh = self.oh_length
    new = flex.vec3_double()
    na = len(acc_dirs)
    for k1 in range(na):
      a1 = acc_dirs[k1]
      for k2 in range(k1 + 1, na):
        a2 = acc_dirs[k2]
        dirs = _split_dirs(a1, a2, self.cos_half_hoh, self.sin_half_hoh)
        if dirs is None:
          continue
        d1, d2 = dirs
        if (d1[0] * a1[0] + d1[1] * a1[1] + d1[2] * a1[2] < _WATER_H1_ACC_ALIGN
            or d2[0] * a2[0] + d2[1] * a2[1] + d2[2] * a2[2]
            < _WATER_H1_ACC_ALIGN):
          continue
        for d in dirs:
          new.append(tuple(o[c] + oh * d[c] for c in range(3)))
    for a in acc_dirs:
      new.extend(flex.vec3_double(self._place_one(wi, fixed, fixed_d1=a)))
    ok = _cat_ok(g["cat"], new, o)
    if self.w_self_ops[wi]:
      # Each candidate against its own equivalents.
      clear = flex.bool()
      for c in range(0, new.size(), 2):
        pair = new[c:c + 2]
        clear.extend(self._clear(wi, pair, fixed,
                                 own_fixed=(pair[0], pair[1])).ok)
      ok &= clear
    else:
      ok &= self._clear(wi, new, fixed).ok

    current = {di: self.placed_coords[slot] for _, slot, di in slots}
    pts = flex.vec3_double([current[1], current[2]])
    point = flex.size_t_range(new.size())
    both = ok.select(point - point % 2) & ok.select(point - point % 2 + 1)
    pts.extend(new.select(both))
    pts = _unrepeated_pairs(pts)

    acc = self.static_xyz.select(flex.size_t(g["acceptors"]))
    k = pts.size()
    m = acc.size()
    ik, im = _tile(k, m)
    dx, dy, dz = (pts.select(ik) - acc.select(im)).parts()
    dist = flex.sqrt(_row_mins(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz),
                               k, m, flex.size_t_range(k)))
    dist.set_selected(dist > _WATER_CLUSTER_CAP, _WATER_CLUSTER_CAP)
    score = (dist.select(flex.size_t_range(0, k, 2))
             + dist.select(flex.size_t_range(1, k, 2)))
    return pts, score

  def _cluster_allowed(self, members, cands, near):
    """Which candidates of each pair of a cluster's waters may coexist.

    Two candidates may not when a proton of one comes within
    ``_WATER_MIN_H_CLEARANCE`` of a proton of the other, or when each aims a
    proton within the donation limit at the other's O, symmetry equivalents
    included. Both waters' current orientations always may. ``members`` are
    the cluster's waters, ``cands`` their candidate protons as consecutive
    pairs and ``near`` the pairs in reach from :meth:`_water_clusters`.

    Returns the pair tables as :func:`_search_exhaustive` takes them.
    """
    pos = {wi: i for i, wi in enumerate(members)}
    n = [c.size() // 2 for c in cands]
    allowed = {}

    def forbid(i, a, j, b):
      """Rule out candidate ``a[k]`` of member ``i`` with ``b[k]`` of ``j``,
      for every ``k``."""
      if i > j:
        i, a, j, b = j, b, i, a
      table = allowed.get((i, j))
      if table is None:
        table = allowed[(i, j)] = flex.bool(n[i] * n[j], True)
      table.set_selected(a * n[j] + b, False)

    def equiv(b, wj):
      """Candidate protons and O of member ``wj`` under block ``b``."""
      pts = cands[pos[wj]]
      o = self.w_o_col[wj]
      if b:
        rot, trn = self.sym_tf[b]
        pts = rot.elems * pts + trn.elems
        o = rot * o + trn
      return pts, o

    thr = _WATER_MIN_H_CLEARANCE ** 2
    for i, wi in enumerate(members):
      pts = cands[i]
      partners = [(0, wj) for wj in self.w_wnbr[wi]
                  if wj in pos and pos[wj] > i]
      partners.extend((b, wj) for b, wj in self.w_sym_nbr[wi] if wj in pos)
      for b, wj in partners:
        other = equiv(b, wj)[0]
        k = other.size()
        ik, im = _tile(pts.size(), k)
        dx, dy, dz = (pts.select(ik) - other.select(im)).parts()
        bad = ((flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)) < thr
               ).iselection()
        if bad.size():
          forbid(i, bad / k / 2, pos[wj], bad % k / 2)
      oi = self.w_o_col[wi]
      for b, wj in near[wi]:
        other, oj = equiv(b, wj)
        gives = _donating(pts, oi, oj, self.oh_length).iselection()
        takes = _donating(other, oj, oi, self.oh_length).iselection()
        if gives.size() and takes.size():
          ia, ib = _tile(gives.size(), takes.size())
          forbid(i, gives.select(ia), pos[wj], takes.select(ib))
    for table in allowed.values():
      table[0] = True
    return allowed


def place_water_hydrogens(hier, oh_length=None, element=None, n_refine=5,
                          refine_tol=1, n_basin=0, cluster_search=False,
                          existing_h="keep", on_state=None,
                          crystal_symmetry=None, geometry=None):
  """Place the two H on every bare water, H-bond-aware.

  For each water residue missing H (any common water alias: HOH, DOD, H2O,
  WAT, OH2, ...): H1 along O -> the nearest acceptor giving a clash-free H
  (``_WATER_ACCEPTOR_RADIUS``, ``_WATER_ACCEPTOR_ELEMENTS``, excluding an N
  that carries an H or is bonded to three heavy atoms in a plane), else the
  max-clearance direction over a dense sphere;
  H2 on the H-O-H cone about O-H1, at a clash-free angle toward a second
  acceptor when there is one (outside the waters if it can reach both
  kinds), else the clearest; with both on acceptors, the water turns so
  each is equally far off its acceptor, where that stays clash-free.
  Candidates are scored against every other atom and every other water's
  placed H; waters go most-crowded first and new H inherit the parent O's
  B factor and their conformer's occupancy, divided by n for an O at
  occupancy 1 on a site of order n.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Single-model hierarchy; modified in place. Every atom must carry an
      element symbol.
  oh_length : float or None, optional
      O-H bond length in A, positive, overriding the one ``geometry`` sets.
  element : str or None, optional
      Element of the placed atoms, ``"H"`` or ``"D"``, forced on every water;
      None (default) takes the element of the water's own H/D, those
      ``existing_h="reorient"`` strips included, else ``"D"`` if the model
      has D atoms or DOD waters and ``"H"`` otherwise. Names follow it:
      ``H1``/``H2`` or ``D1``/``D2``.
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
  cluster_search : bool, optional
      After refinement and basin-hopping, re-choose the orientations of each
      cluster of waters whose O lie within H-bonding reach of each other
      jointly (default False). Each water's candidates are its current
      orientation, the even splits over two acceptors and H1 on each
      acceptor; the search keeps the combination whose protons lie closest to
      acceptors, without a new H-H clash or two waters donating to each
      other.
  existing_h : str, optional
      What to do with a water that already carries H: ``"keep"`` (default,
      leave it untouched, single-H waters included), ``"complete"`` (a water
      with one H gets its partner on that proton's cone, inheriting its
      element) or ``"reorient"`` (strip all water H and re-place both). An H
      of an O at occupancy 1 on a site of order n counts as n times its
      occupancy.
  on_state : callable, optional
      Called as ``on_state(label, stats)`` at each state reached:
      ``"initial"`` after the greedy pass, ``"sweep N"`` after each refinement
      sweep, ``"basin N.M"`` during basin-hopping and ``"cluster search"``
      after the cluster search, on the kept state; ``stats`` is the
      ``_water_clash_stats`` tuple, under ``crystal_symmetry``.
  crystal_symmetry : cctbx.crystal.symmetry or None, optional
      Symmetry-equivalent atoms within reach of a water join its
      environment, placed protons included, and contacts with those protons
      count towards the kept state. None (default) treats the model as
      isolated.
  geometry : str or None, optional
      O-H length and H-O-H angle, a key of ``_WATER_GEOMETRY``: ``"xray"``
      (0.850 A, 103.91 deg) or ``"neutron"`` (0.980 A, 103.91 deg), cctbx's
      restraint targets, or ``"gas_phase"`` (0.957 A, 104.5 deg). None
      (default) picks neutron if the model contains D or a DOD water, else
      X-ray.

  Returns
  -------
  libtbx.group_args
      ``kept_label``: label of the kept state (``"initial"``, ``"sweep N"`` or
      ``"basin N.M"``) when refinement or basin-hopping ran, else None.
      ``partial_waters``: one ``(residue_id, metal, action)`` per water that
      carried one H on input, counted as for ``existing_h``, ``metal`` being
      ``(element, distance)`` for a coordinating cation or None and
      ``action`` one of ``"kept"``, ``"completed"``, ``"stripped"``.
      ``basin_rounds_run``: basin-hopping rounds run; fewer than ``n_basin``
      once no clashing water in the best state can move.

  Raises
  ------
  AssertionError
      If ``hier`` has more than one model, or an atom without an element
      symbol.
  libtbx.utils.Sorry
      If an atom group of ``hier`` has the improper altloc ``" "``.
  """
  placer = _WaterHydrogenPlacer(
    hier, oh_length=oh_length, element=element, n_refine=n_refine,
    refine_tol=refine_tol, n_basin=n_basin, cluster_search=cluster_search,
    existing_h=existing_h, on_state=on_state,
    crystal_symmetry=crystal_symmetry, geometry=geometry)
  kept_label = placer.run()
  return group_args(kept_label=kept_label,
                    partial_waters=placer.partial_waters,
                    basin_rounds_run=placer.basin_rounds_run)


def _water_h_contacts(xyz, wid, conf=None):
  """Inter-water H-H contacts within 2.0 A.

  ``xyz`` are water H sites, ``wid`` the water each belongs to and ``conf``
  their conformer indices, None when no water H is in an altloc. Returns
  ``(i, j, d)``: the two H of each contact, H on the same water and H that
  never coexist (see :func:`_compatible`) excluded, and their distance.
  """
  from scipy.spatial import KDTree
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
  keep = wid.select(i) != wid.select(j)
  if conf is not None:
    ci, cj = conf.select(i), conf.select(j)
    keep &= (ci == 0) | (cj == 0) | (ci == cj)
  keep = keep.iselection()
  i = i.select(keep)
  j = j.select(keep)
  dx, dy, dz = (xyz.select(i) - xyz.select(j)).parts()
  return i, j, flex.sqrt(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz))


def _water_h_sym_pairs(xyz, wid, anchor, crystal_symmetry, conf=None):
  """Water H pairs that crystal symmetry can bring within 2.0 A.

  ``xyz`` are water H sites, ``wid`` the water each belongs to, ``anchor``
  the O site each is bound to (see :func:`_h_anchor`), None for an H left
  out, and ``conf`` their conformer indices, None when none is in an altloc.
  Two waters pair up when an operator brings the equivalent of one's O
  within ``2.0 + 2 * reach`` A of the other's, ``reach`` being the longest
  H-anchor distance; a water split between altlocs has an O per conformer.
  Each contact is listed once: a pair of waters and its mirror, the second
  against the first under the inverse operator, are one contact seen from
  either end, of which the first in ``(i, j, op)`` order is kept, and a
  water paired with itself by an operator that is its own inverse lists
  each pair of its H once. H that never coexist (see :func:`_compatible`)
  are not paired.

  Returns ``[(op, i, j)]``, one entry per operator: H ``i[k]`` against the
  equivalent of H ``j[k]`` under ``op``.
  """
  sel = flex.size_t([k for k in range(wid.size()) if anchor[k] is not None])
  if not sel.size():
    return []
  dx, dy, dz = (xyz.select(sel)
                - flex.vec3_double([anchor[k] for k in sel])).parts()
  reach = math.sqrt(flex.max(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)))
  h_of = {}
  for k in sel:
    w = wid[k]
    if w not in h_of:
      h_of[w] = []
    h_of[w].append(k)
  sites = sorted({(wid[k], anchor[k]) for k in sel})
  found = _sym_equiv_pairs(
    flex.vec3_double([o for _, o in sites]), crystal_symmetry,
    2.0 + 2.0 * reach + 0.01)
  # Water pairs; a split water's two O can pair it twice under one operator.
  pairs = []
  listed = set()
  for si, sj, op in found:
    pair = (sites[si][0], sites[sj][0], op)
    if pair not in listed:
      listed.add(pair)
      pairs.append(pair)
  by_op = {}
  for wi, wj, op in pairs:
    mirror = (wj, wi, op.inverse().new_denominators(op))
    if mirror in listed and mirror < (wi, wj, op):
      continue
    key = str(op)
    if key not in by_op:
      by_op[key] = (op, [], [])
    _, i, j = by_op[key]
    own_inverse = mirror == (wi, wj, op)
    for a in h_of[wi]:
      for b in h_of[wj]:
        if own_inverse and b < a:
          continue   # the same contact as (b, a)
        if conf is not None and not _compatible(conf[a], conf[b]):
          continue
        i.append(a)
        j.append(b)
  return [(op, flex.size_t(i), flex.size_t(j))
          for _, (op, i, j) in sorted(by_op.items())]


def _water_h_sym_contacts(xyz, sym_pairs, unit_cell):
  """The pairs of :func:`_water_h_sym_pairs` within 2.0 A.

  Returns ``[(op, i, j, d)]`` for the operators with any: H ``i[k]`` lies
  ``d[k]`` from the equivalent of H ``j[k]`` under ``op``.
  """
  contacts = []
  for op, i, j in sym_pairs:
    dx, dy, dz = (xyz.select(i) - super_cell.sym_equiv_sites_cart(
      sites_cart=xyz, unit_cell=unit_cell, rt_mx=op, selection=j)).parts()
    d2 = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
    near = (d2 <= 4.0).iselection()
    if near.size():
      contacts.append((op, i.select(near), j.select(near),
                       flex.sqrt(d2.select(near))))
  return contacts


def _contact_stats(n, d):
  """``(n, n_lt_20, n_lt_18, n_lt_15, closest)`` for ``n`` water H with
  contact distances ``d``; ``closest`` is None without a contact."""
  if not d.size():
    return n, 0, 0, 0, None
  return (n, d.size(), (d < 1.8).count(True), (d < 1.5).count(True),
          flex.min(d))


def _water_h_sites(hier):
  """``(xyz, wid, atoms, anchor, conf)`` for every water H/D in ``hier``:
  ``wid`` numbers the water residue each belongs to, ``anchor`` holds the O
  site each is bound to (see :func:`_h_anchor`), None without one, and
  ``conf`` the H's conformer indices, None when none is in an altloc."""
  altloc_index = hier.get_conformer_indices().index_altloc_mapping
  xyz = flex.vec3_double()
  wid = flex.size_t()
  conf = flex.size_t()
  atoms = []
  anchor = []
  for w, (rg, ags) in enumerate(_water_residue_groups(hier)):
    o_sites = _water_o_sites(ags)
    for ag in ags:
      ats, hd = _hd_flags(ag)
      sel = hd.iselection()
      xyz.extend(ats.extract_xyz().select(sel))
      wid.extend(flex.size_t(sel.size(), w))
      conf.extend(flex.size_t(sel.size(), altloc_index.get(ag.altloc, 0)))
      for k in range(ats.size()):
        if hd[k]:
          a = ats[k]
          atoms.append(a)
          anchor.append(_h_anchor(a.xyz, ag.altloc, o_sites))
  if (conf == 0).all_eq(True):
    conf = None
  return xyz, wid, atoms, anchor, conf


def _all_water_h_contacts(hier, crystal_symmetry):
  """Inter-water H-H contacts within 2.0 A, symmetry equivalents included.

  Returns ``(atoms, contacts)``: the water H atoms, and one ``(d, i, j, op)``
  per contact, between ``atoms[i]`` and the equivalent of ``atoms[j]`` under
  ``op``, None for a contact within the model. Without a crystal symmetry
  the model is isolated.
  """
  xyz, wid, atoms, anchor, conf = _water_h_sites(hier)
  i, j, d = _water_h_contacts(xyz, wid, conf)
  contacts = [(d[k], i[k], j[k], None) for k in range(d.size())]
  if crystal_symmetry is not None:
    for op, si, sj, sd in _water_h_sym_contacts(
        xyz, _water_h_sym_pairs(xyz, wid, anchor, crystal_symmetry, conf),
        crystal_symmetry.unit_cell()):
      contacts.extend((sd[k], si[k], sj[k], op) for k in range(sd.size()))
  return atoms, contacts


def _water_clash_stats(hier, crystal_symmetry=None):
  """Count H-H contacts between the H of different waters.

  Given a crystal symmetry, contacts with symmetry-equivalent water H count
  too, each once.

  Returns ``(n_placed, n_lt_20, n_lt_18, n_lt_15, closest)``: the number of
  water H, the counts of inter-water H-H contacts below 2.0/1.8/1.5 A, and
  the closest such distance (None if no pair is within 2.0 A).
  """
  atoms, contacts = _all_water_h_contacts(hier, crystal_symmetry)
  return _contact_stats(len(atoms), flex.double([c[0] for c in contacts]))


def _atom_id(a):
  """Compact atom identity, e.g. ``"HOH A 863 H2"`` or
  ``"HOH A 863 H2 altloc B"``."""
  L = a.fetch_labels()
  alt = L.altloc.strip()
  return (f"{L.resname.strip()} {L.chain_id.strip()} "
          f"{L.resseq.strip()} {L.name.strip()}"
          + (f" altloc {alt}" if alt else ""))


def _water_id(ag):
  """Compact water identity, e.g. ``"HOH A 863"`` or
  ``"HOH A 863 altloc B"``."""
  rg = ag.parent()
  alt = ag.altloc.strip()
  return (f"{ag.resname.strip()} {rg.parent().id.strip()} {rg.resseq.strip()}"
          + (f" altloc {alt}" if alt else ""))


def _worst_water_clashes(hier, crystal_symmetry=None):
  """Inter-water H-H contacts within 2.0 A, closest first.

  One ``(distance, id_a, id_b)`` per contact, as counted by
  :func:`_water_clash_stats`; ``id_b`` of a symmetry equivalent ends in its
  operator, in brackets.
  """
  atoms, contacts = _all_water_h_contacts(hier, crystal_symmetry)
  return [(d, _atom_id(atoms[i]),
           _atom_id(atoms[j]) + (f" [{op}]" if op is not None else ""))
          for d, i, j, op in sorted(contacts, key=lambda c: c[0])]

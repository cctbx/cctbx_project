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
``mmtbx.development.naiad`` command line.

Needs ``scipy`` (KDTree); the vectorised work is
``scitbx.array_family.flex``.
"""

from __future__ import absolute_import, division, print_function

import math
import random

import iotbx.pdb
from cctbx.crystal import super_cell
from iotbx.pdb.utils import check_for_missing_elements
from libtbx import group_args
from scitbx import matrix
from scitbx.array_family import flex
from scipy.spatial import KDTree


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
_EMPTY_SLOTS = flex.size_t()


def _as_vec3(tuples):
  """Build a flex.vec3_double, tolerating the empty case."""
  return flex.vec3_double(tuples) if tuples else flex.vec3_double()


# Index arrays for the (candidate x neighbour) block, keyed by its shape and
# shared by every water with that shape. The block is flattened row-major:
# candidate i occupies ``[i * m, (i+1) * m)``, so a row is a plain slice.
_TILE_CACHE = {}


def _tile(k, m):
  """``(cand_index, nbr_index)`` for a flattened (k, m) block."""
  t = _TILE_CACHE.get((k, m))
  if t is None:
    r = flex.size_t_range(k * m)
    t = (r / m, r % m)
    _TILE_CACHE[(k, m)] = t
  return t


def _row_mins(d, k, m, rows):
  """Per-row minima of a flattened (k, m) block, for the rows in ``rows``."""
  n = rows.size()
  if n == 0:
    return flex.double()
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
  """Rows of maximal ``2 * cat + ok`` rank, as ``(top, rank)``."""
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


# The H2 cone is sampled at fixed azimuths, so its cosines and sines are
# constants, as is the normalized fallback sphere.
_CONE_COS = flex.double([math.cos(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_CONE_SIN = flex.double([math.sin(2.0 * math.pi * k / _WATER_CONE_SAMPLES)
                         for k in range(_WATER_CONE_SAMPLES)])
_FALLBACK_DIRS = [matrix.col(v).normalize().elems
                  for v in _WATER_FALLBACK_DIRECTIONS]
_FALLBACK_FX = _as_vec3(_FALLBACK_DIRS)


def _ortho_frame(d1):
  """Two unit vectors ``p, q`` making ``(d1, p, q)`` an orthonormal frame.

  Mirrors ``scitbx.matrix.col.ortho()`` operation for operation.
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
  """A random unit vector, from Gaussian deviates of ``rng``."""
  while True:
    v = (rng.gauss(0, 1), rng.gauss(0, 1), rng.gauss(0, 1))
    n = math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    if n > 1e-6:
      return (v[0] / n, v[1] / n, v[2] / n)


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


def _compatible(ca, cb):
  """Whether atoms of conformer indices ``ca`` and ``cb`` coexist.

  cctbx's rule: atoms in two different altlocs never meet, and a blank atom
  (index 0) meets every altloc.
  """
  return not ca or not cb or ca == cb


def _sp2_plane_normal(atoms, static_tree, c, exclude, conf=None, view=0):
  """Unit normal of the sp2 plane around atom ``c``.

  Built from ``c``'s heavy substituents other than ``exclude``, the bonded
  acceptor O that supplies one in-plane vector, taking only substituents
  that coexist with conformer ``view`` (any, for view 0; ``conf`` holds the
  per-atom conformer indices). None if underdetermined.
  """
  C = matrix.col(atoms[c].xyz)
  oc = matrix.col(atoms[exclude].xyz) - C   # C -> O
  for k in static_tree.query_ball_point(atoms[c].xyz, _WATER_BOND_HEAVY):
    if k == c or k == exclude:
      continue
    if view and not _compatible(view, conf[k]):
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


def _bond_lobes(atoms, static_tree, i, el, nbrs, conf, view):
  """Lone-pair lobes of acceptor ``i`` (element ``el``) from its bonded
  neighbours ``nbrs``, the sp2 plane taken in conformer ``view``; see
  :func:`_acceptor_lobes`."""
  A = matrix.col(atoms[i].xyz)
  bond_dirs = [(matrix.col(atoms[j].xyz) - A).normalize() for j in nbrs]
  if not bond_dirs:
    return []
  if el == "O" and len(nbrs) == 1:
    away = bond_dirs[0] * -1.0
    n = _sp2_plane_normal(atoms, static_tree, nbrs[0], i, conf, view)
    if n is None:
      return [away]
    return [
      away.rotate_around_origin(axis=n, angle=_WATER_SP2_LOBE_DEG, deg=True),
      away.rotate_around_origin(axis=n, angle=-_WATER_SP2_LOBE_DEG, deg=True)]
  bsum = matrix.col((0.0, 0.0, 0.0))
  for b in bond_dirs:
    bsum = bsum + b
  return [(bsum * -1.0).normalize()] if bsum.length() > 1e-6 else []


def _acceptor_lobes(atoms, static_tree, donor_n, indices, conf=None):
  """Lone-pair lobe unit vectors per acceptor atom among ``indices``, from
  bonded geometry.

  - terminal O (1 bond, carbonyl/carboxylate): two in-plane sp2 lobes
    ``2 * _WATER_SP2_LOBE_DEG`` apart, straddling the direction away from the
    bonded atom.
  - otherwise: one lobe opposite the sum of the bond directions.

  ``donor_n`` holds the indices of N that carry an H (donors, not acceptors).
  ``conf`` holds per-atom conformer indices, None for a model without
  altlocs: a bond joins only atoms that coexist, and a blank acceptor whose
  bonded neighbours are split between altlocs has one geometry per altloc.

  Returns acceptor index -> ``{view: lobe list}``, a lobe list empty where
  the geometry is underdetermined. ``view`` is the conformer the lobes were
  built in: the acceptor's own, unless it is blank with split neighbours,
  which gives one view per altloc among them and view 0 from its blank
  neighbours alone. :func:`_view_lobes` picks the lobes a water sees.
  """
  lobes = {}
  for i in indices:
    a = atoms[i]
    el = a.element.strip().upper()
    if el not in _WATER_ACCEPTOR_ELEMENTS:
      continue
    if el == "N" and i in donor_n:
      continue  # protonated N is a donor, not an acceptor
    ci = conf[i] if conf is not None else 0
    A = matrix.col(a.xyz)
    nbrs = []
    for j in static_tree.query_ball_point(a.xyz, _WATER_BOND_HEAVY):
      if j == i:
        continue
      if conf is not None and not _compatible(ci, conf[j]):
        continue  # another altloc's copy of a neighbour, or of this atom
      d = (matrix.col(atoms[j].xyz) - A).length()
      lim = (_WATER_NH_BOND if atoms[j].element_is_hydrogen()
             else _WATER_BOND_HEAVY)
      # An atom on top of this one contributes no bond direction.
      if 1e-3 < d <= lim:
        nbrs.append(j)
    split = ({conf[j] for j in nbrs} - {0}) if conf is not None and not ci \
        else ()
    if not split:
      lobes[i] = {ci: _bond_lobes(atoms, static_tree, i, el, nbrs, conf, ci)}
    else:
      lobes[i] = {
        v: _bond_lobes(atoms, static_tree, i, el,
                       [j for j in nbrs if conf[j] in (0, v)], conf, v)
        for v in sorted(split | {0})}
  return lobes


def _view_lobes(views, c):
  """The lobes of ``views`` (see :func:`_acceptor_lobes`) that a water of
  conformer index ``c`` sees: its own altloc's, else the blank view's; a
  blank water sees every altloc's."""
  if not views:
    return ()
  if not c:
    split = [lobe for v in sorted(views) if v for lobe in views[v]]
    return split or views.get(0, ())
  return views.get(c, views.get(0, ()))


def _symmetry_environment(hier, sites_cart, crystal_symmetry, radius,
                          min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Copies of the atoms crystal symmetry places within ``radius`` of a water.

  Each copy carries its whole residue group, so an acceptor keeps the bonded
  neighbours its lobe geometry needs and an N keeps the H that marks it a
  donor. Copies are environment only: the waters to protonate come from the
  asymmetric unit, and a mate's placed protons are tracked separately (see
  :func:`_water_sym_equiv_neighbours`). An equivalent closer than
  ``min_distance_sym_equiv`` to its own site counts as that same atom, which
  the asymmetric unit already holds.

  Returns ``(hierarchies, atoms, xyz, source)``: the sub-hierarchies, which
  own the atom objects and must be kept alive, their atoms in coordinate
  order, their sites, and the index of the model atom each was copied from.
  All four are empty when symmetry places nothing in range.
  """
  nothing = [], [], flex.vec3_double(), flex.size_t()
  o_sel = flex.size_t([a.i_seq for a in hier.atoms()
                       if _is_water(a.parent().resname)
                       and a.element.strip().upper() == "O"])
  if not o_sel.size():
    return nothing
  siiu, _ = super_cell.get_siiu(
    sites_cart=sites_cart, crystal_symmetry=crystal_symmetry,
    select_within_radius=radius, selection=o_sel, buffer=0,
    min_distance_sym_equiv=min_distance_sym_equiv)
  if not siiu:
    return nothing
  # Group the residue groups by operator, keyed on str(op): get_siiu gives
  # every operator the same denominators, so equal operators print alike.
  hier_atoms = hier.atoms()
  rg_of = {}
  by_op = {}
  for j_seq, ops in siiu.items():
    rg = hier_atoms[j_seq].parent().parent()
    grp = rg_of.get(rg.memory_id())
    if grp is None:
      grp = rg_of[rg.memory_id()] = rg.atoms().extract_i_seq()
    for op in ops:
      by_op.setdefault(str(op), (op, set()))[1].update(grp)
  unit_cell = crystal_symmetry.unit_cell()
  hiers = []
  atoms = []
  xyz = flex.vec3_double()
  source = flex.size_t()
  for key in sorted(by_op):
    op, grp = by_op[key]
    sel = flex.size_t(sorted(grp))
    equiv = super_cell.sym_equiv_sites_cart(
      sites_cart=sites_cart, unit_cell=unit_cell, rt_mx=op, selection=sel)
    # A grown atom on a symmetry element of op lands on itself.
    dx, dy, dz = (equiv - sites_cart.select(sel)).parts()
    keep = (flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
            >= min_distance_sym_equiv ** 2).iselection()
    # copy_atoms: without it set_xyz moves the model's own atoms.
    sub = hier.select(sel.select(keep), copy_atoms=True)
    sub_atoms = sub.atoms()
    sub_atoms.set_xyz(equiv.select(keep))
    hiers.append(sub)
    atoms.extend(list(sub_atoms))
    xyz.extend(sub_atoms.extract_xyz())
    source.extend(sel.select(keep))
  return hiers, atoms, xyz, source


def _sym_equiv_pairs(sites, crystal_symmetry, radius,
                     min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Every ``(i, j, op)`` whose equivalent of site ``j`` lies within
  ``radius`` of site ``i``.

  ``op`` acts on fractional coordinates. A contact is listed from both ends,
  ``(j, i)`` under the inverse operator, except that a site on a special
  position can see more equivalents of its partner than the partner sees of
  it. An equivalent closer than ``min_distance_sym_equiv`` to its own site is
  that same site, and is dropped.
  """
  if not sites.size():
    return []
  # One table build serves the query for every site.
  tables = super_cell.get_sym_equiv_tables(
    sites_cart=sites, crystal_symmetry=crystal_symmetry, radius=radius,
    min_distance_sym_equiv=min_distance_sym_equiv)
  found = []
  for i in range(sites.size()):
    siiu, _ = super_cell.get_siiu(
      crystal_symmetry=crystal_symmetry, selection=[i],
      symmetry_tables=tables)
    for j, ops in siiu.items():
      for op in ops:
        found.append((i, j, op))
  return found


def _water_sym_equiv_neighbours(sites, crystal_symmetry, radius,
                            min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL,
                            conf=None):
  """Symmetry equivalents of ``sites`` that fall within ``radius`` of a site.

  ``conf`` holds per-site conformer indices; a pair whose sites never
  coexist (see :func:`_compatible`) is dropped. None keeps every pair.

  Returns ``(by_site, transforms, own)``: per site, the sorted ``(block, other
  site)`` pairs whose equivalent is in range; per block the cartesian
  ``(rotation, translation)`` that produces it, block 0 being an unused
  placeholder so a block index is never zero; and per site the blocks whose
  operator brings the site's own equivalent into range. An equivalent closer
  than ``min_distance_sym_equiv`` to its own site is that same site, and is
  dropped.
  """
  n = sites.size()
  by_site = [[] for _ in range(n)]
  own = [[] for _ in range(n)]
  transforms = [None]
  found = _sym_equiv_pairs(sites, crystal_symmetry, radius,
                           min_distance_sym_equiv)
  if conf is not None:
    found = [(wi, wj, op) for wi, wj, op in found
             if _compatible(conf[wi], conf[wj])]
  if not found:
    return by_site, transforms, own
  # Number the blocks by operator rather than by discovery order, so the slot
  # order a site sees does not depend on the order the tables list pairs in.
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


class _Clearance(object):
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
               min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL,
               geometry=None):
    self.hier = hier
    self.oh_length = oh_length
    self.geometry = geometry
    self.element = element
    self.n_refine = n_refine
    self.refine_tol = refine_tol
    self.n_basin = n_basin
    self.existing_h = existing_h
    self.lone_pair_directed = lone_pair_directed
    self.on_state = on_state
    self.crystal_symmetry = crystal_symmetry
    self.min_distance_sym_equiv = min_distance_sym_equiv

    # Placed-H coordinates, one slot per proton, filled in run(); placed_xyz
    # is the same data as a flex.vec3_double.
    self.placed_coords = []
    self.placed_xyz = None
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
    # Squared distances go through parts() and pow2, never .dot() or .norms():
    # those contract to fused multiply-add in C++, breaking exactness and
    # flipping threshold tests. Same for every squared distance in this file.
    dx, dy, dz = (cands.select(ik) - P.select(im)).parts()
    d = flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz)
    bad = (d < self.w_sthr[wi].select(im)).iselection()
    if bad.size():
      ok.set_selected(bad / m, False)
    return d, m, ok

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

    Every candidate lies exactly ``oh_length`` from the water O, so water
    ``wi``'s cached static block (a ball of
    ``oh_length + _WATER_CLEARANCE_RADIUS`` about the O, own atoms already
    dropped) covers every candidate's own clearance ball. The test is one dense
    (candidate x neighbour) block.

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
    # Squared distances throughout: the thresholds and the cap are exact
    # squares and sqrt is monotone.
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
          X, W, anchor, self.crystal_symmetry,
          self.min_distance_sym_equiv, C)
      for _op, _i, _j, d_sym in _water_h_sym_contacts(
          X, self.h_sym_pairs, self.crystal_symmetry.unit_cell()):
        d = d.concatenate(d_sym)
    return _contact_stats(X.size(), d)

  def _nearest_cation(self, o_xyz, own_idx, c):
    """Closest metal cation coordinating a water O, if any.

    Returns ``(element, distance)`` or None, over
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

  def _water_geom(self, wi):
    """Fixed per-water geometry, computed once.

    Returns the dict of ``o`` (O coordinates), ``acceptors`` (atom indices,
    nearest first), ``acc_dirs``/``acc_pts`` (O-H unit directions as tuples and
    the H positions built from them), ``cat`` (O->cation unit directions),
    ``sph`` (the lazily built fallback-sphere H positions) and
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
        for lobe in _view_lobes(self.acc_lobes.get(i), self.w_conf[wi]):
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
    """The two H positions for one water O, clash-aware.

    ``nbr_slots`` are the pool slots that may lie near this water, its own
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

    p, q = _ortho_frame(d1)

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
    # The alignment term only ever separates candidates that are both
    # cation-ok and clash-free (rank 3); below that it cannot break a tie.
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

    True if any of them moved, which is what :meth:`_dirty` reads.
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

    A water is placed against its static surroundings, which never move, and
    the placed H of the waters in ``w_wnbr[wi]`` and the symmetry
    equivalents in ``w_sym_nbr[wi]``, an equivalent moving exactly when its
    own water does. If none of those H has moved since this water was last
    placed, the placement re-derives the two positions it already holds. The
    monotonic tick is bumped once per water per sweep, so a neighbour that
    moves earlier in the same sweep still counts.
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
    """
    wi, slots, fixed_d1 = self.records[ri]
    o = self._water_geom(wi)["o"]
    d1 = fixed_d1 if fixed_d1 is not None else _rand_unit(rng)
    p, q = _ortho_frame(d1)
    theta = rng.uniform(0.0, 2.0 * math.pi)
    ct, st = math.cos(theta), math.sin(theta)
    d2 = tuple(self.cos_hoh * d1[c] + self.sin_hoh * (ct * p[c] + st * q[c])
               for c in range(3))
    self.tick += 1
    if self._store(wi, slots,
                   tuple(o[c] + self.oh_length * d1[c] for c in range(3)),
                   tuple(o[c] + self.oh_length * d2[c] for c in range(3))):
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
    if self.geometry is None:
      self.geometry = "neutron" if _has_deuterium(hier) else "xray"
    oh, self.hoh_deg = _WATER_GEOMETRY[self.geometry]
    if self.oh_length is None:
      self.oh_length = oh
    # What each water carried is read off the strip; the walk below sees the
    # same protons in every other mode.
    stripped = {}
    if self.existing_h == "reorient":
      stripped = _strip_water_hydrogens(hier)

    sel = hier.atoms()
    atoms = list(sel)
    if not atoms:
      return None
    sel.reset_i_seq()
    # Conformer index per atom, 0 for blank; which atoms coexist follows
    # cctbx's rule (see _compatible). An improper altloc (one atom name both
    # blank and in an altloc) would read as an altloc of its own.
    ci = hier.get_conformer_indices()
    if " " in ci.index_altloc_mapping:
      hier.overall_counts().raise_improper_alt_conf_if_necessary()
    conf = ci.conformer_indices
    altloc_index = dict(ci.index_altloc_mapping)

    # Gather the waters to protonate: one per water conformer, a residue's
    # blank atoms plus one altloc's (see _water_conformers), so a blank O
    # whose H sit in altlocs A and B is two complete waters. Every water H
    # is also listed once for the clash statistics, the residue being the
    # water it belongs to.
    waters = []
    # Single-H waters to report, once per such H, annotated once the tree
    # exists.
    single = {}
    wh_xyz = []
    wh_wid = []
    wh_conf = []
    wh_anchor = []
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
        carried = [(el, occ, g)
                   for g in ([target] if target is blank else [blank, target])
                   if g is not None
                   for el, occ in stripped.get(g.memory_id(), ())]
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
        is_single = (len(carried) == 1 if self.existing_h == "reorient"
                     else len(existing) == 1)
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
          action = ("stripped" if self.existing_h == "reorient"
                    else "completed" if fixed_d1 is not None else "kept")
          single.setdefault((rg.memory_id(), holder.memory_id()),
                            (_water_id(holder), o.xyz, own_idx, c, action))
        if skip:
          continue
        # New H take the conformer's occupancy: its own altloc's atoms', or
        # what the strip took from them, else the O's.
        occ = o.occ
        if altloc:
          tats = target.atoms()
          if tats.size():
            occ = tats[0].occ
          else:
            occ = next((q for _, q, g in carried if g is target), occ)
        waters.append((target, o, own_idx, fixed_d1, existing, names, wid, c,
                       occ, [el for el, _, _ in carried]))
    single = list(single.values())
    self.wh_xyz = _as_vec3(wh_xyz)
    self.wh_wid = flex.size_t(wh_wid) if wh_wid else flex.size_t()
    self.wh_conf = flex.size_t(wh_conf) if wh_conf else flex.size_t()
    self.wh_anchor = wh_anchor
    # Water H pairs across symmetry for the clash counts, built by _stats.
    self.h_sym_pairs = None

    # Neighbouring asymmetric units, as environment atoms appended after the
    # model's own. They keep the asymmetric unit's indices 0..n-1 valid as
    # both tree and i_seq, and the waters come from the walk above, so only
    # the model's waters are protonated. Only placement and the single-H
    # report read them.
    self.sym_hier = []
    sym_atoms = []
    sym_xyz = flex.vec3_double()
    sym_src = flex.size_t()
    if self.crystal_symmetry is not None and (waters or single):
      self.sym_hier, sym_atoms, sym_xyz, sym_src = _symmetry_environment(
        hier, sel.extract_xyz(), self.crystal_symmetry,
        max(self.oh_length + _WATER_CLEARANCE_RADIUS + 0.01,
            _WATER_ACCEPTOR_RADIUS), self.min_distance_sym_equiv)
    atoms = atoms + sym_atoms
    self.atoms = atoms
    # A copy keeps its source's altloc, whatever the operator, as in cctbx.
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

    # N atoms that already carry an H are donors, not acceptors (amide,
    # ammonium, guanidinium, protonated His ring N, ...). O always accepts,
    # so only N is filtered.
    self.donor_n = set()
    h_idx = self.static_is_h.iselection()
    if h_idx.size():
      for i, nbrs in zip(h_idx, self.static_tree.query_ball_point(
          self.static_xyz.select(h_idx).as_numpy_array(), _WATER_NH_BOND,
          return_sorted=False)):
        for j in nbrs:
          if j != i and atoms[j].element.strip().upper() == "N":
            self.donor_n.add(j)

    # Lone-pair lobe directions per acceptor (opt-in; empty when off), filled
    # once the waters' acceptors are known.
    self.acc_lobes = {}

    self.cos_hoh = math.cos(math.radians(self.hoh_deg))
    self.sin_hoh = math.sin(math.radians(self.hoh_deg))

    # The cation coordinating each single-H water, symmetry copies included.
    self.partial_waters = [
      (rid, self._nearest_cation(o_xyz, own_idx, c), action)
      for rid, o_xyz, own_idx, c, action in single]

    # Per-water constants: every neighbour list the placement needs, built
    # once here and reused by the greedy pass, every relaxation sweep and
    # every basin round.
    n = len(waters)
    self.w_geom = [None] * n
    self.placed_coords = []
    self.placed_xyz = flex.vec3_double(2 * n, (0.0, 0.0, 0.0))
    # Protons the clearance test may draw on: the placed ones, plus one
    # transformed block per symmetry operator that brings another water's
    # protons into reach. Block b slot s lives at b * n_slots + s, so a
    # symmetry-equivalent proton is just another slot. Without them the pool
    # is the placed array itself and nothing extra is paid.
    self.n_slots = 2 * n
    self.pool_xyz = self.placed_xyz
    self.sym_tf = [None]
    self.w_sym_nbr = [[] for _ in range(n)]
    self.w_blocks = [()] * n
    # Cartesian (rotation, translation) per operator that brings the water's
    # own equivalent into range.
    self.w_self_ops = [()] * n
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
      crowd = [sum(1 for j in nb if j not in waters[k][2]
                   and (not w_conf[k] or _compatible(w_conf[k], self.conf[j])))
               for k, nb in enumerate(self.static_tree.query_ball_point(
                 o_pts.as_numpy_array(), _WATER_CLEARANCE_RADIUS,
                 return_sorted=False))]
      order = sorted(range(n), key=lambda k: crowd[k], reverse=True)
      waters = [waters[k] for k in order]
      # Conformer index per water; the filters below run only for the
      # waters in an altloc, a blank water meeting every atom.
      self.w_conf = [w_conf[k] for k in order]
      o_pts = o_pts.select(flex.size_t(order))
      o_np = o_pts.as_numpy_array()
      self.w_own = [w[2] for w in waters]
      self.w_o_col = [matrix.col(w[1].xyz) for w in waters]
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
      if self.lone_pair_directed:
        near = set()
        for nb in self.w_acc_raw:
          near.update(nb)
        has_altlocs = not (self.conf == 0).all_eq(True)
        self.acc_lobes = _acceptor_lobes(
          atoms, self.static_tree, self.donor_n, near,
          self.conf if has_altlocs else None)
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
      self.w_wnbr = [[j for j in nb if j != wi and (
                        not split or _compatible(self.w_conf[wi],
                                                 self.w_conf[j]))]
                     for wi, nb in enumerate(KDTree(o_np).query_ball_point(
                       o_np, r_wh, return_sorted=False))]
      # The same test against the symmetry equivalents of these waters.
      if self.crystal_symmetry is not None:
        self.w_sym_nbr, self.sym_tf, own_blocks = _water_sym_equiv_neighbours(
          o_pts, self.crystal_symmetry, r_wh,
          self.min_distance_sym_equiv, self.w_conf if split else None)
        # matrix_cart's raw tuple is what multiplies a vec3_double array;
        # the matrix.sqr wrapper _pool_write uses does not.
        self.w_self_ops = [
          tuple((self.sym_tf[b][0].elems, self.sym_tf[b][1].elems)
                for b in own_blocks[wi]) for wi in range(n)]
        if len(self.sym_tf) > 1:
          in_block = [set() for _ in range(n)]
          for wi in range(n):
            for b, wj in self.w_sym_nbr[wi]:
              in_block[wj].add(b)
          self.w_blocks = [tuple((b,) + self.sym_tf[b] for b in sorted(bs))
                           for bs in in_block]
          self.pool_xyz = flex.vec3_double(
            len(self.sym_tf) * self.n_slots, (0.0, 0.0, 0.0))

    # Initial greedy pass over the ordered waters, each avoiding the H
    # already placed on earlier ones. ``records`` keeps per-water (index,
    # placed-H slots, fixed_d1) for the refinement sweeps.
    slots_of = [()] * n

    def nbr_slots(wi):
      """Pool slots of the waters and equivalents neighbouring water ``wi``."""
      nbr = [s for wj in self.w_wnbr[wi] for s in slots_of[wj]]
      for b, wj in self.w_sym_nbr[wi]:
        off = b * self.n_slots
        nbr.extend(off + s for s in slots_of[wj])
      return flex.size_t(nbr) if nbr else _EMPTY_SLOTS

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
        proton_element = "D" if ag.resname.strip().upper() == "DOD" else "H"
      h1, h2 = self._place_one(wi, nbr_slots(wi), fixed_d1)

      # Completing builds only the cone proton; the name takes the free slot.
      slots = []
      for di in ((2,) if fixed_d1 is not None else (1, 2)):
        proton_name = _free_proton_name(existing_names, proton_element)
        if proton_name is None:
          continue
        existing_names.add(proton_name.strip())
        xyz = h1 if di == 1 else h2
        atom = _new_h_atom(proton_name, proton_element, xyz, occ, o.b,
                           o.hetero)
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
                          min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL,
                          geometry=None):
  """Place the two H on every bare water, H-bond-aware.

  For each water residue missing H (any common water alias: HOH, DOD, H2O,
  WAT, OH2, ...): H1 along O -> the nearest acceptor giving a clash-free H
  (``_WATER_ACCEPTOR_RADIUS``, ``_WATER_ACCEPTOR_ELEMENTS``, N carrying an H
  excluded as donors), else the max-clearance direction over a dense sphere;
  H2 on the H-O-H cone about O-H1, at a clash-free angle toward
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
      O-H bond length in A, positive, overriding the one ``geometry`` sets.
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
      ``_water_clash_stats`` tuple, under ``crystal_symmetry``.
  crystal_symmetry : cctbx.crystal.symmetry or None, optional
      Honour crystal packing: atoms that symmetry places within reach of a
      water join its environment, so H at a lattice contact avoid the
      neighbouring asymmetric units instead of pointing into them. None
      (default) treats the model as isolated. A symmetry mate contributes its
      O and its placed protons both, and contacts with its protons count
      towards the kept state.
  min_distance_sym_equiv : float, optional
      Distance in A under which a symmetry equivalent counts as coincident
      with its own site, and so as that same atom rather than a second copy
      (default 0.5). A water refined a little off a symmetry element needs a
      larger value to be recognised as sitting on it.
  geometry : str or None, optional
      O-H length and H-O-H angle, a key of ``_WATER_GEOMETRY``: ``"xray"``
      (0.850 A, 103.91 deg) or ``"neutron"`` (0.980 A, 103.91 deg), cctbx's
      restraint targets, or ``"gas_phase"`` (0.957 A, 104.5 deg). None
      (default) picks neutron if the model contains D, else X-ray.

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
    min_distance_sym_equiv=min_distance_sym_equiv,
    geometry=geometry)
  kept_label = placer.run()
  return group_args(kept_label=kept_label,
                    partial_waters=placer.partial_waters)


def _water_h_contacts(xyz, wid, conf=None):
  """Inter-water H-H contacts within 2.0 A.

  ``xyz`` are water H sites, ``wid`` the water each belongs to and ``conf``
  their conformer indices, None when no water H is in an altloc. Returns
  ``(i, j, d)``: the two H of each contact, H on the same water and H that
  never coexist (see :func:`_compatible`) excluded, and their distance.
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
  keep = wid.select(i) != wid.select(j)
  if conf is not None:
    ci, cj = conf.select(i), conf.select(j)
    keep &= (ci == 0) | (cj == 0) | (ci == cj)
  keep = keep.iselection()
  i = i.select(keep)
  j = j.select(keep)
  dx, dy, dz = (xyz.select(i) - xyz.select(j)).parts()
  return i, j, flex.sqrt(flex.pow2(dx) + flex.pow2(dy) + flex.pow2(dz))


def _water_h_sym_pairs(xyz, wid, anchor, crystal_symmetry,
                       min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL,
                       conf=None):
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
    h_of.setdefault(wid[k], []).append(k)
  sites = sorted({(wid[k], anchor[k]) for k in sel})
  found = _sym_equiv_pairs(
    flex.vec3_double([o for _, o in sites]), crystal_symmetry,
    2.0 + 2.0 * reach + 0.01, min_distance_sym_equiv)
  # Water pairs; a split water's two O can pair it twice under one operator.
  pairs = list(dict.fromkeys((sites[si][0], sites[sj][0], op)
                             for si, sj, op in found))
  listed = set(pairs)
  by_op = {}
  for wi, wj, op in pairs:
    mirror = (wj, wi, op.inverse().new_denominators(op))
    if mirror in listed and mirror < (wi, wj, op):
      continue
    _, i, j = by_op.setdefault(str(op), (op, [], []))
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


def _all_water_h_contacts(hier, crystal_symmetry, min_distance_sym_equiv):
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
        xyz, _water_h_sym_pairs(xyz, wid, anchor, crystal_symmetry,
                                min_distance_sym_equiv, conf),
        crystal_symmetry.unit_cell()):
      contacts.extend((sd[k], si[k], sj[k], op) for k in range(sd.size()))
  return atoms, contacts


def _water_clash_stats(hier, crystal_symmetry=None,
                       min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Count H-H contacts between the H of different waters.

  Given a crystal symmetry, contacts with symmetry-equivalent water H count
  too, each once.

  Returns ``(n_placed, n_lt_20, n_lt_18, n_lt_15, closest)``: the number of
  water H, the counts of inter-water H-H contacts below 2.0/1.8/1.5 A, and
  the closest such distance (None if no pair is within 2.0 A).
  """
  atoms, contacts = _all_water_h_contacts(hier, crystal_symmetry,
                                          min_distance_sym_equiv)
  return _contact_stats(len(atoms), flex.double([c[0] for c in contacts]))


def _clash_row(label, stats, log):
  """Print one row of the per-sweep clash table for ``stats`` to ``log``."""
  _, n20, n18, n15, worst = stats
  w = f"{worst:.2f}" if worst is not None else ">2.0"
  print(f"  {label:<9} <2.0={n20:<5} <1.8={n18:<5} <1.5={n15:<5} closest={w}",
        file=log)


def _atom_id(a):
  """Compact atom identity, e.g. ``"HOH A 863 H2"`` (altloc in
  parentheses)."""
  L = a.fetch_labels()
  alt = L.altloc.strip()
  return (f"{L.resname.strip()} {L.chain_id.strip()} "
          f"{L.resseq.strip()} {L.name.strip()}"
          + (f" ({alt})" if alt else ""))


def _water_id(ag):
  """Compact water identity, e.g. ``"HOH A 863"`` (altloc in parentheses)."""
  rg = ag.parent()
  alt = ag.altloc.strip()
  return (f"{ag.resname.strip()} {rg.parent().id.strip()} {rg.resseq.strip()}"
          + (f" ({alt})" if alt else ""))


def _worst_water_clashes(hier, crystal_symmetry=None,
                         min_distance_sym_equiv=_WATER_SYM_EQUIV_TOL):
  """Inter-water H-H contacts within 2.0 A, closest first.

  One ``(distance, id_a, id_b)`` per contact, as counted by
  :func:`_water_clash_stats`; ``id_b`` of a symmetry equivalent ends in its
  operator.
  """
  atoms, contacts = _all_water_h_contacts(hier, crystal_symmetry,
                                          min_distance_sym_equiv)
  return [(d, _atom_id(atoms[i]),
           _atom_id(atoms[j]) + (f" ({op})" if op is not None else ""))
          for d, i, j, op in sorted(contacts, key=lambda c: c[0])]


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

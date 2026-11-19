"""H on a parent atom that has no heavy neighbour (NH4+, H3O+, H2S, CH4).

Riding H needs a heavy neighbour of the parent to define a frame, so
reduce2's placement (mmtbx.hydrogens.reduce_hydrogen) dropped these H and
reported them as "riding H refused"; an ammonium ion then left the model as
a bare N, which the charge calculation counts as N3- instead of NH4+.

Here the H of one parent are a rigid ideal group (bond lengths and the
H-X-H angle from the restraints: a tetrahedron for four, a pyramid for
three, a bent pair for two) whose orientation is searched exhaustively
(golden-spiral directions times azimuths) for the best score against the
environment within ENVIRONMENT_RADIUS of the parent, symmetry mates
included: one H-bond gain per H for an acceptor (O, S, halide, or an N
that carries no H; cations are not acceptors) at H...A within HBOND_RANGE,
growing linearly with the X-H...A angle from HBOND_MIN_ANGLE (0) to
straight (1), minus a clearance penalty for every atom closer than
MIN_CLEARANCE (heavy) or MIN_H_CLEARANCE (H). The H stay free atoms (not
riding): the restraints keep their geometry.
"""
from __future__ import absolute_import, division, print_function
import math
from libtbx import group_args
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.math import golden_spiral

ENVIRONMENT_RADIUS = 5.0
MIN_CLEARANCE = 1.5     # A, to a heavy atom
MIN_H_CLEARANCE = 2.0   # A, to an H
CLEARANCE_WEIGHT = 10.0 # score per A of shortfall
HBOND_RANGE = (1.5, 2.7)
HBOND_MIN_ANGLE = 110.0
NH_BOND = 1.3           # an N with an H this close is a donor, not an acceptor
N_DIRECTIONS = 120
N_AZIMUTHS = 24
ACCEPTOR_ELEMENTS = frozenset(("O", "S", "F", "CL", "BR", "I", "N"))
TETRAHEDRAL = 109.47

def ideal_directions(n, angle):
  """n unit vectors around the origin with the pairwise angle (deg): one
  along z; two in the xz plane; three on a cone around z (a pyramid);
  four a tetrahedron (angle ignored); more spread on a golden spiral."""
  if n == 1: return [matrix.col((0, 0, 1))]
  a = math.radians(angle)
  if n == 2:
    s, c = math.sin(a / 2), math.cos(a / 2)
    return [matrix.col((s, 0, c)), matrix.col((-s, 0, c))]
  if n == 3:
    # cos(angle) = 1 - 1.5 sin^2(alpha) for three vectors at polar alpha
    s2 = min(1.0, max(0.0, (1 - math.cos(a)) / 1.5))
    alpha = math.asin(math.sqrt(s2))
    return [matrix.col((math.sin(alpha) * math.cos(phi),
                        math.sin(alpha) * math.sin(phi), math.cos(alpha)))
            for phi in (0, 2 * math.pi / 3, 4 * math.pi / 3)]
  if n == 4:
    r = 1 / math.sqrt(3)
    return [matrix.col(v) * r for v in ((1, 1, 1), (1, -1, -1),
                                        (-1, 1, -1), (-1, -1, 1))]
  return [matrix.col(p) for p in list(golden_spiral(n).points)]

def environment(model, i_parent, h_i_seqs, ignore=None,
                radius=ENVIRONMENT_RADIUS):
  """
  The atoms within radius of atom i_parent of model, symmetry mates
  included: group_args(xyz: flex.vec3_double, clearance: flex.double per
  atom, acceptor: list of indices into xyz). h_i_seqs (the H being placed)
  are left out, so are the atoms selected by ignore (flex.bool).
  """
  hierarchy = model.get_hierarchy()
  atoms = hierarchy.atoms()
  sites = atoms.extract_xyz()
  p = sites[i_parent]
  skip = set(h_i_seqs)
  xyz, elements = [], []
  def take(atom, site, i_seq):
    if i_seq in skip or i_seq == i_parent: return
    if ignore is not None and ignore[i_seq]: return
    xyz.append(site)
    elements.append(atom.element.strip().upper())
  d = (sites - flex.vec3_double(sites.size(), p)).norms()
  for i_seq in (d < radius).iselection():
    take(atoms[i_seq], sites[i_seq], i_seq)
  cs = model.crystal_symmetry()
  if cs is not None and cs.unit_cell() is not None and \
     cs.space_group() is not None:
    for equiv in hierarchy.sym_equiv_residue_groups(cs, radius,
        selection=flex.size_t([i_parent]), min_distance_sym_equiv=0.5):
      sub = equiv.hierarchy.atoms()
      for atom, site, src in zip(sub, sub.extract_xyz(), equiv.i_seqs):
        if src in skip: continue
        if ignore is not None and ignore[src]: continue
        if (matrix.col(site) - matrix.col(p)).length() >= radius: continue
        xyz.append(site)
        elements.append(atom.element.strip().upper())
  xyz = flex.vec3_double(xyz)
  is_h = [e in ("H", "D") for e in elements]
  clearance = flex.double([MIN_H_CLEARANCE if h else MIN_CLEARANCE
                           for h in is_h])
  from mmtbx.hydrogens.water_protonation import _WATER_CATION_ELEMENTS
  acceptor = []
  for i, e in enumerate(elements):
    if e not in ACCEPTOR_ELEMENTS or e in _WATER_CATION_ELEMENTS: continue
    if e == "N":
      dn = (xyz - flex.vec3_double(xyz.size(), xyz[i])).norms()
      if any(is_h[j] and dn[j] < NH_BOND for j in range(xyz.size())):
        continue
    acceptor.append(i)
  return group_args(xyz=xyz, clearance=clearance, acceptor=acceptor)

def _rotations():
  for d in list(golden_spiral(N_DIRECTIONS).points):
    d = matrix.col(d)
    e1 = d.ortho().normalize()
    e2 = d.cross(e1)
    for k in range(N_AZIMUTHS):
      phi = 2 * math.pi * k / N_AZIMUTHS
      yield (e1 * math.cos(phi) + e2 * math.sin(phi),
             -e1 * math.sin(phi) + e2 * math.cos(phi), d)

def score(parent, sites, env):
  """Score of the H at sites (list of matrix.col) on parent against env
  (environment()): H-bond gains minus clearance penalties."""
  if env.xyz.size() == 0: return 0.0
  total = 0.0
  p = matrix.col(parent)
  for h in sites:
    dv = env.xyz - flex.vec3_double(env.xyz.size(), h.elems)
    d = dv.norms()
    short = env.clearance - d
    total -= CLEARANCE_WEIGHT * flex.sum(short.select(short > 0))
    best = 0.0
    hp = p - h
    for i in env.acceptor:
      if not HBOND_RANGE[0] <= d[i] <= HBOND_RANGE[1]: continue
      ha = matrix.col(dv[i])
      cos_angle = hp.dot(ha) / (hp.length() * d[i])
      angle = math.degrees(math.acos(max(-1.0, min(1.0, cos_angle))))
      gain = (angle - HBOND_MIN_ANGLE) / (180.0 - HBOND_MIN_ANGLE)
      if gain > best: best = gain
    total += best
  return total

def place(parent_xyz, bond_lengths, angle, environment):
  """
  Sites (tuples) of len(bond_lengths) H on the parent at parent_xyz, at
  the bond lengths and the H-X-H angle (deg), in the orientation of the best
  score against environment (environment()); the first orientation on a
  tie. One H points along the first direction tried, so even an isolated
  group is placed deterministically.
  """
  local = ideal_directions(len(bond_lengths), angle)
  p = matrix.col(parent_xyz)
  best, best_sites = None, None
  for e1, e2, e3 in _rotations():
    sites = [p + (e1 * v[0] + e2 * v[1] + e3 * v[2]) * length
             for v, length in zip(local, bond_lengths)]
    s = score(p, sites, environment)
    if best is None or s > best + 1e-9:
      best, best_sites = s, sites
  return [tuple(s) for s in best_sites]

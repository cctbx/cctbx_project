"""
Side-chain geometry for rotamer_bayes: the chi angles of a residue (atoms as
mmtbx/rotamer/sidechain_angles.props), the atoms each chi moves (those on the
far side of its axis bond, found from the residue's own covalent graph, so
riding hydrogens move with their parents), setting the chi angles, and the
residue's environment under the crystal symmetry for clash checks.
"""
from __future__ import absolute_import, division, print_function
import math
import numpy as np

def dihedral(p):
  """ Dihedral angle (degrees) of four points, as scitbx.math.dihedral_angle """
  b0 = p[0] - p[1]
  b1 = p[2] - p[1]
  b2 = p[3] - p[2]
  b1n = b1 / np.linalg.norm(b1)
  v = b0 - np.dot(b0, b1n) * b1n
  w = b2 - np.dot(b2, b1n) * b1n
  x = np.dot(v, w)
  y = np.dot(np.cross(b1n, v), w)
  return math.degrees(math.atan2(y, x))

def dihedral_gradients(p):
  """ d(dihedral, degrees)/d(points) (4, 3), for dihedral() above
  (Blondel and Karplus, 1996) """
  f = p[0] - p[1]; g = p[1] - p[2]; h = p[3] - p[2]
  a = np.cross(f, g); b = np.cross(h, g)
  aa = a.dot(a); bb = b.dot(b); gl = math.sqrt(g.dot(g))
  d0 = -gl / aa * a
  d3 = gl / bb * b
  d1 = gl / aa * a + f.dot(g) / (aa * gl) * a - h.dot(g) / (bb * gl) * b
  d2 = -d0 - d1 - d3
  return np.degrees(np.array([d0, d1, d2, d3]))

def rotate(xyz, origin, axis, angle_deg):
  """ Rotate points about the axis through origin (right-handed, degrees) """
  k = axis / np.linalg.norm(axis)
  a = math.radians(angle_deg)
  v = xyz - origin
  c, s = math.cos(a), math.sin(a)
  return origin + v * c + np.cross(k, v) * s + np.outer(v.dot(k), k) * (1 - c)

def rotation_matrix(axis, angle_deg):
  k = axis / np.linalg.norm(axis)
  a = math.radians(angle_deg)
  c, s = math.cos(a), math.sin(a)
  x, y, z = k
  return np.array([
    [c + x*x*(1-c), x*y*(1-c) - z*s, x*z*(1-c) + y*s],
    [y*x*(1-c) + z*s, c + y*y*(1-c), y*z*(1-c) - x*s],
    [z*x*(1-c) - y*s, z*y*(1-c) + x*s, c + z*z*(1-c)]])

def _bonded(names, elements, xyz):
  is_h = np.array([e in ("H", "D") for e in elements])
  d = np.sqrt(((xyz[:, None, :] - xyz[None, :, :]) ** 2).sum(axis=2))
  cut = np.where(is_h[:, None] | is_h[None, :], 1.3, 2.0)
  b = (d < cut) & (d > 0.1)
  return b

class side_chain(object):
  """ One residue (single conformer): chi definitions and movers.

  ok is False (with reason) when the residue cannot be handled: no chi
  definition, missing chi atoms, or a cyclic side chain (PRO). """

  def __init__(self, residue, sidechain_angles):
    self.residue = residue
    self.resname = residue.resname.strip().upper()
    atoms = residue.atoms()
    self.i_seqs = np.array(list(atoms.extract_i_seq()), dtype=int)
    self.names = [a.name.strip() for a in atoms]
    self.elements = [a.element.strip().upper() for a in atoms]
    self.xyz = np.array(list(atoms.extract_xyz()))
    self.is_h = np.array([e in ("H", "D") for e in self.elements])
    self.ok = False
    self.reason = None
    aa = self.resname.lower()
    if(aa == "mse"): aa = "met"
    if(aa == "pro"):
      self.reason = "proline"; return
    n_chi = sidechain_angles.chisPerAA.get(aa)
    if(n_chi is None):
      self.reason = "no chi definition"; return
    index = dict((n, i) for i, n in enumerate(self.names))
    self.chi_atoms = []
    for k in range(int(n_chi)):
      key = "%s.chi%d" % (aa, k + 1)
      quad = []
      for alternatives in sidechain_angles.atomsForAngle[key]:
        found = None
        for name in alternatives.split(";"):
          name = name.strip()
          if(self.resname == "MSE" and name == "SD"): name = "SE"
          if(name in index):
            found = index[name]; break
        if(found is None):
          self.reason = "missing atom for %s" % key; return
        quad.append(found)
      self.chi_atoms.append(quad)
    bonds = _bonded(self.names, self.elements, self.xyz)
    self.movers = []
    for quad in self.chi_atoms:
      j, k = quad[1], quad[2]
      # atoms reached from k without crossing the bond j-k
      seen = set([k])
      stack = [k]
      while stack:
        a = stack.pop()
        for b in np.nonzero(bonds[a])[0]:
          b = int(b)
          if((a == k and b == j) or b in seen): continue
          seen.add(b); stack.append(b)
      if(j in seen):
        self.reason = "cyclic side chain"; return
      seen.discard(k)
      mask = np.zeros(len(self.names), dtype=bool)
      mask[list(seen)] = True
      self.movers.append(mask)
    # the atoms any chi moves: the ones whose density changes
    self.moving = np.zeros(len(self.names), dtype=bool)
    for m in self.movers: self.moving |= m
    self.ok = len(self.movers) > 0
    if(not self.ok): self.reason = "no movable atoms"

  def chis(self, xyz=None):
    if(xyz is None): xyz = self.xyz
    return [dihedral(xyz[q]) for q in self.chi_atoms]

  def set_chis(self, chis, xyz=None):
    """ Coordinates with the chi angles set to chis (degrees; None: keep) """
    if(xyz is None): xyz = self.xyz
    xyz = xyz.copy()
    for quad, mask, target in zip(self.chi_atoms, self.movers, chis):
      if(target is None): continue
      delta = target - dihedral(xyz[quad])
      delta = (delta + 180.) % 360. - 180.
      origin = xyz[quad[2]]
      axis = xyz[quad[2]] - xyz[quad[1]]
      xyz[mask] = rotate(xyz[mask], origin, axis, delta)
    return xyz

  def rotate_chi(self, xyz, k, delta):
    """ Coordinates with chi k changed by delta (degrees) """
    quad, mask = self.chi_atoms[k], self.movers[k]
    origin = xyz[quad[2]]
    m = rotation_matrix(xyz[quad[2]] - xyz[quad[1]], delta)
    out = xyz.copy()
    out[mask] = (xyz[mask] - origin).dot(m.T) + origin
    return out

class environment(object):
  """ Heavy atoms of the structure near a residue, under the crystal
  symmetry (nearest lattice translation of every symmetry copy), excluding
  the residue's own atoms in the identity copy. """

  def __init__(self, crystal_symmetry, sites_cart, keep):
    uc = crystal_symmetry.unit_cell()
    self.f = np.array(uc.fractionalization_matrix()).reshape(3, 3)
    self.o = np.array(uc.orthogonalization_matrix()).reshape(3, 3)
    self.ops = []
    for op in crystal_symmetry.space_group().all_ops():
      self.ops.append((np.array(op.r().as_double()).reshape(3, 3),
                       np.array(op.t().as_double())))
    self.keep = np.nonzero(keep)[0]
    self.frac = np.array(sites_cart)[self.keep].dot(self.f.T)

  def near(self, centre_cart, radius, exclude_i_seqs, return_i_seqs=False):
    """ Cartesian sites within radius of centre (identity copies of
    exclude_i_seqs left out); with return_i_seqs also their i_seqs """
    c = self.f.dot(centre_cart)
    out, idx = [], []
    excl = np.isin(self.keep, exclude_i_seqs)
    for n, (r, t) in enumerate(self.ops):
      d = self.frac.dot(r.T) + t - c
      d -= np.round(d)
      x = d.dot(self.o.T)
      sel = (x * x).sum(axis=1) <= radius * radius
      if(n == 0): sel &= ~excl
      out.append(x[sel] + centre_cart)
      idx.append(self.keep[sel])
    xyz = np.concatenate(out) if out else np.zeros((0, 3))
    if(return_i_seqs):
      return xyz, (np.concatenate(idx) if idx else np.zeros(0, dtype=int))
    return xyz

def _turn_to(p, origin, axis, target):
  """ angle (degrees) about axis through origin that brings p closest to
  target """
  k = axis / np.linalg.norm(axis)
  v = p - origin; w = target - origin
  v = v - v.dot(k) * k; w = w - w.dot(k) * k
  return math.degrees(math.atan2(np.cross(k, v).dot(w), v.dot(w)))

def backrub(prev, res, nxt, theta):
  """ Backrub of residue i by theta (degrees): rotation of C, O (and H of
  N) of the peptides on both sides and of the residue about the
  CA(i-1)-CA(i+1) axis, then counter-rotations of the two peptides (about
  CA(i-1)-CA(i) and CA(i)-CA(i+1)) that bring both carbonyl O back as close
  as possible to where they were (Davis et al., Structure 14, 265, 2006).
  prev, res, nxt: dicts atom name -> site (numpy) of residues i-1, i, i+1;
  returns new dicts for the atoms that move (prev C, O; all of res; nxt N,
  H). """
  ca_p, ca_n = prev["CA"], nxt["CA"]
  axis = ca_n - ca_p
  m = rotation_matrix(axis, theta)
  turn = lambda x: (x - ca_p).dot(m.T) + ca_p
  p1 = dict((n, turn(prev[n])) for n in ("C", "O"))
  r1 = dict((n, turn(x)) for n, x in res.items())
  n1 = dict((n, turn(nxt[n])) for n in ("N", "H") if n in nxt)
  # peptide i-1: C(i-1), O(i-1), N(i) and its H, about CA(i-1)-CA(i)
  ca_i = r1["CA"]
  t1 = _turn_to(p1["O"], ca_p, ca_i - ca_p, prev["O"])
  m1 = rotation_matrix(ca_i - ca_p, t1)
  for d, names in ((p1, ("C", "O")), (r1, ("N", "H"))):
    for n in names:
      if(n in d): d[n] = (d[n] - ca_p).dot(m1.T) + ca_p
  # peptide i: C(i), O(i), N(i+1) and its H, about CA(i)-CA(i+1)
  t2 = _turn_to(r1["O"], ca_i, ca_n - ca_i, res["O"])
  m2 = rotation_matrix(ca_n - ca_i, t2)
  for d, names in ((r1, ("C", "O")), (n1, ("N", "H"))):
    for n in names:
      if(n in d): d[n] = (d[n] - ca_i).dot(m2.T) + ca_i
  return p1, r1, n1

from cctbx import geometry_restraints
import scitbx.math
from scitbx import matrix
import libtbx.phil
from libtbx.math_utils import normalize_angle

master_phil = libtbx.phil.parse("""\
  bond_min_distance = 1.2
    .type = float
  bond_max_distance = 1.8
    .type = float
  epsilon_range_not_2p_min = 155
    .type = float
  epsilon_range_not_2p_max = 310
    .type = float
  delta_range_2p_min = 115
    .type = float
  delta_range_2p_max = 180
    .type = float
  p_distance_c1_n_line_2p_max = 2.9
    .type = float
""")

c1p_n_name = {
  "A": " N9 ",
  "G": " N9 ",
  "C": " N1 ",
  "U": " N1 "}

def _build_required_atom_indices():
  result = {}
  # attention: indices hard-wired in evaluate.__init__ for efficiency
  common = ["P", "C1'", "O2'", "O3'", "C3'", "C4'", "C5'"]
  def as_dict(l): return dict(zip(l, range(len(l))))
  result["A"] = as_dict(common + [c1p_n_name["A"].strip()])
  result["G"] = result["A"]
  result["C"] = as_dict(common + [c1p_n_name["C"].strip()])
  result["U"] = result["C"]
  return result

required_atom_indices = _build_required_atom_indices()

def extract_required_atoms(resname, atoms, atom_names=None):
  assert atom_names is None or len(atom_names) == len(atoms)
  if (atom_names is None):
    atom_names = [atom.name for atom in atoms]
  indices = required_atom_indices.get(resname.strip())
  if (indices is None):
    raise RuntimeError('Unexpected residue name: "%s"' % resname)
  result = [None] * len(indices)
  n_found = 0
  for atom, atom_name in zip(atoms, atom_names):
    if (atom_name is None): continue
    i = indices.get(atom_name.strip())
    if (i is None): continue
    if (result[i] is not None):
      raise RuntimeError('Duplicate atom name: "%s"' % atom_name)
    result[i] = atom
    n_found += 1
  if (n_found != len(result)):
    result = None
  return result

class evaluate(object):

  def given_atoms(params, atoms_1, atoms_2):
    return evaluate(
      params=params,
      sites_cart_1=[atom.xyz for atom in atoms_1],
      sites_cart_2=[atom.xyz for atom in atoms_2])
  given_atoms = staticmethod(given_atoms)

  def __init__(self, params, sites_cart_1, sites_cart_2):
    assert len(sites_cart_1) == 8
    assert len(sites_cart_2) == 8
    p = matrix.col(sites_cart_2[0])
    def v(i): return matrix.col(sites_cart_1[i])
    c1p = v(1)
    o3p = v(3)
    c3p = v(4)
    c4p = v(5)
    c5p = v(6)
    n = v(7)
    bonded_sites = [
      (p, o3p),
      (o3p, c3p),
      (c3p, c4p),
      (c1p, n),
      (c5p, c4p)]
    distances = [abs(s2-s1) for s1,s2 in bonded_sites]
    if (   min(distances) < params.bond_min_distance
        or max(distances) > params.bond_max_distance):
      epsilon = None
      delta = None
      p_distance_c1_n_line = None
      is_2p_epsilon = None
      is_2p_delta = None
      is_2p_p_distance_c1_n_line = None
      is_2p = None
    else:
      def a(sites): return normalize_angle(geometry_restraints.dihedral(
        sites=sites, angle_ideal=0, weight=1).angle_model, deg=True)
      epsilon = a([p, o3p, c3p, c4p])
      delta = a([c5p, c4p, c3p, o3p])
      p_distance_c1_n_line = scitbx.math.line_given_points(
        points=(c1p, n)).distance_sq(point=p)**0.5
      n_decisions = 0
      if (   params.epsilon_range_not_2p_min is None
          or params.epsilon_range_not_2p_max is None):
        is_2p_epsilon = None
      else:
        is_2p_epsilon = (   epsilon < params.epsilon_range_not_2p_min
                         or epsilon > params.epsilon_range_not_2p_max)
        n_decisions += 1
      if (   params.delta_range_2p_min is None
          or params.delta_range_2p_max is None):
        is_2p_delta = None
      else:
        is_2p_delta = (    delta >= params.delta_range_2p_min
                       and delta <= params.delta_range_2p_max)
        n_decisions += 1
      if (params.p_distance_c1_n_line_2p_max is None):
        is_2p_p_distance_c1_n_line = None
      else:
        is_2p_p_distance_c1_n_line = (
          p_distance_c1_n_line < params.p_distance_c1_n_line_2p_max)
        n_decisions += 1
      if (n_decisions == 0):
        is_2p = None
      else:
        is_2p = (   is_2p_epsilon
                 or is_2p_delta
                 or is_2p_p_distance_c1_n_line)
    self.epsilon = epsilon
    self.delta = delta
    self.p_distance_c1_n_line = p_distance_c1_n_line
    self.is_2p_epsilon = is_2p_epsilon
    self.is_2p_delta = is_2p_delta
    self.is_2p_p_distance_c1_n_line = is_2p_p_distance_c1_n_line
    self.is_2p = is_2p

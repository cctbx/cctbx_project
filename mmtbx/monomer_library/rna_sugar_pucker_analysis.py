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
  p_distance_c1p_outbound_line_2p_max = 2.9
    .type = float
  bond_detection_distance_tolerance = 0.5
    .type = float
""")

class evaluate(object):

  def __init__(self,
        params,
        residue_1_deoxy_ribo_atom_dict,
        residue_1_c1p_outbound_atom,
        residue_2_p_atom):
    if (   residue_1_c1p_outbound_atom is None
        or residue_2_p_atom is None):
      distances = None
    else:
      p = matrix.col(residue_2_p_atom.xyz)
      def v(key): return matrix.col(residue_1_deoxy_ribo_atom_dict[key].xyz)
      c1p = v("C1'")
      o3p = v("O3'")
      c3p = v("C3'")
      c4p = v("C4'")
      c5p = v("C5'")
      c1p_outbound = matrix.col(residue_1_c1p_outbound_atom.xyz)
      bonded_sites = [
        (p, o3p),
        (o3p, c3p),
        (c3p, c4p),
        (c1p, c1p_outbound),
        (c5p, c4p)]
      distances = [abs(s2-s1) for s1,s2 in bonded_sites]
    if (   distances is None
        or min(distances) < params.bond_min_distance
        or max(distances) > params.bond_max_distance):
      epsilon = None
      delta = None
      p_distance_c1p_outbound_line = None
      is_2p_epsilon = None
      is_2p_delta = None
      is_2p_p_distance_c1p_outbound_line = None
      is_2p = None
    else:
      def a(sites): return normalize_angle(geometry_restraints.dihedral(
        sites=sites, angle_ideal=0, weight=1).angle_model, deg=True)
      epsilon = a([p, o3p, c3p, c4p])
      delta = a([c5p, c4p, c3p, o3p])
      p_distance_c1p_outbound_line = scitbx.math.line_given_points(
        points=(c1p, c1p_outbound)).distance_sq(point=p)**0.5
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
      if (params.p_distance_c1p_outbound_line_2p_max is None):
        is_2p_p_distance_c1p_outbound_line = None
      else:
        is_2p_p_distance_c1p_outbound_line = (
            p_distance_c1p_outbound_line
          < params.p_distance_c1p_outbound_line_2p_max)
        n_decisions += 1
      if (n_decisions == 0):
        is_2p = None
      else:
        is_2p = (   is_2p_epsilon
                 or is_2p_delta
                 or is_2p_p_distance_c1p_outbound_line)
    self.epsilon = epsilon
    self.delta = delta
    self.p_distance_c1p_outbound_line = p_distance_c1p_outbound_line
    self.is_2p_epsilon = is_2p_epsilon
    self.is_2p_delta = is_2p_delta
    self.is_2p_p_distance_c1p_outbound_line =is_2p_p_distance_c1p_outbound_line
    self.is_2p = is_2p

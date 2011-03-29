from cctbx import crystal, sgtbx
from iotbx.cif import model
from libtbx.utils import format_float_with_standard_uncertainty \
     as format_float_with_su


class distances_as_cif_loop(object):

  def __init__(self,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               eps=2e-16):
    assert [sites_frac, sites_cart].count(None) == 1
    fmt = "%.4f"
    asu_mappings = pair_asu_table.asu_mappings()
    space_group_info = sgtbx.space_group_info(group=asu_mappings.space_group())
    unit_cell = asu_mappings.unit_cell()
    if sites_cart is not None:
      sites_frac = unit_cell.fractionalize(sites_cart)
    self.loop = model.loop(header=(
      "_geom_bond_atom_site_label_1",
      "_geom_bond_atom_site_label_2",
      "_geom_bond_distance",
      "_geom_bond_site_symmetry_2"
    ))
    distances = crystal.calculate_distances(
      pair_asu_table, sites_frac,
      covariance_matrix=covariance_matrix,
      cell_covariance_matrix=cell_covariance_matrix,
      parameter_map=parameter_map)
    for d in distances:
      if site_labels[d.i_seq].startswith('H') or site_labels[d.j_seq].startswith('H'):
        continue
      if d.variance is not None and d.variance > eps:
        distance = format_float_with_su(d.distance, math.sqrt(d.variance))
      else:
        distance = fmt % d.distance
      sym_code = space_group_info.cif_symmetry_code(d.rt_mx_ji)
      if sym_code == "1": sym_code = "."
      self.loop.add_row((site_labels[d.i_seq],
                         site_labels[d.j_seq],
                         distance,
                         sym_code))
    self.distances = distances.distances
    self.variances = distances.variances
    self.pair_counts = distances.pair_counts

class angles_as_cif_loop(object):

  def __init__(self,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               eps=2e-16):
    assert [sites_frac, sites_cart].count(None) == 1
    fmt = "%.1f"
    asu_mappings = pair_asu_table.asu_mappings()
    space_group_info = sgtbx.space_group_info(group=asu_mappings.space_group())
    unit_cell = asu_mappings.unit_cell()
    if sites_cart is not None:
      sites_frac = unit_cell.fractionalize(sites_cart)
    self.loop = model.loop(header=(
      "_geom_angle_atom_site_label_1",
      "_geom_angle_atom_site_label_2",
      "_geom_angle_atom_site_label_3",
      "_geom_angle",
      "_geom_angle_site_symmetry_1",
      "_geom_angle_site_symmetry_3"
    ))
    angles = crystal.calculate_angles(
      pair_asu_table, sites_frac,
      covariance_matrix=covariance_matrix,
      cell_covariance_matrix=cell_covariance_matrix,
      parameter_map=parameter_map)
    for a in angles:
      i_seq, j_seq, k_seq = a.i_seqs
      if site_labels[i_seq].startswith('H') or site_labels[k_seq].startswith('H'):
        continue
      sym_code_ji = space_group_info.cif_symmetry_code(a.rt_mx_ji)
      sym_code_ki = space_group_info.cif_symmetry_code(a.rt_mx_ki)
      if sym_code_ji == "1": sym_code_ji = "."
      if sym_code_ki == "1": sym_code_ki = "."
      if a.variance is not None and a.variance > eps:
        angle = format_float_with_su(a.angle, math.sqrt(a.variance))
      else:
        angle = fmt % a.angle
      self.loop.add_row((site_labels[i_seq],
                         site_labels[j_seq],
                         site_labels[k_seq],
                         angle,
                         sym_code_ji,
                         sym_code_ki,
                         ))
    self.angles = angles.angles
    self.variances = angles.variances


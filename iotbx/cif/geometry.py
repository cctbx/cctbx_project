from cctbx.array_family import flex
from cctbx import crystal, sgtbx
from cctbx import covariance, geometry
from iotbx.cif import model
from libtbx.utils import format_float_with_standard_uncertainty \
     as format_float_with_su
from libtbx import adopt_init_args

import math


class distances_as_cif_loop(object):

  def __init__(self,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               include_bonds_to_hydrogen=False,
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
      if (not include_bonds_to_hydrogen
          and (site_labels[d.i_seq].startswith('H') or
               site_labels[d.j_seq].startswith('H'))):
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
               include_bonds_to_hydrogen=False,
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
      if (not include_bonds_to_hydrogen
          and (site_labels[i_seq].startswith('H') or
               site_labels[k_seq].startswith('H'))):
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


class hbond(object):
  def __init__(self, d_seq, a_seq, rt_mx=None):
    # rt_mx is the optional symmetry operator for the acceptor atom
    adopt_init_args(self, locals())


unit_mx = sgtbx.rt_mx()

class hbonds_as_cif_loop(object):

  def __init__(self,
               hbonds,
               pair_asu_table,
               site_labels,
               sites_frac=None,
               sites_cart=None,
               min_dha_angle=150, # degrees
               max_da_distance=2.9, # angstrom
               covariance_matrix=None,
               cell_covariance_matrix=None,
               parameter_map=None,
               eps=2e-16):
    assert [sites_frac, sites_cart].count(None) == 1
    fmt_a = "%.1f"
    pair_asu_table = pair_asu_table
    asu_mappings = pair_asu_table.asu_mappings()
    space_group_info = sgtbx.space_group_info(group=asu_mappings.space_group())
    self.unit_cell = asu_mappings.unit_cell()
    if sites_cart is not None:
      sites_frac = self.unit_cell.fractionalize(sites_cart)
    if sites_frac is not None:
      sites_cart = self.unit_cell.orthogonalize(sites_frac)
    if covariance_matrix is not None:
      assert parameter_map is not None
      self.covariance_matrix_cart = covariance.orthogonalize_covariance_matrix(
        covariance_matrix, self.unit_cell, parameter_map)
    else: self.covariance_matrix_cart = None
    self.cell_covariance_matrix = cell_covariance_matrix
    self.eps = eps
    self.loop = model.loop(header=(
      "_geom_hbond_atom_site_label_D",
      "_geom_hbond_atom_site_label_H",
      "_geom_hbond_atom_site_label_A",
      "_geom_hbond_distance_DH",
      "_geom_hbond_distance_HA",
      "_geom_hbond_distance_DA",
      "_geom_hbond_angle_DHA",
      "_geom_hbond_site_symmetry_A",
    ))
    for hbond in hbonds:
      d_seq, a_seq = hbond.d_seq, hbond.a_seq
      site_cart_d = sites_cart[d_seq]
      if hbond.rt_mx is not None:
        site_frac_a = sites_frac[a_seq]
        site_frac_a = hbond.rt_mx * site_frac_a
        site_cart_a = self.unit_cell.orthogonalize(site_frac_a)
      else:
        site_cart_a = sites_cart[a_seq]
      distance_da = geometry.distance((site_cart_d, site_cart_a))
      for h_seq, h_sym_groups in pair_asu_table.table()[hbond.d_seq].items():
        if site_labels[h_seq][0] not in ('H','D'):
          # XXX better to pass scattering types instead?
          continue
        site_cart_h = sites_cart[h_seq]
        distance_dh = geometry.distance((site_cart_d, site_cart_h))
        distance_ha = geometry.distance((site_cart_h, site_cart_a))
        angle_dha = geometry.angle((site_cart_d, site_cart_h, site_cart_a))
        if (angle_dha.angle_model < min_dha_angle or
            distance_da.distance_model > max_da_distance):
          continue
        if hbond.rt_mx is not None:
          sym_code = space_group_info.cif_symmetry_code(hbond.rt_mx)
        else: sym_code = '.'
        self.loop.add_row((
          site_labels[d_seq],
          site_labels[h_seq],
          site_labels[a_seq],
          self.formatted_distance(d_seq, h_seq, distance_dh, unit_mx),
          self.formatted_distance(h_seq, a_seq, distance_ha, unit_mx),
          self.formatted_distance(d_seq, a_seq, distance_da, hbond.rt_mx),
          self.formatted_angle(d_seq, h_seq, a_seq, angle_dha, hbond.rt_mx),
          sym_code
        ))

  def formatted_distance(self, i_seq, j_seq, distance, rt_mx_ji):
    if rt_mx_ji is None: rt_mx_ji = unit_mx
    if self.covariance_matrix_cart is not None:
      cov = covariance.extract_covariance_matrix_for_sites(
        flex.size_t((i_seq,j_seq)),
        self.covariance_matrix_cart,
        self.parameter_map)
      if self.cell_covariance_matrix is not None:
        var = distance.variance(
          cov, self.cell_covariance_matrix, self.unit_cell, rt_mx_ji)
      else:
        var = distance.variance(cov, self.unit_cell, rt_mx_ji)
      if var > self.eps:
        return format_float_with_su(distance.distance_model, math.sqrt(var))
    return "%.4f" %distance.distance_model

  def formatted_angle(self, i_seq, j_seq, k_seq, angle, rt_mx_ki):
    if rt_mx_ki is None: rt_mx_ki = unit_mx
    if self.covariance_matrix_cart is not None:
      cov = covariance.extract_covariance_matrix_for_sites(
        flex.size_t((i_seq,j_seq,k_seq)),
        self.covariance_matrix_cart,
        self.parameter_map)
      if self.cell_covariance_matrix is not None:
        var = angle.variance(
          cov, self.cell_covariance_matrix, self.unit_cell, rt_mx_ki)
      else:
        var = angle.variance(cov, self.unit_cell, rt_mx_ji)
      if var > self.eps:
        return format_float_with_su(angle.distance_model, math.sqrt(var))
    return "%.1f" %angle.angle_model

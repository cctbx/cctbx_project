from __future__ import division
from cctbx import sgtbx
from cctbx import adp_restraints, geometry_restraints
from cctbx.adp_restraints import adp_restraint_params
from cctbx.array_family import flex
from iotbx.cif import model
import math

# http://www.iucr.org/__data/iucr/cifdic_html/1/cif_core_restraints.dic/index.html

def add_to_cif_block(cif_block, xray_structure,
                     bond_proxies=None,
                     angle_proxies=None,
                     dihedral_proxies=None,
                     bond_similarity_proxies=None,
                     rigid_bond_proxies=None,
                     adp_similarity_proxies=None,
                     isotropic_adp_proxies=None):
  if bond_proxies is not None:
    cif_block.add_loop(distances_as_cif_loop(xray_structure, bond_proxies))
  if angle_proxies is not None:
    cif_block.add_loop(angles_as_cif_loop(xray_structure, angle_proxies))
  if dihedral_proxies is not None:
    cif_block.add_loop(dihedrals_as_cif_loop(xray_structure, dihedral_proxies))
  if bond_similarity_proxies is not None:
    loops = bond_similarity_as_cif_loops(xray_structure, bond_similarity_proxies)
    for loop in loops: cif_block.add_loop(loop)
  if rigid_bond_proxies is not None:
    cif_block.add_loop(rigid_bond_as_cif_loop(xray_structure, rigid_bond_proxies))
  if adp_similarity_proxies is not None:
    cif_block.add_loop(
      adp_similarity_as_cif_loop(xray_structure, adp_similarity_proxies))
  if isotropic_adp_proxies is not None:
    cif_block.add_loop(
      isotropic_adp_as_cif_loop(xray_structure, isotropic_adp_proxies))

def distances_as_cif_loop(xray_structure, proxies):
  space_group_info = sgtbx.space_group_info(group=xray_structure.space_group())
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_distance_atom_site_label_1",
    "_restr_distance_atom_site_label_2",
    "_restr_distance_site_symmetry_2",
    "_restr_distance_target",
    "_restr_distance_target_weight_param",
    "_restr_distance_diff"
  ))
  for proxy in proxies:
    restraint = geometry_restraints.bond(
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      proxy=proxy)
    i_seqs = proxy.i_seqs
    sym_op = proxy.rt_mx_ji
    if sym_op is None: sym_op = sgtbx.rt_mx()
    loop.add_row((site_labels[i_seqs[0]],
                  site_labels[i_seqs[1]],
                  space_group_info.cif_symmetry_code(sym_op),
                  fmt % restraint.distance_ideal,
                  fmt % math.sqrt(1/restraint.weight),
                  fmt % restraint.delta))
  return loop

def angles_as_cif_loop(xray_structure, proxies):
  space_group_info = sgtbx.space_group_info(group=xray_structure.space_group())
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_angle_atom_site_label_1",
    "_restr_angle_atom_site_label_2",
    "_restr_angle_atom_site_label_3",
    "_restr_angle_site_symmetry_1",
    "_restr_angle_site_symmetry_2",
    "_restr_angle_site_symmetry_3",
    "_restr_angle_target",
    "_restr_angle_target_weight_param",
    "_restr_angle_diff",
  ))
  unit_mxs = [sgtbx.rt_mx()]*3
  for proxy in proxies:
    restraint = geometry_restraints.angle(
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      proxy=proxy)
    sym_ops = proxy.sym_ops
    if sym_ops is None: sym_ops = unit_mxs
    i_seqs = proxy.i_seqs
    loop.add_row((site_labels[i_seqs[0]],
                  site_labels[i_seqs[1]],
                  site_labels[i_seqs[2]],
                  space_group_info.cif_symmetry_code(sym_ops[0]),
                  space_group_info.cif_symmetry_code(sym_ops[1]),
                  space_group_info.cif_symmetry_code(sym_ops[2]),
                  fmt % restraint.angle_ideal,
                  fmt % math.sqrt(1/restraint.weight),
                  fmt % restraint.delta))
  return loop

def dihedrals_as_cif_loop(xray_structure, proxies):
  space_group_info = sgtbx.space_group_info(group=xray_structure.space_group())
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_torsion_atom_site_label_1",
    "_restr_torsion_atom_site_label_2",
    "_restr_torsion_atom_site_label_3",
    "_restr_torsion_atom_site_label_4",
    "_restr_torsion_site_symmetry_1",
    "_restr_torsion_site_symmetry_2",
    "_restr_torsion_site_symmetry_3",
    "_restr_torsion_site_symmetry_4",
    "_restr_torsion_angle_target",
    "_restr_torsion_weight_param",
    "_restr_torsion_diff",
  ))
  unit_mxs = [sgtbx.rt_mx()]*4
  for proxy in proxies:
    restraint = geometry_restraints.dihedral(
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      proxy=proxy)
    sym_ops = proxy.sym_ops
    if sym_ops is None: sym_ops = unit_mxs
    i_seqs = proxy.i_seqs
    loop.add_row((site_labels[i_seqs[0]],
                  site_labels[i_seqs[1]],
                  site_labels[i_seqs[2]],
                  site_labels[i_seqs[3]],
                  space_group_info.cif_symmetry_code(sym_ops[0]),
                  space_group_info.cif_symmetry_code(sym_ops[1]),
                  space_group_info.cif_symmetry_code(sym_ops[2]),
                  space_group_info.cif_symmetry_code(sym_ops[3]),
                  fmt % restraint.angle_ideal,
                  fmt % math.sqrt(1/restraint.weight),
                  fmt % restraint.delta))
  return loop

def bond_similarity_as_cif_loops(xray_structure, proxies):
  space_group_info = sgtbx.space_group_info(group=xray_structure.space_group())
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_equal_distance_atom_site_label_1",
    "_restr_equal_distance_atom_site_label_2",
    "_restr_equal_distance_site_symmetry_2",
    "_restr_equal_distance_class_id",
  ))
  class_loop = model.loop(header=(
    "_restr_equal_distance_class_class_id",
    "_restr_equal_distance_class_target_weight_param",
    "_restr_equal_distance_class_average",
    "_restr_equal_distance_class_esd",
    "_restr_equal_distance_class_diff_max",
  ))
  class_id = 0
  for proxy in proxies:
    restraint = geometry_restraints.bond_similarity(
      unit_cell=unit_cell,
      sites_cart=sites_cart,
      proxy=proxy)
    class_id += 1
    esd = math.sqrt(flex.sum(flex.pow2(restraint.deltas())) *
                    (1./proxy.i_seqs.size()))
    class_loop.add_row((class_id,
                        fmt % math.sqrt(1/proxy.weights[0]),# assume equal weights
                        fmt % restraint.mean_distance(),
                        fmt % esd,
                        fmt % flex.max_absolute(restraint.deltas())))
    for i in range(proxy.i_seqs.size()):
      i_seq, j_seq = proxy.i_seqs[i]
      if proxy.sym_ops is None:
        sym_op = sgtbx.rt_mx()
      else:
        sym_op = proxy.sym_ops[i]
      loop.add_row((site_labels[i_seq],
                    site_labels[j_seq],
                    space_group_info.cif_symmetry_code(sym_op),
                    class_id))
  return class_loop, loop

def rigid_bond_as_cif_loop(xray_structure, proxies):
  unit_cell = xray_structure.unit_cell()
  sites_cart = xray_structure.sites_cart()
  u_cart = xray_structure.scatterers().extract_u_cart(unit_cell)
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_U_rigid_atom_site_label_1",
    "_restr_U_rigid_atom_site_label_2",
    "_restr_U_rigid_target_weight_param",
    "_restr_U_rigid_U_parallel",
    "_restr_U_rigid_diff",
  ))
  for proxy in proxies:
    restraint = adp_restraints.rigid_bond(
      adp_restraint_params(sites_cart=sites_cart, u_cart=u_cart),
      proxy=proxy)
    loop.add_row((site_labels[proxy.i_seqs[0]],
                  site_labels[proxy.i_seqs[1]],
                  fmt % math.sqrt(1/proxy.weight),
                  fmt % (0.5*(restraint.z_12()+restraint.z_21())),
                  fmt % restraint.delta_z()))
  return loop

def adp_similarity_as_cif_loop(xray_structure, proxies):
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_U_similar_atom_site_label_1",
    "_restr_U_similar_atom_site_label_2",
    "_restr_U_similar_weight_param",
  ))
  for proxy in proxies:
    loop.add_row((site_labels[proxy.i_seqs[0]],
                  site_labels[proxy.i_seqs[1]],
                  fmt % math.sqrt(1/proxy.weight)))
  return loop

def isotropic_adp_as_cif_loop(xray_structure, proxies):
  site_labels = xray_structure.scatterers().extract_labels()
  fmt = "%.4f"
  loop = model.loop(header=(
    "_restr_U_iso_atom_site_label",
    "_restr_U_iso_weight_param",
  ))
  for proxy in proxies:
    loop.add_row((site_labels[proxy.i_seqs[0]], fmt % math.sqrt(1/proxy.weight)))
  return loop

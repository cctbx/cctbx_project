from __future__ import division
import cctbx.geometry_restraints
from cctbx.array_family import flex
import scitbx.graph.tardy_tree
from scitbx import matrix
from libtbx.math_utils import nested_loop

def expand_model_or_conformer_indices(
      indices,
      x_n_seq,
      related_x_i_seqs):
  result = flex.size_t(x_n_seq, x_n_seq)
  for i_seq,i in enumerate(indices):
    result.set_selected(related_x_i_seqs[i_seq], i)
  assert result.count(x_n_seq) == 0
  return result

def expand_site_symmetry_table(
      site_symmetry_table,
      x_n_seq,
      related_x_i_seqs):
  x_indices = flex.size_t(x_n_seq, x_n_seq)
  for i_seq,i in enumerate(site_symmetry_table.indices()):
    x_indices.set_selected(related_x_i_seqs[i_seq], i)
  assert x_indices.count(x_n_seq) == 0
  return cctbx.sgtbx.site_symmetry_table(
    indices=x_indices,
    table=site_symmetry_table.table(),
    special_position_indices=(x_indices == 0).iselection())
  return result

def expand_bond_params_table(
      bond_params_table,
      x_n_seq,
      related_x_i_seqs):
  result = cctbx.geometry_restraints.bond_params_table(x_n_seq)
  for i_seq,bond_params_dict in enumerate(bond_params_table):
    x_i_seqs = related_x_i_seqs[i_seq]
    for j_seq,bond_params in bond_params_dict.items():
      x_j_seqs = related_x_i_seqs[j_seq]
      bond_params_scaled = bond_params.scale_weight(factor=1/
        (len(x_i_seqs)*len(x_j_seqs)))
      for x_i_seq in x_i_seqs:
        for x_j_seq in x_j_seqs:
          if (x_i_seq <= x_j_seq):
            assert x_j_seq not in result[x_i_seq]
            result[x_i_seq][x_j_seq] = bond_params_scaled
          else:
            assert x_i_seq not in result[x_j_seq]
            result[x_j_seq][x_i_seq] = bond_params_scaled
  return result

def expand_pair_sym_table(
      pair_sym_table,
      x_n_seq,
      related_x_i_seqs):
  result = cctbx.crystal.pair_sym_table(x_n_seq)
  for i_seq,pair_sym_dict in enumerate(pair_sym_table):
    x_i_seqs = related_x_i_seqs[i_seq]
    for j_seq,sym_ops in pair_sym_dict.items():
      x_j_seqs = related_x_i_seqs[j_seq]
      for x_i_seq in x_i_seqs:
        for x_j_seq in x_j_seqs:
          if (x_i_seq <= x_j_seq):
            x_sym_ops = result[x_i_seq].setdefault(
              x_j_seq, cctbx.sgtbx.stl_vector_rt_mx())
            for sym_op in sym_ops:
              x_sym_ops.append(sym_op)
          else:
            x_sym_ops = result[x_j_seq].setdefault(
              x_i_seq, cctbx.sgtbx.stl_vector_rt_mx())
            for sym_op in sym_ops:
              x_sym_ops.append(sym_op.inverse())
  return result

def expand_nonbonded_types(
      nonbonded_types,
      x_n_seq,
      related_x_i_seqs):
  result = flex.std_string(x_n_seq)
  for i_seq,s in enumerate(nonbonded_types):
    result.set_selected(related_x_i_seqs[i_seq], s)
  return result

def expand_angle_proxies(
      angle_proxies,
      x_n_seq,
      related_x_i_seqs):
  result = cctbx.geometry_restraints.shared_angle_proxy()
  angle_proxy_t = cctbx.geometry_restraints.angle_proxy
  for proxy in angle_proxies:
    i,j,k = proxy.i_seqs
    proxy_scaled = proxy.scale_weight(factor=1/(
        len(related_x_i_seqs[i])
      * len(related_x_i_seqs[j])
      * len(related_x_i_seqs[k])))
    for x_i in related_x_i_seqs[i]:
      for x_j in related_x_i_seqs[j]:
        for x_k in related_x_i_seqs[k]:
          result.append(angle_proxy_t(
            i_seqs=(x_i,x_j,x_k), proxy=proxy_scaled).sort_i_seqs())
  return result

def expand_dihedral_or_chirality_proxies(
      proxies,
      proxy_type,
      proxy_array_type,
      x_n_seq,
      related_x_i_seqs):
  result = proxy_array_type()
  for proxy in proxies:
    i,j,k,l = proxy.i_seqs
    proxy_scaled = proxy.scale_weight(factor=1/(
        len(related_x_i_seqs[i])
      * len(related_x_i_seqs[j])
      * len(related_x_i_seqs[k])
      * len(related_x_i_seqs[l])))
    for x_i in related_x_i_seqs[i]:
      for x_j in related_x_i_seqs[j]:
        for x_k in related_x_i_seqs[k]:
          for x_l in related_x_i_seqs[l]:
            result.append(proxy_type(
              i_seqs=(x_i,x_j,x_k,x_l), proxy=proxy_scaled).sort_i_seqs())
  return result

def expand_planarity_proxies(
      planarity_proxies,
      x_n_seq,
      related_x_i_seqs):
  result = cctbx.geometry_restraints.shared_planarity_proxy()
  planarity_proxy_t = cctbx.geometry_restraints.planarity_proxy
  for proxy in planarity_proxies:
    x_i_seqs_list = [tuple(related_x_i_seqs[i_seq]) for i_seq in proxy.i_seqs]
    loop_n = [len(x_i_seqs) for x_i_seqs in x_i_seqs_list]
    sym_ops = proxy.sym_ops
    weights = proxy.weights / matrix.col(loop_n).product()
    for loop_i in nested_loop(loop_n):
      i_seqs = flex.size_t([x_i_seqs[i]
        for x_i_seqs,i in zip(x_i_seqs_list, loop_i)])
      result.append(planarity_proxy_t(
        i_seqs=i_seqs, sym_ops=sym_ops, weights=weights).sort_i_seqs())
  return result

class expand(object):

  def __init__(O, labels, sites_cart, masses, geo_manager):
    from scitbx.array_family import shared
    edge_list = geo_manager.simple_edge_list(omit_slack_greater_than=0)
    tt = scitbx.graph.tardy_tree.construct(
      n_vertices=len(sites_cart),
      edge_list=edge_list,
      external_clusters=geo_manager.rigid_clusters_due_to_dihedrals_and_planes(
        constrain_dihedrals_with_sigma_less_than=10))
    orcs = tt.cluster_manager.overlapping_rigid_clusters(
      edge_sets=tt.edge_sets)
    O.labels = flex.std_string()
    O.sites_cart = flex.vec3_double()
    O.masses = flex.double()
    O.indices = []
    O.external_clusters = []
    O.related_x_i_seqs = shared.stl_vector_unsigned(sites_cart.size())
    orca_dof = 0
    for i_orc,orc in enumerate(orcs):
      external_cluster = []
      for i_seq in orc:
        x_i_seq = len(O.sites_cart)
        O.related_x_i_seqs[i_seq].append(x_i_seq)
        external_cluster.append(x_i_seq)
        O.labels.append("%s:%d" % (labels[i_seq], i_orc))
        O.sites_cart.append(sites_cart[i_seq])
        O.masses.append(masses[i_seq])
        O.indices.append((i_orc, i_seq))
      O.external_clusters.append(external_cluster)
      if (len(orc) == 1):
        orca_dof += 3
      elif (len(orc) == 2):
        orca_dof += 5
      else:
        orca_dof += 6
    print "sites_cart.size():", sites_cart.size()
    print "O.sites_cart.size():", O.sites_cart.size()
    print "original Cartesian dof:", sites_cart.size()*3
    print "orca dof:", orca_dof
    for i_seq,x_i_seqs in enumerate(O.related_x_i_seqs):
      O.masses /= len(x_i_seqs)
    # XXX find more direct way to get x_edge_list
    x_edge_list = []
    for i,j in edge_list:
      for x_i in O.related_x_i_seqs[i]:
        i_orc = O.indices[x_i][0]
        for x_j in O.related_x_i_seqs[j]:
          j_orc = O.indices[x_j][0]
          if (i_orc == j_orc):
            x_edge_list.append(tuple(sorted((x_i,x_j))))
    O.tardy_tree = scitbx.graph.tardy_tree.construct(
      sites=O.sites_cart,
      edge_list=x_edge_list,
      external_clusters=O.external_clusters)
    x_model_indices = expand_model_or_conformer_indices(
      indices=geo_manager.model_indices,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_conformer_indices = expand_model_or_conformer_indices(
      indices=geo_manager.conformer_indices,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_site_symmetry_table = expand_site_symmetry_table(
      site_symmetry_table=geo_manager.site_symmetry_table,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_bond_params_table = expand_bond_params_table(
      bond_params_table=geo_manager.bond_params_table,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_shell_sym_tables = []
    for shell_sym_table in geo_manager.shell_sym_tables:
      x_shell_sym_tables.append(expand_pair_sym_table(
        pair_sym_table=shell_sym_table,
        x_n_seq=O.sites_cart.size(),
        related_x_i_seqs=O.related_x_i_seqs))
    x_nonbonded_types = expand_nonbonded_types(
      nonbonded_types=geo_manager.nonbonded_types,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_angle_proxies = expand_angle_proxies(
      angle_proxies=geo_manager.angle_proxies,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_dihedral_proxies = expand_dihedral_or_chirality_proxies(
      proxies=geo_manager.dihedral_proxies,
      proxy_type=cctbx.geometry_restraints.dihedral_proxy,
      proxy_array_type=cctbx.geometry_restraints.shared_dihedral_proxy,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_chirality_proxies = expand_dihedral_or_chirality_proxies(
      proxies=geo_manager.chirality_proxies,
      proxy_type=cctbx.geometry_restraints.chirality_proxy,
      proxy_array_type=cctbx.geometry_restraints.shared_chirality_proxy,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    x_planarity_proxies = expand_planarity_proxies(
      planarity_proxies=geo_manager.planarity_proxies,
      x_n_seq=O.sites_cart.size(),
      related_x_i_seqs=O.related_x_i_seqs)
    O.geo_manager = cctbx.geometry_restraints.manager.manager(
      crystal_symmetry=geo_manager.crystal_symmetry,
      model_indices=x_model_indices,
      conformer_indices=x_conformer_indices,
      site_symmetry_table=x_site_symmetry_table,
      bond_params_table=x_bond_params_table,
      shell_sym_tables=x_shell_sym_tables,
      nonbonded_params=geo_manager.nonbonded_params,
      nonbonded_types=x_nonbonded_types,
      nonbonded_function=geo_manager.nonbonded_function,
      nonbonded_distance_cutoff=geo_manager.nonbonded_distance_cutoff,
      nonbonded_buffer=geo_manager.nonbonded_buffer,
      angle_proxies=x_angle_proxies,
      dihedral_proxies=x_dihedral_proxies,
      chirality_proxies=x_chirality_proxies,
      planarity_proxies=x_planarity_proxies,
      plain_pairs_radius=geo_manager.plain_pairs_radius,
      max_reasonable_bond_distance=geo_manager.max_reasonable_bond_distance,
      min_cubicle_edge=geo_manager.min_cubicle_edge)
    O.tardy_tree_rmsd_calculator = None

  def rmsd_calculator(O, tardy_tree_rmsd_calculator):
    O.tardy_tree_rmsd_calculator = tardy_tree_rmsd_calculator
    return O.rmsd_calculation

  def rmsd_calculation(O, sites_cart_1, sites_cart_2):
    return O.tardy_tree_rmsd_calculator(
      average_sites_cart(
        related_x_i_seqs=O.related_x_i_seqs, x_sites_cart=sites_cart_1),
      average_sites_cart(
        related_x_i_seqs=O.related_x_i_seqs, x_sites_cart=sites_cart_2))

def average_sites_cart(related_x_i_seqs, x_sites_cart):
  result = flex.vec3_double()
  for x_i_seqs in related_x_i_seqs:
    result.append(x_sites_cart.select(x_i_seqs).mean())
  return result

from cctbx.crystal import minimization
from cctbx import restraints
from cctbx import xray
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx import matrix as mx
from scitbx.python_utils.misc import adopt_init_args
import scitbx.lbfgs
from libtbx.itertbx import count

if (1):
  flex.set_random_seed(0)

class restraint_parameters:

  def __init__(self, distance_ideal, weight):
    adopt_init_args(self, locals())

restraint_parameters_si_o = restraint_parameters(1.61, 2.0)
restraint_parameters_o_si_o = restraint_parameters(2.629099, 0.41)
restraint_parameters_si_o_si = restraint_parameters(3.070969, 0.2308)

class add_oxygen:

  def __init__(self, si_structure, si_pair_asu_table):
    self.structure = si_structure.deep_copy_scatterers()
    self.bond_sym_table = crystal.pair_sym_table(
      si_pair_asu_table.table().size())
    sites_frac = si_structure.sites_frac()
    si_asu_mappings = si_pair_asu_table.asu_mappings()
    i_oxygen = count(1)
    for i_seq,asu_dict in enumerate(si_pair_asu_table.table()):
      rt_mx_i_inv = si_asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = sites_frac[i_seq]
      for j_seq,j_sym_groups in asu_dict.items():
        for i_group,j_sym_group in enumerate(j_sym_groups):
          if (j_seq < i_seq): continue
          j_sym = j_sym_group[0]
          rt_mx_ji = rt_mx_i_inv.multiply(
            si_asu_mappings.get_rt_mx(j_seq, j_sym))
          site_frac_ji = rt_mx_ji * sites_frac[j_seq]
          bond_center = (mx.col(site_frac_i) + mx.col(site_frac_ji)) / 2
          i_seq_o = self.structure.scatterers().size()
          self.structure.add_scatterer(xray.scatterer(
            label="O%d"%i_oxygen.next(),
            site=bond_center))
          self.bond_sym_table[i_seq].setdefault(i_seq_o).append(
            sgtbx.rt_mx(1,1))
          self.bond_sym_table[j_seq].setdefault(i_seq_o).append(
            rt_mx_ji.inverse_cancel())

def make_o_si_o_asu_table(si_o_structure, si_o_bond_asu_table):
  scatterers = si_o_structure.scatterers()
  asu_mappings = si_o_bond_asu_table.asu_mappings()
  o_si_o_asu_table = crystal.pair_asu_table(
    asu_mappings=asu_mappings)
  for i_seq,asu_dict in enumerate(si_o_bond_asu_table.table()):
    if (scatterers[i_seq].scattering_type != "Si"): continue
    pair_list = []
    for j_seq,j_sym_groups in asu_dict.items():
      if (scatterers[j_seq].scattering_type != "O"): continue
      for i_group,j_sym_group in enumerate(j_sym_groups):
        for j_sym in j_sym_group:
          pair_list.append((j_seq,j_sym))
    for i_jj1 in xrange(0,len(pair_list)-1):
      jj1 = pair_list[i_jj1]
      rt_mx_jj1_inv = asu_mappings.get_rt_mx(*jj1).inverse()
      for i_jj2 in xrange(i_jj1+1,len(pair_list)):
        jj2 = pair_list[i_jj2]
        rt_mx_jj21 = rt_mx_jj1_inv.multiply(asu_mappings.get_rt_mx(*jj2))
        o_si_o_asu_table.add_pair(
          i_seq=jj1[0],
          j_seq=jj2[0],
          rt_mx_ji=rt_mx_jj21)
  return o_si_o_asu_table

def get_all_proxies(
      structure,
      bond_asu_table,
      bonded_distance_cutoff,
      nonbonded_distance_cutoff,
      minimal=00000):
  bond_asu_proxies = restraints.shared_bond_asu_proxy()
  repulsion_asu_proxies = restraints.shared_repulsion_asu_proxy()
  pair_generator = crystal.neighbors_fast_pair_generator(
    asu_mappings=bond_asu_table.asu_mappings(),
    distance_cutoff=max(bonded_distance_cutoff, nonbonded_distance_cutoff),
    minimal=minimal)
  for pair in pair_generator:
    if (pair in bond_asu_table):
      bond_asu_proxies.append(restraints.bond_asu_proxy(
        pair=pair, distance_ideal=0, weight=0))
    elif (pair.dist_sq**.5 <= nonbonded_distance_cutoff):
      repulsion_asu_proxies.append(restraints.repulsion_asu_proxy(
        pair=pair, vdw_radius=-1))
  return bond_asu_proxies, repulsion_asu_proxies

def edit_bond_asu_proxies(structure, asu_mappings, bond_asu_proxies):
  scatterers = structure.scatterers()
  for proxy in bond_asu_proxies:
    edit_bond_proxy(
      scatterers=scatterers,
      i_seqs=(proxy.i_seq, proxy.j_seq),
      proxy=proxy)

def edit_bond_proxy(scatterers, i_seqs, proxy):
  scattering_types = [scatterers[i].scattering_type for i in i_seqs]
  scattering_types.sort()
  if (scattering_types == ["Si", "Si"]):
    proxy.distance_ideal = restraint_parameters_si_o_si.distance_ideal
    proxy.weight = restraint_parameters_si_o_si.weight
  elif (scattering_types == ["O", "Si"]):
    proxy.distance_ideal = restraint_parameters_si_o.distance_ideal
    proxy.weight = restraint_parameters_si_o.weight
  elif (scattering_types == ["O", "O"]):
    proxy.distance_ideal = restraint_parameters_o_si_o.distance_ideal
    proxy.weight = restraint_parameters_o_si_o.weight
  else:
    raise AssertionError("Unknown scattering type pair.")

def edit_repulsion_asu_proxies(structure, asu_mappings, repulsion_asu_proxies):
  scatterers = structure.scatterers()
  for proxy in repulsion_asu_proxies:
    i_seqs = proxy.i_seq, proxy.j_seq
    scattering_types = [scatterers[i].scattering_type for i in i_seqs]
    scattering_types.sort()
    if (scattering_types == ["Si", "Si"]):
      proxy.vdw_radius = 3.1
    elif (scattering_types == ["O", "Si"]):
      proxy.vdw_radius = 1.5
    elif (scattering_types == ["O", "O"]):
      proxy.vdw_radius = 2.0
    else:
      raise AssertionError("Unknown scattering type pair.")

class show_pairs:

  def __init__(self, structure, pair_asu_table):
    self.distances = flex.double()
    self.pair_counts = flex.size_t()
    unit_cell = structure.unit_cell()
    scatterers = structure.scatterers()
    sites_frac = structure.sites_frac()
    asu_mappings = pair_asu_table.asu_mappings()
    for i_seq,asu_dict in enumerate(pair_asu_table.table()):
      rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
      site_frac_i = sites_frac[i_seq]
      pair_count = 0
      dists = flex.double()
      j_seq_i_group = []
      for j_seq,j_sym_groups in asu_dict.items():
        site_frac_j = sites_frac[j_seq]
        for i_group,j_sym_group in enumerate(j_sym_groups):
          pair_count += j_sym_group.size()
          j_sym = j_sym_group[0]
          rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
          distance = unit_cell.distance(site_frac_i, rt_mx_ji * site_frac_j)
          dists.append(distance)
          j_seq_i_group.append((j_seq,i_group))
      s = "%s(%d):" % (scatterers[i_seq].label, i_seq+1)
      s = "%-15s pair count: %3d" % (s, pair_count)
      print "%-32s"%s, "<<"+",".join([" %7.4f" % x for x in site_frac_i])+">>"
      permutation = flex.sort_permutation(dists)
      for j_seq,i_group in flex.select(j_seq_i_group, permutation):
        site_frac_j = sites_frac[j_seq]
        j_sym_groups = asu_dict[j_seq]
        j_sym_group = j_sym_groups[i_group]
        for i_j_sym,j_sym in enumerate(j_sym_group):
          rt_mx_ji = rt_mx_i_inv.multiply(
            asu_mappings.get_rt_mx(j_seq, j_sym))
          site_frac_ji = rt_mx_ji * site_frac_j
          distance = unit_cell.distance(site_frac_i, site_frac_ji)
          self.distances.append(distance)
          print "  %-10s" % ("%s(%d):" % (scatterers[j_seq].label, j_seq+1)),
          print "%8.4f" % distance,
          if (i_j_sym != 0):
            s = "sym. equiv."
          else:
            s = "           "
          s += " (" + ",".join([" %7.4f" % x for x in site_frac_ji]) +")"
          print s
      if (pair_count == 0):
        print "  no neighbors"
      self.pair_counts.append(pair_count)

def show_nonbonded_interactions(structure, asu_mappings, nonbonded_proxies):
  distances = flex.double()
  unit_cell = structure.unit_cell()
  scatterers = structure.scatterers()
  sites_frac = structure.sites_frac()
  for proxy in nonbonded_proxies:
    i_seq, j_seq, j_sym = proxy.i_seq, proxy.j_seq, proxy.j_sym
    rt_mx_i_inv = asu_mappings.get_rt_mx(i_seq, 0).inverse()
    rt_mx_ji = rt_mx_i_inv.multiply(asu_mappings.get_rt_mx(j_seq, j_sym))
    pair_labels = "%s(%d) - %s(%d):" % (
      scatterers[i_seq].label, i_seq+1,
      scatterers[j_seq].label, j_seq+1)
    distance = unit_cell.distance(
      sites_frac[i_seq], rt_mx_ji*sites_frac[j_seq])
    distances.append(distance)
    print "%-20s %8.4f" % (pair_labels, distance)
  return distances

def get_distances(unit_cell, sites_cart, pair_sym_table):
  sites_frac = unit_cell.fractionalization_matrix() * sites_cart
  distances = flex.double()
  for i_seq,pair_sym_dict in enumerate(pair_sym_table):
    site_i = sites_frac[i_seq]
    for j_seq,rt_mx_list in pair_sym_dict.items():
      site_j = sites_frac[j_seq]
      for rt_mx in rt_mx_list:
        distances.append(unit_cell.distance(site_i, rt_mx*site_j))
  return distances

def distance_and_repulsion_least_squares(
      si_structure,
      distance_cutoff,
      nonbonded_distance_cutoff,
      n_trials=1,
      connectivities=None):
  assert n_trials > 0
  si_structure.show_summary().show_scatterers()
  print
  si_asu_mappings = si_structure.asu_mappings(
    buffer_thickness=distance_cutoff)
  si_pair_asu_table = crystal.pair_asu_table(
    asu_mappings=si_asu_mappings)
  si_pair_asu_table.add_all_pairs(distance_cutoff=distance_cutoff)
  si_pairs = show_pairs(
    structure=si_structure,
    pair_asu_table=si_pair_asu_table)
  if (connectivities is not None):
    assert list(si_pairs.pair_counts) == connectivities
  print
  si_o = add_oxygen(
    si_structure=si_structure,
    si_pair_asu_table=si_pair_asu_table)
  si_o.structure.show_summary().show_scatterers()
  print
  assert nonbonded_distance_cutoff \
       > flex.max(si_pairs.distances)/2.*(1-1.e-6)
  si_o_asu_mappings = si_o.structure.asu_mappings(
    buffer_thickness=nonbonded_distance_cutoff)
  si_o_bond_asu_table = crystal.pair_asu_table(
    asu_mappings=si_o_asu_mappings)
  si_o_bond_asu_table.add_pair_sym_table(sym_table=si_o.bond_sym_table)
  si_o_bonds = show_pairs(
    structure=si_o.structure,
    pair_asu_table=si_o_bond_asu_table)
  n_si = si_pairs.pair_counts.size()
  n_si_o = si_o_bonds.pair_counts.size()
  assert si_o_bonds.pair_counts[:n_si].all_eq(si_pairs.pair_counts)
  assert si_o_bonds.pair_counts[n_si:].count(2) == n_si_o-n_si
  print
  o_si_o_asu_table = make_o_si_o_asu_table(
    si_o_structure=si_o.structure,
    si_o_bond_asu_table=si_o_bond_asu_table)
  o_si_o_pairs = show_pairs(
    structure=si_o.structure,
    pair_asu_table=o_si_o_asu_table)
  assert o_si_o_pairs.pair_counts[:n_si].all_eq(0)
  if (si_pairs.pair_counts.count(4) == n_si):
    assert o_si_o_pairs.pair_counts[n_si:].all_eq(6)
  print
  if (1):
    si_o_bond_asu_table.add_pair_sym_table(
      sym_table=si_pair_asu_table.extract_pair_sym_table())
  if (1):
    si_o_bond_asu_table.add_pair_sym_table(
      sym_table=o_si_o_asu_table.extract_pair_sym_table())
  bond_sym_table = si_o_bond_asu_table.extract_pair_sym_table()
  minimized = None
  for i_trial in xrange(n_trials):
    trial_structure = si_o.structure.deep_copy_scatterers()
    if (i_trial > 0):
      n_scatterers = trial_structure.scatterers().size()
      trial_structure.set_sites_frac(flex.vec3_double(flex.random_double(
        size=n_scatterers*3)))
      trial_structure.apply_symmetry_sites()
    trial_minimized = minimization.lbfgs(
      structure=trial_structure,
      bond_sym_table=bond_sym_table,
      nonbonded_distance_cutoff=nonbonded_distance_cutoff,
      lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=100))
    print "i_trial, target value: %d, %.6g" % (
      i_trial, trial_minimized.final_target_result.target())
    if (minimized is None or       minimized.final_target_result.target()
                           > trial_minimized.final_target_result.target()):
      minimized = trial_minimized
      minimized_structure = trial_structure
      best_i_trial = i_trial
  assert minimized is not None
  print
  print "Energies at start:"
  minimized.first_target_result.show()
  print
  print "Energies at end:"
  minimized.final_target_result.show()
  print
  print "Final target value (i_trial=%d): %.6g" % (
    best_i_trial, minimized.final_target_result.target())
  if (minimized.final_target_result.target() > 0.1):
    print "WARNING: LARGE final target value: %.6g" % (
      minimized.final_target_result.target())
  print
  show_pairs(
    structure=minimized_structure,
    pair_asu_table=si_o_bond_asu_table)
  print

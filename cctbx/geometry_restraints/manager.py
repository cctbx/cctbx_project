from __future__ import division
from cctbx import geometry_restraints
import cctbx.geometry_restraints.flags
import cctbx.geometry_restraints.energies
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import store
from libtbx import introspection
from libtbx import adopt_init_args
from libtbx import dict_with_default_0
from libtbx.utils import Sorry
import sys, math, StringIO

import boost.python
boost.python.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import reindexing_array

#from mmtbx.geometry_restraints.hbond import get_simple_bonds

class manager(object):

  def __init__(self,
        crystal_symmetry=None,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        site_symmetry_table=None,
        donor_acceptor_excl_groups=None,
        bond_params_table=None,
        shell_sym_tables=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_charges=None,
        nonbonded_function=None,
        nonbonded_distance_cutoff=None,
        nonbonded_buffer=1,
        angle_proxies=None,
        dihedral_proxies=None,
        reference_dihedral_proxies=None,
        ncs_dihedral_proxies=None,
        chirality_proxies=None,
        planarity_proxies=None,
        parallelity_proxies=None,
        generic_restraints_manager=None,
        external_energy_function=None,
        plain_pairs_radius=None,
        max_reasonable_bond_distance=None,
        min_cubicle_edge=5,
        hbonds_in_bond_list=None):
    if (site_symmetry_table is not None): assert crystal_symmetry is not None
    if (bond_params_table is not None and site_symmetry_table is not None):
      assert bond_params_table.size() == site_symmetry_table.indices().size()
    if (shell_sym_tables is not None and site_symmetry_table is not None):
      assert len(shell_sym_tables) > 0
      assert shell_sym_tables[0].size() == site_symmetry_table.indices().size()
    if (nonbonded_types is not None and site_symmetry_table is not None):
      assert nonbonded_types.size() == site_symmetry_table.indices().size()
    if (nonbonded_types is not None) and (nonbonded_charges is not None) :
      assert (nonbonded_charges.size() == nonbonded_types.size())
    adopt_init_args(self, locals())
    self.reset_internals()

  def reset_internals(self):
    self._sites_cart_used_for_pair_proxies = None
    self._flags_bond_used_for_pair_proxies = False
    self._flags_nonbonded_used_for_pair_proxies = False
    self._pair_proxies = None
    self.plain_pair_sym_table = None
    self.nonbonded_distance_cutoff_was_determined_automatically = False
    self.adjusted_nonbonded_distance_cutoff = self.nonbonded_distance_cutoff
    self.effective_nonbonded_buffer = self.nonbonded_buffer
    self.n_updates_pair_proxies = 0

  def replace_site_symmetry(self, new_site_symmetry_table):
    assert self.site_symmetry_table is not None
    self.site_symmetry_table = new_site_symmetry_table
    self.reset_internals()

  def simple_edge_list(self, omit_slack_greater_than=0):
    assert self.bond_params_table is not None
    assert self.shell_sym_tables is not None
    result = []
    bpt = self.bond_params_table
    for i,j in self.shell_sym_tables[0].simple_edge_list():
      assert i < j
      if (bpt[i][j].slack > omit_slack_greater_than): continue
      result.append((i,j))
    return result

  def rigid_clusters_due_to_dihedrals_and_planes(self,
        constrain_dihedrals_with_sigma_less_than):
    result = []
    if (self.dihedral_proxies is not None):
      assert constrain_dihedrals_with_sigma_less_than > 0
      weight_limit = 1.0 / constrain_dihedrals_with_sigma_less_than**2
      for proxy in self.dihedral_proxies:
        if (proxy.weight > weight_limit):
          result.append(proxy.i_seqs)
    if (self.planarity_proxies is not None):
      for proxy in self.planarity_proxies:
        result.append(proxy.i_seqs)
    return result

  def construct_tardy_tree(self,
        sites=None,
        sites_cart=None,
        selection=None,
        omit_bonds_with_slack_greater_than=0,
        constrain_dihedrals_with_sigma_less_than=10,
        near_singular_hinges_angular_tolerance_deg=5):
    assert [sites_cart, sites].count(None) == 1
    from scitbx.graph import tardy_tree
    from scitbx import matrix
    if (sites is None):
      sites = matrix.col_list(sites_cart)
    if (selection is not None):
      if (not isinstance(selection, flex.bool)):
        selection = flex.bool(len(sites), selection)
      fixed_vertices_bool = ~selection
    else:
      fixed_vertices_bool = flex.bool(len(sites), False)
    if (self.site_symmetry_table is not None):
      fixed_vertices_bool.set_selected(
        self.site_symmetry_table.special_position_indices(), True)
    return tardy_tree.construct(
      sites=sites,
      edge_list=self.simple_edge_list(
        omit_slack_greater_than
          =omit_bonds_with_slack_greater_than),
      external_clusters=self.rigid_clusters_due_to_dihedrals_and_planes(
        constrain_dihedrals_with_sigma_less_than
          =constrain_dihedrals_with_sigma_less_than),
      fixed_vertices=fixed_vertices_bool.iselection(),
      near_singular_hinges_angular_tolerance_deg
        =near_singular_hinges_angular_tolerance_deg)

  def reduce_for_tardy(self,
        tardy_tree,
        omit_bonds_with_slack_greater_than=0,
        include_den_restraints=False):
    from cctbx import sgtbx
    from scitbx.graph.utils import construct_edge_sets
    #
    assert tardy_tree.n_vertices == self.bond_params_table.size()
    loop_edge_sets = construct_edge_sets(
      n_vertices=tardy_tree.n_vertices,
      edge_list=tardy_tree.cluster_manager.loop_edges,
      assert_i_lt_j=False)
    #
    def get():
      result = crystal.pair_sym_table(tardy_tree.n_vertices)
      bond_params_table = self.bond_params_table
      for i,pair_sym_dict in enumerate(self.shell_sym_tables[0]):
        reduced_pair_sym_dict = result[i]
        for j,sym_ops in pair_sym_dict.items():
          reduced_sym_ops = sgtbx.stl_vector_rt_mx()
          for sym_op in sym_ops:
            if (sym_op.is_unit_mx()):
              if (  bond_params_table[i][j].slack
                  > omit_bonds_with_slack_greater_than):
                continue
              if (j not in loop_edge_sets[i]):
                continue
            reduced_sym_ops.append(sym_op)
          if (reduced_sym_ops.size() != 0):
            reduced_pair_sym_dict[j] = reduced_sym_ops
      return [result]
    reduced_shell_sym_tables = get()
    #
    def get():
      result = geometry_restraints.shared_angle_proxy()
      for proxy in self.angle_proxies:
        i,j,k = proxy.i_seqs
        if (   j in loop_edge_sets[i]
            or j in loop_edge_sets[k]):
          result.append(proxy)
      if (result.size() == 0):
        return None
      return result
    reduced_angle_proxies = get()
    #
    def get():
      result = geometry_restraints.shared_dihedral_proxy()
      ec = tardy_tree.cluster_manager.edge_classifier()
      for proxy in self.dihedral_proxies:
        if (ec(edge=proxy.i_seqs[1:3]) in ["hinge", "loop"]):
          result.append(proxy)
      if (result.size() == 0):
        return None
      return result
    reduced_dihedral_proxies = get()
    #
    if not include_den_restraints:
      return manager(
        crystal_symmetry=self.crystal_symmetry,
        site_symmetry_table=self.site_symmetry_table,
        bond_params_table=self.bond_params_table,
        shell_sym_tables=reduced_shell_sym_tables,
        angle_proxies=reduced_angle_proxies,
        dihedral_proxies=reduced_dihedral_proxies,
        ncs_dihedral_proxies=self.ncs_dihedral_proxies,
        hbonds_in_bond_list=self.hbonds_in_bond_list)
    else:
      return manager(
        crystal_symmetry=self.crystal_symmetry,
        site_symmetry_table=self.site_symmetry_table,
        bond_params_table=self.bond_params_table,
        shell_sym_tables=reduced_shell_sym_tables,
        angle_proxies=reduced_angle_proxies,
        dihedral_proxies=reduced_dihedral_proxies,
        ncs_dihedral_proxies=self.ncs_dihedral_proxies,
        generic_restraints_manager=self.generic_restraints_manager,
        hbonds_in_bond_list=self.hbonds_in_bond_list)

  def sites_cart_used_for_pair_proxies(self):
    return self._sites_cart_used_for_pair_proxies

  def new_including_isolated_sites(self,
        n_additional_sites,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        donor_acceptor_excl_groups=None,
        site_symmetry_table=None,
        nonbonded_types=None,
        nonbonded_charges=None):
    assert n_additional_sites >= 0
    assert (model_indices is None) == (self.model_indices is None)
    assert (conformer_indices is None) == (self.conformer_indices is None)
    assert (donor_acceptor_excl_groups is None) == \
           (self.donor_acceptor_excl_groups is None)
    assert (sym_excl_indices is None) == (self.sym_excl_indices is None)
    assert (site_symmetry_table is None) == (self.site_symmetry_table is None)
    assert (nonbonded_types is None) == (self.nonbonded_types is None)
    if (self.model_indices is not None):
      assert model_indices.size() == n_additional_sites
      model_indices = self.model_indices.concatenate(
        model_indices)
    if (self.conformer_indices is not None):
      assert conformer_indices.size() == n_additional_sites
      conformer_indices = self.conformer_indices.concatenate(
        conformer_indices)
    if (self.donor_acceptor_excl_groups is not None):
      assert donor_acceptor_excl_groups.size() == n_additional_sites
      donor_acceptor_excl_groups = self.donor_acceptor_excl_groups.concatenate(
        donor_acceptor_excl_groups)
    if (self.sym_excl_indices is not None):
      assert sym_excl_indices.size() == n_additional_sites
      sym_excl_indices = self.sym_excl_indices.concatenate(
        sym_excl_indices)
    if (self.site_symmetry_table is not None):
      assert site_symmetry_table.indices().size() == n_additional_sites
      # XXX should become site_symmetry_table.concatenate()
      new_site_symmetry_table = self.site_symmetry_table.deep_copy()
      new_site_symmetry_table.reserve(new_site_symmetry_table.indices().size()
                                    + n_additional_sites)
      for i_seq in xrange(n_additional_sites):
        new_site_symmetry_table.process(site_symmetry_table.get(i_seq))
      site_symmetry_table = new_site_symmetry_table
    bond_params_table = None
    if (self.bond_params_table is not None):
      bond_params_table = self.bond_params_table.deep_copy()
      bond_params_table.extend(geometry_restraints.bond_params_table(
        n_additional_sites))
    shell_sym_tables = None
    if (self.shell_sym_tables is not None):
      shell_sym_tables = []
      for shell_sym_table in self.shell_sym_tables:
        shell_sym_table = shell_sym_table.deep_copy()
        shell_sym_table.extend(crystal.pair_sym_table(n_additional_sites))
        shell_sym_tables.append(shell_sym_table)
    if (self.nonbonded_types is not None):
      assert nonbonded_types.size() == n_additional_sites
      nonbonded_types = self.nonbonded_types.concatenate(
        nonbonded_types)
    if (self.nonbonded_charges is not None) :
      if (nonbonded_charges is None) :
        nonbonded_charges = flex.int(n_additional_sites, 0)
      assert (nonbonded_charges.size() == n_additional_sites)
      nonbonded_charges = self.nonbonded_charges.concatenate(
        nonbonded_charges)
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      model_indices=model_indices,
      conformer_indices=conformer_indices,
      sym_excl_indices=sym_excl_indices,
      donor_acceptor_excl_groups=donor_acceptor_excl_groups,
      site_symmetry_table=site_symmetry_table,
      bond_params_table=bond_params_table,
      shell_sym_tables=shell_sym_tables,
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=nonbonded_types,
      nonbonded_charges=nonbonded_charges,
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      reference_dihedral_proxies=self.reference_dihedral_proxies,
      generic_restraints_manager=self.generic_restraints_manager,
      ncs_dihedral_proxies=self.ncs_dihedral_proxies,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      parallelity_proxies=self.parallelity_proxies,
      plain_pairs_radius=self.plain_pairs_radius,
      hbonds_in_bond_list=self.hbonds_in_bond_list)

  def select(self, selection=None, iselection=None):
    assert [selection, iselection].count(None) == 1
    n_seqs = dict_with_default_0()
    if (iselection is None):
      iselection = selection.iselection()
      n_seqs[selection.size()] += 1
    selected_model_indices = None
    if (self.model_indices is not None):
      selected_model_indices = self.model_indices.select(
        iselection)
      n_seqs[self.model_indices.size()] += 1
    selected_conformer_indices = None
    if (self.conformer_indices is not None):
      selected_conformer_indices = self.conformer_indices.select(
        iselection)
      n_seqs[self.conformer_indices.size()] += 1
    selected_sym_excl_indices = None
    if (self.sym_excl_indices is not None):
      selected_sym_excl_indices = self.sym_excl_indices.select(
        iselection)
      n_seqs[self.sym_excl_indices.size()] += 1
    selected_donor_acceptor_excl_groups = None
    if (self.donor_acceptor_excl_groups is not None):
      selected_donor_acceptor_excl_groups = \
        self.donor_acceptor_excl_groups.select(iselection)
      n_seqs[self.donor_acceptor_excl_groups.size()] += 1
    selected_site_symmetry_table = None
    if (self.site_symmetry_table is not None):
      selected_site_symmetry_table = self.site_symmetry_table.select(
        iselection)
      n_seqs[self.site_symmetry_table.indices().size()] += 1
    selected_bond_params_table = None
    if (self.bond_params_table is not None):
      selected_bond_params_table = self.bond_params_table.proxy_select(
        iselection)
      n_seqs[self.bond_params_table.size()] += 1
    selected_shell_sym_tables = None
    if (self.shell_sym_tables is not None):
      selected_shell_sym_tables = [shell_sym_table.proxy_select(iselection)
        for shell_sym_table in self.shell_sym_tables]
      if (len(self.shell_sym_tables) > 0):
        n_seqs[self.shell_sym_tables[0].size()] += 1
    selected_nonbonded_types = None
    if (self.nonbonded_types is not None):
      selected_nonbonded_types = self.nonbonded_types.select(
        iselection)
      n_seqs[self.nonbonded_types.size()] += 1
    selected_nonbonded_charges = None
    if (self.nonbonded_charges is not None) :
      selected_nonbonded_charges = self.nonbonded_charges.select(
        iselection)
      n_seqs[self.nonbonded_charges.size()] += 1
    n_seq = None
    def get_n_seq():
      if (len(n_seqs) == 0):
        raise RuntimeError("Cannot determine n_seq.")
      if (len(n_seqs) != 1):
        raise RuntimeError("Selection size mismatches: %s." % str(n_seqs))
      return n_seqs.keys()[0]
    selected_angle_proxies = None
    if (self.angle_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_angle_proxies = self.angle_proxies.proxy_select(
        n_seq, iselection)
    selected_dihedral_proxies = None
    if (self.dihedral_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_dihedral_proxies = self.dihedral_proxies.proxy_select(
        n_seq, iselection)
    selected_reference_dihedral_proxies = None
    if (self.reference_dihedral_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_reference_dihedral_proxies = self.reference_dihedral_proxies.proxy_select(
        n_seq, iselection)
    selected_ncs_dihedral_proxies = None
    if (self.ncs_dihedral_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_ncs_dihedral_proxies = self.ncs_dihedral_proxies.proxy_select(
        n_seq, iselection)
    selected_chirality_proxies = None
    if (self.chirality_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_chirality_proxies = self.chirality_proxies.proxy_select(
        n_seq, iselection)
    selected_planarity_proxies = None
    if (self.planarity_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_planarity_proxies = self.planarity_proxies.proxy_select(
        n_seq, iselection)
    selected_parallelity_proxies = None
    if (self.parallelity_proxies is not None):
      if (n_seq is None): n_seq = get_n_seq()
      selected_parallelity_proxies = self.parallelity_proxies.proxy_select(
        n_seq, iselection)
    generic_restraints_manager = None
    if (self.generic_restraints_manager is not None) :
      generic_restraints_manager = self.generic_restraints_manager.select(
        n_seq, iselection)
    selected_hbonds_in_list=[]
    if self.hbonds_in_bond_list is not None:
      r_a = reindexing_array(n_seq, iselection.as_int())
      for pair in self.hbonds_in_bond_list:
        if pair[0] in iselection and pair[1] in iselection:
          reindexed_pair = (r_a[pair[0]], r_a[pair[1]])
          selected_hbonds_in_list.append(reindexed_pair)
      if len(selected_hbonds_in_list) == 0:
        selected_hbonds_in_list = None
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      model_indices=selected_model_indices,
      conformer_indices=selected_conformer_indices,
      sym_excl_indices=selected_sym_excl_indices,
      donor_acceptor_excl_groups=selected_donor_acceptor_excl_groups,
      site_symmetry_table=selected_site_symmetry_table,
      bond_params_table=selected_bond_params_table,
      shell_sym_tables=selected_shell_sym_tables,
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=selected_nonbonded_types,
      nonbonded_charges=selected_nonbonded_charges,
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=selected_angle_proxies,
      dihedral_proxies=selected_dihedral_proxies,
      reference_dihedral_proxies=selected_reference_dihedral_proxies,
      generic_restraints_manager=generic_restraints_manager,
      ncs_dihedral_proxies=selected_ncs_dihedral_proxies,
      chirality_proxies=selected_chirality_proxies,
      planarity_proxies=selected_planarity_proxies,
      parallelity_proxies=selected_parallelity_proxies,
      plain_pairs_radius=self.plain_pairs_radius,
      hbonds_in_bond_list=selected_hbonds_in_list)

  def discard_symmetry(self, new_unit_cell):
    assert self.site_symmetry_table is not None #XXX lazy
    assert self.shell_sym_tables is not None #XXX lazy
    return manager(
      crystal_symmetry=crystal.symmetry(
        unit_cell=new_unit_cell,
        space_group_symbol="P1"),
      model_indices=self.model_indices,
      conformer_indices=self.conformer_indices,
      sym_excl_indices=None,
      donor_acceptor_excl_groups=self.donor_acceptor_excl_groups,
      site_symmetry_table=self.site_symmetry_table.discard_symmetry(),
      bond_params_table=self.bond_params_table,
      shell_sym_tables=[t.discard_symmetry() for t in self.shell_sym_tables],
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=self.nonbonded_types,
      nonbonded_charges=self.nonbonded_charges,
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      reference_dihedral_proxies=self.reference_dihedral_proxies,
      generic_restraints_manager=self.generic_restraints_manager,
      ncs_dihedral_proxies=self.ncs_dihedral_proxies,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      parallelity_proxies=self.parallelity_proxies,
      plain_pairs_radius=self.plain_pairs_radius,
      hbonds_in_bond_list=self.hbonds_in_bond_list)

  def remove_angles_in_place(self, selection):
    self.angle_proxies = self.angle_proxies.proxy_remove(
      selection=selection)

  def remove_dihedrals_in_place(self, selection):
    self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_reference_dihedrals_in_place(self, selection):
    self.reference_dihedral_proxies = self.reference_dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_ncs_dihedrals_in_place(self, selection):
    self.ncs_dihedral_proxies = self.ncs_dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_chiralities_in_place(self, selection):
    self.chirality_proxies = self.chirality_proxies.proxy_remove(
      selection=selection)

  def remove_planarities_in_place(self, selection):
    self.planarity_proxies = self.planarity_proxies.proxy_remove(
      selection=selection)

  def remove_parallelities_in_place(self, selection):
    self.parallelity_proxies = self.parallelity_proxies.proxy_remove(
      selection=selection)

  def set_generic_restraints (self, manager) :
    self.generic_restraints_manager = manager

  def set_external_energy_function (self, energy_function) :
    self.external_energy_function = energy_function

  def pair_proxies(self,
        sites_cart=None,
        flags=None,
        asu_is_inside_epsilon=None,
        bonded_distance_cutoff_epsilon=None,
        site_labels=None):
    if (bonded_distance_cutoff_epsilon is None):
      bonded_distance_cutoff_epsilon = 1.e-6
    if (self.crystal_symmetry is None):
      orthogonalization_matrix = None
    else:
      orthogonalization_matrix = self.crystal_symmetry.unit_cell() \
        .orthogonalization_matrix()
    bonded_distance_cutoff = -1
    def check_bonded_distance_cutoff(sites_frac=None, sites_cart=None):
      if (    self.max_reasonable_bond_distance is not None
          and bonded_distance_cutoff > self.max_reasonable_bond_distance):
        lines = format_distances_for_error_message(
          pair_sym_table=self.shell_sym_tables[0],
          larger_than=3,
          orthogonalization_matrix=orthogonalization_matrix,
          sites_frac=sites_frac,
          sites_cart=sites_cart,
          site_labels=site_labels)
        msg = "Bond distance > max_reasonable_bond_distance: %.6g > %.6g" % (
          bonded_distance_cutoff, self.max_reasonable_bond_distance)
        if (len(lines) != 0):
          msg += ":\n  "
          msg += "\n  ".join(lines)
        raise RuntimeError(msg)
    def flags_are_different():
      if (flags is None):
        if (not self._flags_bond_used_for_pair_proxies):
          return True
        if (not self._flags_nonbonded_used_for_pair_proxies):
          return True
      else:
        if (self._flags_bond_used_for_pair_proxies != flags.bond):
          return True
        if (self._flags_nonbonded_used_for_pair_proxies != flags.nonbonded):
          return True
      return False
    if (self.nonbonded_types is None):
      if (    self._pair_proxies is None
          and self.shell_sym_tables is not None
          and self.crystal_symmetry is not None
          and self.site_symmetry_table is not None):
        unit_cell = self.crystal_symmetry.unit_cell()
        sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
        for shell_sym_table in self.shell_sym_tables:
          bonded_distance_cutoff = max(bonded_distance_cutoff,
            flex.max_default(
              values=crystal.get_distances(
                pair_sym_table=shell_sym_table,
                orthogonalization_matrix=orthogonalization_matrix,
                sites_frac=sites_frac),
              default=0))
        check_bonded_distance_cutoff(sites_frac=sites_frac)
        bonded_distance_cutoff *= (1 + bonded_distance_cutoff_epsilon)
        asu_mappings = crystal.symmetry.asu_mappings(self.crystal_symmetry,
          buffer_thickness=bonded_distance_cutoff,
          asu_is_inside_epsilon=asu_is_inside_epsilon)
        asu_mappings.process_sites_frac(
          original_sites=sites_frac,
          site_symmetry_table=self.site_symmetry_table)
        shell_asu_tables = None
        shell_asu_tables = [
          crystal.pair_asu_table(asu_mappings=asu_mappings)
          .add_pair_sym_table(sym_table=shell_sym_table)
          for shell_sym_table in self.shell_sym_tables]
        self.n_updates_pair_proxies += 1
        self._pair_proxies = geometry_restraints.pair_proxies(
          flags=flags,
          bond_params_table=self.bond_params_table,
          min_cubicle_edge=self.min_cubicle_edge,
          shell_asu_tables=shell_asu_tables)
      elif (self._pair_proxies is None):
        self.n_updates_pair_proxies += 1
        self._pair_proxies = geometry_restraints.pair_proxies(
          flags=flags,
          bond_params_table=self.bond_params_table,
          min_cubicle_edge=self.min_cubicle_edge)
    elif (sites_cart is not None
            and (self._sites_cart_used_for_pair_proxies is None
              or flags_are_different()
              or self._sites_cart_used_for_pair_proxies.max_distance(
                sites_cart) > self.effective_nonbonded_buffer)):
      self.n_updates_pair_proxies += 1
      self._sites_cart_used_for_pair_proxies = sites_cart.deep_copy()
      if (flags is None):
        self._flags_bond_used_for_pair_proxies = True
        self._flags_nonbonded_used_for_pair_proxies = True
      else:
        self._flags_bond_used_for_pair_proxies = flags.bond
        self._flags_nonbonded_used_for_pair_proxies = flags.nonbonded
      bonded_distance_cutoff = -1
      if (self.nonbonded_distance_cutoff is None):
        max_vdw_dist = self.nonbonded_params.find_max_vdw_distance(
          nonbonded_types=self.nonbonded_types)
        assert max_vdw_dist > 0
        r = self.nonbonded_function.residual(
          vdw_distance=max_vdw_dist, delta=max_vdw_dist*(1+1.e-6))
        if (r != 0):
          raise RuntimeError(
            "Cannot automatically determine nonbonded_distance_cutoff:"
            " nonbonded_function.residual() not zero beyond maximum VDW"
            " distance.")
        self.nonbonded_distance_cutoff_was_determined_automatically = True
        self.nonbonded_distance_cutoff = max_vdw_dist
        self.adjusted_nonbonded_distance_cutoff = max_vdw_dist
      asu_mappings = None
      shell_asu_tables = None
      while True:
        current_nonbonded_distance_cutoff_plus_buffer \
          = self.adjusted_nonbonded_distance_cutoff \
          + self.nonbonded_buffer
        if (self.crystal_symmetry is None):
          if (bonded_distance_cutoff < 0):
            for shell_sym_table in self.shell_sym_tables:
              bonded_distance_cutoff = max(bonded_distance_cutoff,
                flex.max_default(
                  values=crystal.get_distances(
                    pair_sym_table=shell_sym_table,
                    sites_cart=sites_cart),
                  default=0))
            check_bonded_distance_cutoff(sites_cart=sites_cart)
            bonded_distance_cutoff *= (1 + bonded_distance_cutoff_epsilon)
            asu_mappings = \
              crystal.direct_space_asu.non_crystallographic_asu_mappings(
                sites_cart=sites_cart,
                min_unit_cell_length=
                  2*current_nonbonded_distance_cutoff_plus_buffer)
        else:
          if (   bonded_distance_cutoff < 0
              or self.plain_pairs_radius is not None):
            unit_cell = self.crystal_symmetry.unit_cell()
            sites_frac = unit_cell.fractionalize(sites_cart=sites_cart)
          if (self.plain_pairs_radius is not None):
            self.update_plain_pair_sym_table(sites_frac=sites_frac)
          if (bonded_distance_cutoff < 0):
            for shell_sym_table in self.shell_sym_tables:
              bonded_distance_cutoff = max(bonded_distance_cutoff,
                flex.max_default(
                  values=crystal.get_distances(
                    pair_sym_table=shell_sym_table,
                    orthogonalization_matrix=orthogonalization_matrix,
                    sites_frac=sites_frac),
                  default=0))
            check_bonded_distance_cutoff(sites_frac=sites_frac)
            bonded_distance_cutoff *= (1 + bonded_distance_cutoff_epsilon)
          if (asu_mappings is None
              or asu_mappings.buffer_thickness()
                 < current_nonbonded_distance_cutoff_plus_buffer):
            asu_mappings = crystal.symmetry.asu_mappings(self.crystal_symmetry,
              buffer_thickness=max(
                bonded_distance_cutoff,
                current_nonbonded_distance_cutoff_plus_buffer),
              asu_is_inside_epsilon=asu_is_inside_epsilon)
            asu_mappings.process_sites_frac(
              original_sites=sites_frac,
              site_symmetry_table=self.site_symmetry_table)
            shell_asu_tables = None
        if (shell_asu_tables is None):
          shell_asu_tables = [
            crystal.pair_asu_table(asu_mappings=asu_mappings)
              .add_pair_sym_table(sym_table=shell_sym_table)
                for shell_sym_table in self.shell_sym_tables]
        self._pair_proxies = geometry_restraints.pair_proxies(
          flags=flags,
          bond_params_table=self.bond_params_table,
          shell_asu_tables=shell_asu_tables,
          model_indices=self.model_indices,
          conformer_indices=self.conformer_indices,
          sym_excl_indices=self.sym_excl_indices,
          donor_acceptor_excl_groups=self.donor_acceptor_excl_groups,
          nonbonded_params=self.nonbonded_params,
          nonbonded_types=self.nonbonded_types,
          nonbonded_charges=self.nonbonded_charges,
          nonbonded_distance_cutoff_plus_buffer
            =current_nonbonded_distance_cutoff_plus_buffer,
          min_cubicle_edge=self.min_cubicle_edge)
        introspection.virtual_memory_info().update_max()
        if (self._pair_proxies.nonbonded_proxies is None):
          break
        max_vdw_dist = self._pair_proxies.nonbonded_proxies.max_vdw_distance
        if (max_vdw_dist <= 0):
          break
        r = self.nonbonded_function.residual(
          vdw_distance=max_vdw_dist, delta=max_vdw_dist*(1+1.e-6))
        if (r != 0):
          break
        if (self.nonbonded_distance_cutoff < max_vdw_dist):
          if (self.nonbonded_distance_cutoff_was_determined_automatically):
            raise AssertionError("Internal error.")
          raise Sorry(
            "nonbonded_distance_cutoff=%.6g is too small:"
            " max_vdw_distance=%.6g" % (
              self.nonbonded_distance_cutoff,
              max_vdw_dist))
        self.adjusted_nonbonded_distance_cutoff = max_vdw_dist
        self.effective_nonbonded_buffer \
          = current_nonbonded_distance_cutoff_plus_buffer \
          - self.adjusted_nonbonded_distance_cutoff
        if (self.effective_nonbonded_buffer > bonded_distance_cutoff_epsilon):
          break
        self.adjusted_nonbonded_distance_cutoff \
          = self.nonbonded_distance_cutoff
    elif (self._pair_proxies is None):
      raise AssertionError("pair_proxies not defined already.")
    return self._pair_proxies

  def nonbonded_model_distances(self, sites_cart=None):
    pair_proxies = self.pair_proxies(sites_cart=sites_cart)
    if (sites_cart is None):
      sites_cart = self._sites_cart_used_for_pair_proxies
    return pair_proxies.nonbonded_proxies.deltas(sites_cart=sites_cart)

  def update_plain_pair_sym_table(self, sites_frac):
    asu_mappings = crystal.symmetry.asu_mappings(self.crystal_symmetry,
      buffer_thickness=self.plain_pairs_radius)
    asu_mappings.process_sites_frac(
      original_sites=sites_frac,
      site_symmetry_table=self.site_symmetry_table)
    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    pair_asu_table.add_all_pairs(distance_cutoff=self.plain_pairs_radius)
    self.plain_pair_sym_table=pair_asu_table.extract_pair_sym_table()
    introspection.virtual_memory_info().update_max()

  def energies_sites(self,
        sites_cart,
        flags=None,
        custom_nonbonded_function=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False,
        normalization=False,
        external_energy_function=None,
        extension_objects=[],
        site_labels=None):
    if(external_energy_function is not None):
      assert self.external_energy_function is None
    else:
      external_energy_function = self.external_energy_function
    if (flags is None):
      flags = geometry_restraints.flags.flags(default=True)
    pair_proxies = self.pair_proxies(
      flags=flags, sites_cart=sites_cart, site_labels=site_labels)
    (bond_proxies,
     nonbonded_proxies,
     nonbonded_function,
     angle_proxies,
     dihedral_proxies,
     reference_dihedral_proxies,
     ncs_dihedral_proxies,
     chirality_proxies,
     planarity_proxies,
     parallelity_proxies,
     generic_restraints) = [None]*11
    if (flags.bond):
      assert pair_proxies.bond_proxies is not None
      bond_proxies = pair_proxies.bond_proxies
    if (flags.nonbonded and self.nonbonded_types is not None):
      assert pair_proxies.nonbonded_proxies is not None
      nonbonded_proxies = pair_proxies.nonbonded_proxies
      if (custom_nonbonded_function is None):
        nonbonded_function = self.nonbonded_function
      else:
        nonbonded_function = custom_nonbonded_function
    if (flags.angle):     angle_proxies = self.angle_proxies
    if (flags.dihedral):  dihedral_proxies = self.dihedral_proxies
    if (flags.reference_dihedral):  \
      reference_dihedral_proxies = self.reference_dihedral_proxies
    if (flags.ncs_dihedral): ncs_dihedral_proxies = self.ncs_dihedral_proxies
    if (flags.chirality): chirality_proxies = self.chirality_proxies
    if (flags.planarity): planarity_proxies = self.planarity_proxies
    if (flags.parallelity): parallelity_proxies = self.parallelity_proxies
    if (flags.generic_restraints) :
      generic_restraints = self.generic_restraints_manager
    return geometry_restraints.energies.energies(
      sites_cart=sites_cart,
      bond_proxies=bond_proxies,
      nonbonded_proxies=nonbonded_proxies,
      nonbonded_function=nonbonded_function,
      angle_proxies=angle_proxies,
      dihedral_proxies=dihedral_proxies,
      reference_dihedral_proxies=reference_dihedral_proxies,
      ncs_dihedral_proxies=ncs_dihedral_proxies,
      chirality_proxies=chirality_proxies,
      planarity_proxies=planarity_proxies,
      parallelity_proxies=parallelity_proxies,
      generic_restraints_manager=generic_restraints,
      external_energy_function=external_energy_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      extension_objects=extension_objects,
      hbonds_in_bond_list=self.hbonds_in_bond_list)

  def harmonic_restraints(self, variables, type_indices, type_weights):
    assert self.shell_sym_tables is not None
    assert len(self.shell_sym_tables) > 0
    assert variables.size() == self.shell_sym_tables[0].size()
    residual_sum = 0
    gradients = flex.double(variables.size(), 0)
    for pair in self.shell_sym_tables[0].iterator():
      i,j = pair.i_seqs()
      if (type_indices is None):
        weight = type_weights
      else:
        weight = (  type_weights[type_indices[i]]
                  + type_weights[type_indices[j]]) * 0.5
      delta = variables[i] - variables[j]
      term = weight * delta
      residual_sum += term * delta
      gradients[i] += term * 2
      gradients[j] -= term * 2
    return store(residual_sum=residual_sum, gradients=gradients)

  def ta_harmonic_restraints(self, sites_cart, ta_harmonic_restraint_info, weight = 0.001, slack = 0.5):

    def delta(site1,
              site2):
      delta = math.sqrt(((site1[0]-site2[0])**2) + ((site1[1]-site2[1])**2) + ((site1[2]-site2[2])**2))
      if slack > 0:
        if (delta > self.ta_slack):
          delta_slack = delta - slack
        else:
          delta_slack = 0.0
      else: delta_slack = delta
      return delta, delta_slack

    def residual(distance):
      residual = self.ta_harmonic_weight * distance**2
      return residual

    def gradient(site,
                 ref_site,
                 distance_slack,
                 distance):
      site_delta = ((site[0]-ref_site[0]),(site[1]-ref_site[1]),(site[2]-ref_site[2]))
      if self.ta_slack > 0.0:
        if distance < self.ta_slack:
          return (0.0,0.0,0.0)
      gradient = (self.ta_harmonic_weight * 2.0 * distance_slack * site_delta[0],
                  self.ta_harmonic_weight * 2.0 * distance_slack * site_delta[1],
                  self.ta_harmonic_weight * 2.0 * distance_slack * site_delta[2])
      return gradient

    self.ta_harmonic_weight = weight
    self.ta_slack = slack
    gradients = flex.vec3_double(sites_cart.size(), (0,0,0))
    residuals = flex.double(sites_cart.size(), 0)
    for x in ta_harmonic_restraint_info:
      distance, distance_slack = delta(x[1],sites_cart[x[0]])
      residuals[x[0]] = residual(distance = distance)
      gradients[x[0]] = gradient(distance       = distance,
                                 distance_slack = distance_slack,
                                 site           = sites_cart[x[0]],
                                 ref_site       = x[1])
    return gradients

  def update_atom_nonbonded_type (self,
        i_seq,
        nonbonded_type,
        charge=0) :
    if (self.nonbonded_types is not None) :
      self.nonbonded_types[i_seq] = nonbonded_type
    if (self.nonbonded_charges is not None) :
      self.nonbonded_charges[i_seq] = charge

  def show_interactions(self,
        flags=None,
        sites_cart=None,
        site_labels=None,
        i_seq=None,
        f=None):
    if (f is None): f = sys.stdout
    pair_proxies = self.pair_proxies(flags=flags, sites_cart=sites_cart)
    if (sites_cart is None):
      sites_cart = self._sites_cart_used_for_pair_proxies
    if (pair_proxies.bond_proxies is not None):
      for proxy in pair_proxies.bond_proxies.simple:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "bond simple:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  distance_model: %.6g" % geometry_restraints.bond(
              sites_cart=sites_cart, proxy=proxy).distance_model
          print >> f, "  distance_ideal: %.6g" % proxy.distance_ideal
          print >> f, "  weight: %.6g" % proxy.weight
          if (proxy.slack > 0):
            print >> f, "  slack: %.6g" % proxy.slack
      if (pair_proxies.bond_proxies.asu.size() > 0):
        asu_mappings = pair_proxies.bond_proxies.asu_mappings()
        for proxy in pair_proxies.bond_proxies.asu:
          if (i_seq is None or i_seq in [proxy.i_seq, proxy.j_seq]):
            rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
            if (site_labels is None):
              print >> f, "bond asu:", (proxy.i_seq, proxy.j_seq), rt_mx
            else:
              print >> f, "bond asu:", (proxy.i_seq, proxy.j_seq)
              print >> f, " ", site_labels[proxy.i_seq]
              print >> f, " ", site_labels[proxy.j_seq], rt_mx
            if (sites_cart is not None):
              print >> f, "  distance_model: %.6g" % geometry_restraints.bond(
                sites_cart=sites_cart,
                asu_mappings=asu_mappings,
                proxy=proxy).distance_model
            print >> f, "  distance_ideal: %.6g" % proxy.distance_ideal
            print >> f, "  weight: %.6g" % proxy.weight
            if (proxy.slack > 0):
              print >> f, "  slack: %.6g" % proxy.slack
    if (self.angle_proxies is not None):
      for proxy in self.angle_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "angle:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  angle_model: %.6g" % geometry_restraints.angle(
              sites_cart=sites_cart, proxy=proxy).angle_model
          print >> f, "  angle_ideal: %.6g" % proxy.angle_ideal
          print >> f, "  weight: %.6g" % proxy.weight
    if (self.dihedral_proxies is not None):
      for proxy in self.dihedral_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "dihedral:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  angle_model: %.6g" % geometry_restraints.dihedral(
              sites_cart=sites_cart, proxy=proxy).angle_model
          print >> f, "  angle_ideal: %.6g" % proxy.angle_ideal
          print >> f, "  weight: %.6g" % proxy.weight
          print >> f, "  periodicity: %.6g" % proxy.periodicity
    if (self.reference_dihedral_proxies is not None):
      for proxy in self.reference_dihedral_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "reference dihedral:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  angle_model: %.6g" % geometry_restraints.dihedral(
              sites_cart=sites_cart, proxy=proxy).angle_model
          print >> f, "  angle_ideal: %.6g" % proxy.angle_ideal
          print >> f, "  weight: %.6g" % proxy.weight
          print >> f, "  limit: %.6g" % proxy.limit
    if (self.ncs_dihedral_proxies is not None):
      for proxy in self.ncs_dihedral_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "NCS dihedral:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  angle_model: %.6g" % geometry_restraints.dihedral(
              sites_cart=sites_cart, proxy=proxy).angle_model
          print >> f, "  angle_ideal: %.6g" % proxy.angle_ideal
          print >> f, "  weight: %.6g" % proxy.weight
          print >> f, "  limit: %.6g" % proxy.limit
    if (self.chirality_proxies is not None):
      for proxy in self.chirality_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "chirality:", proxy.i_seqs
          if (site_labels is not None):
            for i in proxy.i_seqs:
              print >> f, " ", site_labels[i]
          if (sites_cart is not None):
            print >> f, "  volume_model: %.6g" % geometry_restraints.chirality(
              sites_cart=sites_cart, proxy=proxy).volume_model
          print >> f, "  volume_ideal: %.6g" % proxy.volume_ideal
          print >> f, "  both_signs: %.6g" % proxy.both_signs
          print >> f, "  weight: %.6g" % proxy.weight
    if (self.planarity_proxies is not None):
      for proxy in self.planarity_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "planarity:", tuple(proxy.i_seqs)
          if (sites_cart is not None):
            deltas = geometry_restraints.planarity(
              sites_cart=sites_cart, proxy=proxy).deltas()
          else:
            deltas = [None]*proxy.i_seqs.size()
          for i,weight,delta in zip(proxy.i_seqs, proxy.weights, deltas):
            print >> f, " ",
            if (site_labels is not None):
              print >> f, site_labels[i],
            if (delta is not None):
              print >> f, "delta: %8.5f," % delta,
            print >> f, "weight: %.6g" % weight
    if (self.parallelity_proxies is not None):
      print >> f, "parallelity functional is not implemented"
    if (pair_proxies.nonbonded_proxies is not None):
      simple = pair_proxies.nonbonded_proxies.simple
      asu = pair_proxies.nonbonded_proxies.asu
      simple_size = simple.size()
      if (asu.size() > 0):
        asu_mappings = pair_proxies.nonbonded_proxies.asu_mappings()
      if (sites_cart is not None):
        deltas = geometry_restraints.nonbonded_deltas(
          sites_cart=sites_cart,
          sorted_asu_proxies=pair_proxies.nonbonded_proxies)
        permutation = flex.sort_permutation(data=deltas)
      else:
        deltas = None
        permutation = xrange(simple_size + asu.size())
      for i_proxy in permutation:
        if (i_proxy < simple_size):
          proxy = simple[i_proxy]
          if (i_seq is None or i_seq in proxy.i_seqs):
            print >> f, "nonbonded simple:", proxy.i_seqs
            if (site_labels is not None):
              for i in proxy.i_seqs:
                print >> f, " ", site_labels[i]
            if (deltas is not None):
              print >> f, "  distance_model: %.6g" % deltas[i_proxy]
            print >> f, "  vdw_distance: %.6g" % proxy.vdw_distance
        else:
          proxy = asu[i_proxy-simple_size]
          if (i_seq is None or i_seq in [proxy.i_seq, proxy.j_seq]):
            rt_mx = asu_mappings.get_rt_mx_ji(pair=proxy)
            if (site_labels is None):
              print >> f, "nonbonded asu:", (proxy.i_seq, proxy.j_seq), rt_mx
            else:
              print >> f, "nonbonded asu:", (proxy.i_seq, proxy.j_seq)
              print >> f, " ", site_labels[proxy.i_seq]
              print >> f, " ", site_labels[proxy.j_seq], rt_mx
            if (deltas is not None):
              print >> f, "  distance_model: %.6g" % deltas[i_proxy]
            print >> f, "  vdw_distance: %.6g" % proxy.vdw_distance

  def show_sorted(self,
        flags=None,
        sites_cart=None,
        site_labels=None,
        f=None):
    if (f is None): f = sys.stdout
    pair_proxies = self.pair_proxies(flags=flags, sites_cart=sites_cart)
    if (sites_cart is None):
      sites_cart = self._sites_cart_used_for_pair_proxies
    if pair_proxies.bond_proxies is not None:
      pair_proxies.bond_proxies.show_sorted(
          by_value="residual",
          sites_cart=sites_cart,
          site_labels=site_labels,
          f=f,
          exclude=self.hbonds_in_bond_list)
      print >> f
      tempbuffer = StringIO.StringIO()
      if self.hbonds_in_bond_list is not None:
        pair_proxies.bond_proxies.show_sorted(
            by_value="residual",
            sites_cart=sites_cart,
            site_labels=site_labels,
            f=tempbuffer,
            prefix="",
            exclude=self.hbonds_in_bond_list,
            exclude_output_only=True)
        print >> f, "Bond-like", tempbuffer.getvalue()[5:]
    if (self.angle_proxies is not None):
      self.angle_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
    if (self.dihedral_proxies is not None):
      self.dihedral_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
    if (self.reference_dihedral_proxies is not None):
      self.reference_dihedral_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f, is_reference=True)
      print >> f
    if (self.ncs_dihedral_proxies is not None):
      self.ncs_dihedral_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f, is_ncs=True)
      print >> f
    if (self.generic_restraints_manager is not None):
      if (self.generic_restraints_manager.c_beta_dihedral_proxies is not None):
        self.generic_restraints_manager.c_beta_dihedral_proxies.\
          show_sorted(
            by_value="residual",
            sites_cart=sites_cart,
            site_labels=site_labels,
            f=f, is_c_beta=True)
        print >> f
      if (self.generic_restraints_manager.hydrogen_bond_proxies is not None):
        self.generic_restraints_manager.\
            show_sorted_hbonds(
                by_value="residual",
                sites_cart=sites_cart,
                site_labels=site_labels,
                f=f)
      self.generic_restraints_manager.\
          show_sorted_ramachandran(
              by_value="residual",
              sites_cart=sites_cart,
              site_labels=site_labels,
              f=f)
    if (self.chirality_proxies is not None):
      self.chirality_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
    if (self.planarity_proxies is not None):
      self.planarity_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
    if (self.parallelity_proxies is not None):
      self.parallelity_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
    if (pair_proxies.nonbonded_proxies is not None):
      pair_proxies.nonbonded_proxies.show_sorted(
        by_value="delta",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f

  def get_nonbonded_clashscore(
    self,
    sites_cart,
    hd_sel,
    site_labels=None):
    '''(multiple arguments) -> nonbonded_clashscore obj

    Construct nonbonded_clash_info, the non-bonded clashss lists and scores

    Arguments:
    sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
    site_labels: a list of lables such as " HA  LEU A  38 ", for each atom
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not

    NOTE:
    As of Dec. 2013 manipulation of scatterers can produce scatteres which
    have no lables. In that case, water interaction score will not be accurate

    Return:
       self.nonbonded_clash_info

    NOTE:
    Be aware to the parameters:
       assume_hydrogens_all_missing=False,
       hard_minimum_nonbonded_distance=0.0

    The defult values of these are True and 0.001, which will alter
    the size of the vdw bonds and the clashes that being counted

    As of Dec. 2013 manipulation of scatterers can produce scatteres whish
    have no lables. In that case, water interaction score will not be accurate


    Usage example:
    mon_lib_srv = mmtbx.monomer_library.server.server()
    ener_lib = mmtbx.monomer_library.server.ener_lib()

    pdb = monomer_library.pdb_interpretation.process(
      file_name=None,
      raw_records=records,
      substitute_non_crystallographic_unit_cell_if_necessary=True,
      mon_lib_srv=mon_lib_srv,
      ener_lib=ener_lib)

    grm = pdb.geometry_restraints_manager(
      assume_hydrogens_all_missing=False,
      hard_minimum_nonbonded_distance=0.0))

    xrs = pdb.xray_structure()
    sites_cart = xrs.sites_cart()
    site_labels = xrs.scatterers().extract_labels()
    hd_sel = xrs.hd_selection()

    nb_clash_info = grm.get_nonbonded_clashscore(
      sites_cart=sites_cart,
      site_labels=site_labels,
      hd_sel=hd_sel)


    @author: Youval Dar, LBL 12-2013
    '''
    pair_proxies = self.pair_proxies(
          sites_cart=sites_cart,
          site_labels=site_labels)

    if pair_proxies.nonbonded_proxies.n_unknown_nonbonded_type_pairs != 0:
      msg = "nonbonded clashscore can't be calculated."
      msg += " PDB file contains unknown type pairs. Please provide cif file."
      raise Sorry(msg)

    proxies_info_nonbonded = pair_proxies.nonbonded_proxies.get_sorted(
      by_value="delta",
      sites_cart=sites_cart,
      site_labels=site_labels)

    if proxies_info_nonbonded != None:
      nonbonded_list = proxies_info_nonbonded[0]
    else:
      nonbonded_list = []

    nonbonded_clash_info = nonbonded_clashscore(
      nonbonded_list=nonbonded_list,
      hd_sel=hd_sel,
      full_connectivty_table=self.shell_sym_tables[0].full_simple_connectivity(),
      connectivty_table_2=self.shell_sym_tables[2].full_simple_connectivity(),
      sites_cart=sites_cart)
    return nonbonded_clash_info

def construct_non_crystallographic_conserving_bonds_and_angles(
      sites_cart,
      edge_list_bonds,
      edge_list_angles,
      bond_weight=100,
      angle_weight=50,
      vdw_distance=1.2,
      non_crystallographic_unit_cell_buffer_layer=5,
      asu_mappings_buffer_thickness=5,
      max_reasonable_bond_distance=10):
  import cctbx.crystal.coordination_sequences
  from cctbx import uctbx
  from scitbx import matrix
  bond_proxies = geometry_restraints.bond_sorted_asu_proxies(
    asu_mappings=None)
  for edge_list,weight in [(edge_list_bonds, bond_weight),
                           (edge_list_angles, angle_weight)]:
    for i,j in edge_list:
      distance = abs(matrix.col(sites_cart[i]) - matrix.col(sites_cart[j]))
      bond_proxies.process(geometry_restraints.bond_simple_proxy(
        i_seqs=(i,j), distance_ideal=distance, weight=weight))
  bond_params_table = geometry_restraints.extract_bond_params(
    n_seq=sites_cart.size(),
    bond_simple_proxies=bond_proxies.simple)
  box = uctbx.non_crystallographic_unit_cell_with_the_sites_in_its_center(
    sites_cart=sites_cart,
    buffer_layer=non_crystallographic_unit_cell_buffer_layer)
  asu_mappings = box.crystal_symmetry().special_position_settings() \
    .asu_mappings(
      buffer_thickness=asu_mappings_buffer_thickness,
      sites_cart=box.sites_cart)
  bond_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
  geometry_restraints.add_pairs(bond_asu_table, bond_proxies.simple)
  shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
    pair_asu_table=bond_asu_table,
    max_shell=3)
  shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
    for shell_asu_table in shell_asu_tables]
  nonbonded_types = flex.std_string(bond_params_table.size(), "Default")
  nonbonded_params = geometry_restraints.nonbonded_params()
  nonbonded_params.distance_table.setdefault(
    "Default")["Default"] = vdw_distance
  return box.sites_cart, manager(
    crystal_symmetry=box.crystal_symmetry(),
    site_symmetry_table=asu_mappings.site_symmetry_table(),
    bond_params_table=bond_params_table,
    shell_sym_tables=shell_sym_tables,
    nonbonded_params=nonbonded_params,
    nonbonded_types=nonbonded_types,
    nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
    max_reasonable_bond_distance=max_reasonable_bond_distance)

def format_distances_for_error_message(
      pair_sym_table,
      larger_than,
      orthogonalization_matrix,
      sites_frac,
      sites_cart,
      site_labels):
  assert [sites_frac, sites_cart].count(None) == 1
  if (site_labels is not None):
    if (sites_frac is not None):
      assert len(site_labels) == len(sites_frac)
    else:
      assert len(site_labels) == len(sites_cart)
  result = []
  from scitbx.matrix import col, sqr
  if (orthogonalization_matrix is not None):
    orthogonalization_matrix = sqr(orthogonalization_matrix)
  for i_seq,pair_sym_dict in enumerate(pair_sym_table):
    for j_seq,sym_ops in pair_sym_dict.items():
      for rt_mx_ji in sym_ops:
        if (sites_cart is not None):
          assert rt_mx_ji.is_unit_mx()
          d_cart = col(sites_cart[i_seq]) - col(sites_cart[j_seq])
        else:
          d_frac = col(sites_frac[i_seq]) - col(rt_mx_ji * sites_frac[j_seq])
          d_cart = orthogonalization_matrix * d_frac
        dist = d_cart.length()
        if (dist > larger_than):
          if (rt_mx_ji.is_unit_mx()):
            ss = ""
          else:
            ss = " " + str(rt_mx_ji)
          if (site_labels is None):
            si, sj = str(i_seq), str(j_seq)
          else:
            si, sj = site_labels[i_seq], site_labels[j_seq]
          if (dist < 1000):
            sd = "%7.3f" % dist
          else:
            sd = "%.6g" % dist
          result.append("distance: %s - %s: %s%s" % (si, sj, sd, ss))
  return result

class nonbonded_clashscore(object):
  """
  Object processing information on non-bonded steric clashes. When the
  distance between atoms is smaller than the Van Der Waals distance.

  @author: Youval Dar (LBL 2013)
  """

  def __init__(
    self,
    nonbonded_list,
    hd_sel,
    full_connectivty_table,
    connectivty_table_2,
    sites_cart):
    """
    Arguments:
    nonbonded_list: a list with items in the following format
                    ([pdb labels], i_seq, j_seq, model, vdw_distance, sym_op_j, rt_mx)
                    i_seq,j_seq: position of residues in the pdb file
                    model: The pdb or other model distance
                    vdw_distance: Van Der Waals distance
                    sym_op_j: is this a product of a symmetry opeation
                    rt_mx: Rotation matrix for symmetry operation
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not
    full_connectivty_table: full_connectivty_table[i] is a dictinary constaining a
                            list of all atoms connected to atom i
    connectivty_table_2: connectivty_table[i] is a dictinary constaining a
                         list of all atoms connected to atom i
    sites_cart: sites_cart[i] tuple containing the x,y,z coordinates of atom i
    """
    self.nonbonded_list = nonbonded_list
    self.hd_sel = hd_sel
    self.full_connectivty_table = full_connectivty_table
    self.connectivty_table_2 = connectivty_table_2
    self.sites_cart = sites_cart
    if nonbonded_list != []:
      try:
        clashlist = self.clashscore_clash_list()
      except TypeError as e:
        # When proxies_info_nonbonded are not available
        if e == "vec3_double' object is not callable":
          raise Sorry(e)
        else:
          clashlist = [[],[],[]]
          print e
      except Sorry as e:
        raise Sorry(e)
      except Exception as e:
        raise Sorry('Unexpected error processing proxies_info_nonbonded in clashscore_clash_list()')
      # Collect, seperatly, clashes due to symmetry operation and those that are not
      self.nb_clash_proxies_due_to_sym_op = clashlist[0]
      self.nb_clash_proxies_solvent_solvent = clashlist[1]
      self.nb_clash_proxies_simple = clashlist[2]
      self.nb_clash_proxies_all_clashes = clashlist[0] + clashlist[1] + clashlist[2]
      # clashscore is the number of serious steric overlaps (>0.4A) per 1000 atoms
      n_atoms = len(self.sites_cart)
      clashscore_due_to_sym_op = len(clashlist[0])*1000/n_atoms
      clashscore_solvent_solvent = len(clashlist[1])*1000/n_atoms
      clashscore_simple = len(clashlist[2])*1000/n_atoms
      clashscore_all_clashes = clashscore_due_to_sym_op + clashscore_solvent_solvent + clashscore_simple
      #
      self.nb_clashscore_due_to_sym_op = clashscore_due_to_sym_op
      self.nb_clashscore_solvent_solvent = clashscore_solvent_solvent
      self.nb_clashscore_simple = clashscore_simple
      self.nb_clashscore_all_clashes = clashscore_all_clashes
    else:
      self.nb_clashscore_due_to_sym_op = 0
      self.nb_clashscore_solvent_solvent = 0
      self.nb_clashscore_simple = 0
      self.nb_clashscore_all_clashes = 0
      #
      self.nb_clash_proxies_due_to_sym_op = []
      self.nb_clash_proxies_solvent_solvent = []
      self.nb_clash_proxies_simple = []
      self.nb_clash_proxies_all_clashes = []

  @staticmethod
  def is_1_5_interaction(i_seq,j_seq,hd_sel,full_connectivty_table):
    '''(int,int,bool array,list of lists of int) -> bool
    Check if we have 1-5 interaction between two hydrogen and heavy atom

    arguments:
    i_seq,j_seq:  are the number of the atoms we are checking in the pdb table
    hd_sel: hd_sel[i] retruns True of False, indicating whether an atom i is a Hydrogen or not
    full_connectivty_table: full_connectivty_table[i] is a dictinary constaining a
                            list of all atoms connected to atom

    retruns:
    True if we have 1-5 hydrogen and heavy atom interaction
    '''
    # check if we have hydrogen - heavy atom interaction
    xor = lambda a,b: (a or b) and not (a and b)
    if xor(hd_sel[i_seq],hd_sel[j_seq]):
      # starting with hydrogen will make process shorter
      if not hd_sel[i_seq]:
        i_seq,j_seq = j_seq,i_seq
      # build connection table of i_seq, 4 steps deep
      atoms_numbers = dict([(i_seq, 0)]) # i_seq is in zero distance
      used_connections = {i_seq}
      new_connections = {i_seq}
      for i in range(2,6):
        connections = set()
        for key in new_connections:
          # add all new connetions for the current step
          connections = connections.union(set(full_connectivty_table[key]))
        # Remove the connection that were already used
        new_connections = connections - used_connections
        # Add the new connection to the used once
        used_connections = used_connections.union(connections)
        # Add the new atoms with their distance from key
        for new_atom in new_connections:
          atoms_numbers[new_atom] = i
      # return true if j_seq in the is 1-5 connection
      return (j_seq in atoms_numbers) and (atoms_numbers[j_seq] == 5)
    else:
      return False

  def get_clash_keys(self,record):
    '''(str) -> str,str,str

    Collect atoms keys and coordinates

    Argumets:
    record : a clash record, for example:
    (['pdb=" HA  LEU A  38 "', 'pdb="HD23 LEU A  38 "'], 523, 532, 1.8108674716831155, 2.44, '', None)

    Returns:
    clash_vec: a unique key for a clash. for example: '0.144,0.323,1.776,0.154,0.327,1.786'
    key1,key2: are the atoms keys. for example ' HA  LEU A  38 ::120::' or ' CB ASN A  55 ::52::sym.op.'
               "atom1 label::atom1 number::sym.op""
    '''
    # get atoms keys
    record = list(record)
    key = '{0}::{1}::{2}'
    key1 = key.format(record[0][0][5:-1],record[1],record[5])
    key2 = key.format(record[0][1][5:-1],record[2],record[5])
    # get coordinates
    atomi_xyz = self.sites_cart[record[1]]
    atomj_xyz = self.sites_cart[record[2]]
    # make tupe containing both atom coordinates
    vec = atomi_xyz + atomj_xyz
    # creat a string from vec
    clash_vec = '{0:.3f},{1:.3f},{2:.3f},{3:.3f},{4:.3f},{5:.3f}'.format(*vec)
    # make key1 always smaller than key2
    if key1 > key2:
      key1,key2 = key2,key1
    return clash_vec,key1,key2

  @staticmethod
  def cos_vec(u,v):
    '''(tuple,tuple) -> float

    Calculate the cosine used to evaluate wheather the atoms
    coliding are inline or not

    Arguments:
    u,v: lists containing x1,y1,z1,x2,y2,z2 coordinates
    Returns:
    cos_angle: the cosine of the angle between center of the common atom
    and the mid point between the other two atoms.
    '''
    # Calculate dot product
    u_dot_v = lambda u,v: (u[0]*v[0]+u[1]*v[1]+u[2]*v[2])
    # Calculate mid point between to atoms
    u_mid_v = lambda u,v: [(x+y)/2 for (x,y) in zip(u,v)]
    # find the length of a vector
    u_len = lambda u: math.sqrt(u[0]**2 + u[1]**2 + u[2]**2)
    # Move vector to the origin
    make_vec = lambda u1,u2: [(x-y) for (x,y) in zip(u1,u2)]
    v1 = v[0:3]
    v2 = v[3:6]
    u1 = u[0:3]
    u2 = u[3:6]
    # find the common atoms
    if v1 == u1:
      mid_uv = u_mid_v(u2,v2)
      vec1 = make_vec(v1,mid_uv)
      vec2 = make_vec(u2,v2)
    elif v1 == u2:
      mid_uv = u_mid_v(u1,v2)
      vec1 = make_vec(v1,mid_uv)
      vec2 = make_vec(u1,v2)
    elif v2 == u1:
      mid_uv = u_mid_v(u2,v1)
      vec1 = make_vec(v2,mid_uv)
      vec2 = make_vec(u2,v1)
    else: #v2 == u2
      mid_uv = u_mid_v(u1,v1)
      vec1 = make_vec(v2,mid_uv)
      vec2 = make_vec(u1,v1)
    # return the cosine of the angle between the two vectors
    try:
      cos_angle = abs(u_dot_v(vec1,vec2)/u_len(vec1)/u_len(vec2))
    except ZeroDivisionError:
      cos_angle = 1
    return cos_angle

  @staticmethod
  def are_both_solvent(key,connectivty_table_2):
    '''(str,dict) -> bool
    retrun true if the two atoms, which number is specified in "key" are single, double or tripple atom molecule

    Arguments:
    connectivty_table_2: connectivty_table[i] is a dictinary constaining a
                         list of all atoms connected to atom i
    key: a list "atom1 lable::atom1 number::sym.op.::atom2 label::atom2 number::sym.op."

    Return True is both atoms have no connection to atoms two bonds apart
    '''
    key = key.split('::')
    i = int(key[1])
    j = int(key[4])
    i_is_solvent = list(connectivty_table_2[i])==[]
    j_is_solvent = list(connectivty_table_2[j])==[]
    return i_is_solvent and j_is_solvent

  def clashscore_clash_list(self):
    '''(self) -> [list,list,list]
    Collect inforamtion of clashing nonbonded proxies (neighboring atoms) for clash
    score calculation.
    - proxies considered to clash when model_distance - van_der_waals_distance < -0.41
    - For every atom consider only worst clash
    - Do not double count in-line clashes.
       for example, C-H1  H2-...
       if both C and H1 are clashing with H2 count only the worst of them
    - Exclude 1-5 interaction of Hydrogen and heavy atom

    Retruns:
    clashing_atoms_list: [[list of clashes due to sym operation],[list of solvent-solvent clashes],[list of all other clashes]]
    '''
    clashing_atoms_list = [[],[],[]]
    clashing_atoms_dict = {}
    # clashing_atoms_dict[atom_key] = [all clashes information for this atom]
    # clashing_atoms_dict[' HA  ASP A  44 '] = [[atom1,vec1, i_seq, j_seq1,model)],
    #                                           [atom2, vec2,i_seq, j_seq2,model]...]
    # vec = [delta_x,delta_y,delta_z] , the difference between the atoms coordinates
    clashes_dict = {}
    # clashes_dict[vec] = [atom1, atom2, clash_record]
    # vec uniquly define a clash, regardless of atoms order
    for rec in self.nonbonded_list:
      i_seq = rec[1]
      j_seq = rec[2]
      model = rec[3]
      vdw = rec[4]
      symop = rec[5]
      delta = model - vdw
      # check for clash
      if (delta < -0.41):
        # Check of 1-5 interaction
        if not self.is_1_5_interaction(i_seq, j_seq,self.hd_sel,self.full_connectivty_table):
          clash_vec,key1,key2 = self.get_clash_keys(rec)
          # clash_key is a string a string that uniquly identify each clash
          if key1>key2:
            clash_key = '::'.join([key2,key1])
          else:
            clash_key = '::'.join([key1,key2])
          if clash_key not in clashes_dict:
            # record new clash
            clashes_dict[clash_key] = [key1,key2,rec]
            if key1 not in clashing_atoms_dict: clashing_atoms_dict[key1] = []
            if key2 not in clashing_atoms_dict: clashing_atoms_dict[key2] = []
            clashing_atoms_dict[key1].append([key2,clash_vec,clash_key,symop,model])
            clashing_atoms_dict[key2].append([key1,clash_vec,clash_key,symop,model])
    for key in clashing_atoms_dict:
      # for atoms that clash more than once, check for inline clashes
      if len(clashing_atoms_dict[key]) > 1:
        temp_clash_list = []
        n_clashes = len(clashing_atoms_dict[key])
        # itereate over all clash combination
        for i in range(n_clashes-1):
          for j in range(i+1,n_clashes):
            vec_i = clashing_atoms_dict[key][i][1]
            vec_j = clashing_atoms_dict[key][j][1]
            u = map(float,vec_i.split(','))
            v = map(float,vec_j.split(','))
            cos_angle = 0
            # test inline only if the two atoms, clashing with the
            # common atom, are connected
            clashing_atom_1 = int(clashing_atoms_dict[key][i][0].split('::')[1])
            clashing_atom_2 = int(clashing_atoms_dict[key][j][0].split('::')[1])
            if clashing_atom_1 in self.full_connectivty_table[clashing_atom_2]:
              if not [0.0,0.0,0.0] in [u,v]:
                # Do not calculate cos_angle when atoms are overlapping
                if clashing_atoms_dict[key][i][3] == '':
                  # ignore clashes that are cause by symmetry operation
                  cos_angle = self.cos_vec(u,v)
            # atoms consider to be inline if cosine of the angle between vectors > 0.866
            if cos_angle > 0.867 and (vec_i != vec_j):
              # for inline clashes keep the closer two(compare models)
              if clashing_atoms_dict[key][i][4] < clashing_atoms_dict[key][j][4]:
                temp_clash_list.append(clashing_atoms_dict[key][i])
                # remove clash from clashes_dict
                remove_key = clashing_atoms_dict[key][j][2]
                if remove_key in clashes_dict: del clashes_dict[remove_key]
              else:
                temp_clash_list.append(clashing_atoms_dict[key][j])
                # remove clash from clashes_dict
                remove_key = clashing_atoms_dict[key][i][2]
                if remove_key in clashes_dict: del clashes_dict[remove_key]
            else:
              # clashes are not inline, keep both
              temp_clash_list.append(clashing_atoms_dict[key][j])
              temp_clash_list.append(clashing_atoms_dict[key][i])
        clashing_atoms_dict[key] = temp_clash_list
    for (key,val) in clashes_dict.iteritems():
      if key.split('::')[2] != '':
        # not to symmetry operation
        clashing_atoms_list[0].append(val[2])
      elif self.are_both_solvent(key,self.connectivty_table_2):
        # solvent-solvent clashes
        clashing_atoms_list[1].append(val[2])
      else:
        # not due to symmetry operation or solvent-solvent, simple clashes
        clashing_atoms_list[2].append(val[2])
    return clashing_atoms_list

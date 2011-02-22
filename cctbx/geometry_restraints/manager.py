from cctbx import geometry_restraints
import cctbx.geometry_restraints.flags
import cctbx.geometry_restraints.energies
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import store
from libtbx import introspection
from libtbx import adopt_init_args
from libtbx import dict_with_default_0
import sys, math

class manager(object):

  def __init__(self,
        crystal_symmetry=None,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        site_symmetry_table=None,
        bond_params_table=None,
        shell_sym_tables=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_function=None,
        nonbonded_distance_cutoff=None,
        nonbonded_buffer=1,
        angle_proxies=None,
        dihedral_proxies=None,
        reference_dihedral_proxies=None,
        chirality_proxies=None,
        planarity_proxies=None,
        generic_restraints_manager=None,
        external_energy_function=None,
        plain_pairs_radius=None,
        max_reasonable_bond_distance=None,
        min_cubicle_edge=5):
    if (site_symmetry_table is not None): assert crystal_symmetry is not None
    if (bond_params_table is not None and site_symmetry_table is not None):
      assert bond_params_table.size() == site_symmetry_table.indices().size()
    if (shell_sym_tables is not None and site_symmetry_table is not None):
      assert len(shell_sym_tables) > 0
      assert shell_sym_tables[0].size() == site_symmetry_table.indices().size()
    if (nonbonded_types is not None and site_symmetry_table is not None):
      assert nonbonded_types.size() == site_symmetry_table.indices().size()
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
        omit_bonds_with_slack_greater_than=0):
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
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      site_symmetry_table=self.site_symmetry_table,
      bond_params_table=self.bond_params_table,
      shell_sym_tables=reduced_shell_sym_tables,
      angle_proxies=reduced_angle_proxies,
      dihedral_proxies=reduced_dihedral_proxies)

  def sites_cart_used_for_pair_proxies(self):
    return self._sites_cart_used_for_pair_proxies

  def new_including_isolated_sites(self,
        n_additional_sites,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        site_symmetry_table=None,
        nonbonded_types=None):
    assert n_additional_sites >= 0
    assert (model_indices is None) == (self.model_indices is None)
    assert (conformer_indices is None) == (self.conformer_indices is None)
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
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      model_indices=model_indices,
      conformer_indices=conformer_indices,
      sym_excl_indices=sym_excl_indices,
      site_symmetry_table=site_symmetry_table,
      bond_params_table=bond_params_table,
      shell_sym_tables=shell_sym_tables,
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=nonbonded_types,
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      reference_dihedral_proxies=self.reference_dihedral_proxies,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      plain_pairs_radius=self.plain_pairs_radius)

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
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      model_indices=selected_model_indices,
      conformer_indices=selected_conformer_indices,
      sym_excl_indices=selected_sym_excl_indices,
      site_symmetry_table=selected_site_symmetry_table,
      bond_params_table=selected_bond_params_table,
      shell_sym_tables=selected_shell_sym_tables,
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=selected_nonbonded_types,
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=selected_angle_proxies,
      dihedral_proxies=selected_dihedral_proxies,
      reference_dihedral_proxies=selected_reference_dihedral_proxies,
      chirality_proxies=selected_chirality_proxies,
      planarity_proxies=selected_planarity_proxies,
      plain_pairs_radius=self.plain_pairs_radius)

  def remove_angles_in_place(self, selection):
    self.angle_proxies = self.angle_proxies.proxy_remove(
      selection=selection)

  def remove_dihedrals_in_place(self, selection):
    self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_reference_dihedrals_in_place(self, selection):
    self.reference_dihedral_proxies = self.reference_dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_chiralities_in_place(self, selection):
    self.chirality_proxies = self.chirality_proxies.proxy_remove(
      selection=selection)

  def remove_planarities_in_place(self, selection):
    self.planarity_proxies = self.planarity_proxies.proxy_remove(
      selection=selection)

  def set_generic_restraints (self, manager) :
    self.generic_restraints_manager = manager

  def set_external_energy_function (self, energy_function) :
    self.external_energy_function = energy_function

  def pair_proxies(self,
        sites_cart=None,
        flags=None,
        asu_is_inside_epsilon=None,
        bonded_distance_cutoff_epsilon=None):
    if (bonded_distance_cutoff_epsilon is None):
      bonded_distance_cutoff_epsilon = 1.e-6
    bonded_distance_cutoff = -1
    def check_bonded_distance_cutoff():
        if (    self.max_reasonable_bond_distance is not None
            and bonded_distance_cutoff > self.max_reasonable_bond_distance):
          raise RuntimeError(
            "Bond distance > max_reasonable_bond_distance: %.6g > %.6g" % (
              bonded_distance_cutoff, self.max_reasonable_bond_distance))
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
                orthogonalization_matrix=
                self.crystal_symmetry.unit_cell().orthogonalization_matrix(),
                sites_frac=sites_frac),
              default=0))
        check_bonded_distance_cutoff()
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
            check_bonded_distance_cutoff()
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
                    orthogonalization_matrix=
                      unit_cell.orthogonalization_matrix(),
                    sites_frac=sites_frac),
                  default=0))
            check_bonded_distance_cutoff()
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
          nonbonded_params=self.nonbonded_params,
          nonbonded_types=self.nonbonded_types,
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
          raise AssertionError(
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
        extension_objects=[]):
    if (flags is None):
      flags = geometry_restraints.flags.flags(default=True)
    pair_proxies = self.pair_proxies(flags=flags, sites_cart=sites_cart)
    (bond_proxies,
     nonbonded_proxies,
     nonbonded_function,
     angle_proxies,
     dihedral_proxies,
     reference_dihedral_proxies,
     chirality_proxies,
     planarity_proxies,
     generic_restraints) = [None]*9
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
    if (flags.reference_dihedral):  reference_dihedral_proxies = self.reference_dihedral_proxies
    if (flags.chirality): chirality_proxies = self.chirality_proxies
    if (flags.planarity): planarity_proxies = self.planarity_proxies
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
      chirality_proxies=chirality_proxies,
      planarity_proxies=planarity_proxies,
      generic_restraints_manager=generic_restraints,
      external_energy_function=self.external_energy_function,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization,
      extension_objects=extension_objects)

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
    if (pair_proxies.bond_proxies is not None):
      pair_proxies.bond_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f
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
    if (pair_proxies.nonbonded_proxies is not None):
      pair_proxies.nonbonded_proxies.show_sorted(
        by_value="delta",
        sites_cart=sites_cart, site_labels=site_labels, f=f)
      print >> f

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

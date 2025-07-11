from __future__ import absolute_import, division, print_function
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
import sys, math
from six import StringIO
import iotbx.pdb

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
bp.import_ext("scitbx_array_family_flex_ext")
from scitbx_array_family_flex_ext import reindexing_array

from cctbx.geometry_restraints.linking_class import linking_class
origin_ids = linking_class()
from cctbx.geometry_restraints.base_geometry import Base_geometry

def create_side_chain_restraints(hierarchy, grm, log=None):
  proxies=[]
  for rg in hierarchy.residue_groups():
    ags=rg.atom_groups()
    names = []
    for atom in rg.atoms():
      if atom.element_is_hydrogen(): continue
      if atom.name not in names:
        names.append(atom.name)
    pairs=[[],[]]
    for name in names:
      if name in [' CA ', ' N  ', ' C  ', ' O  ']: continue
      for i, ag in enumerate(ags):
        if not ag.altloc: continue
        pairs[i].append(ag.get_atom(name.strip()).i_seq)
    sigma=.2
    if not pairs[0] or not pairs[1]: continue
    if len(pairs[0])<3 or len(pairs[1])<3: continue
    proxy=geometry_restraints.parallelity_proxy(
              i_seqs=flex.size_t(pairs[0]),
              j_seqs=flex.size_t(pairs[1]),
              weight=1/(sigma**2),
              target_angle_deg=0,
              slack=0,
              top_out=False,
              limit=1,
              origin_id=origin_ids.get_origin_id('side-chain parallelity'))
    proxies.append(proxy)
  return proxies

class manager(Base_geometry):
# This class is documented in
# http://www.phenix-online.org/papers/iucrcompcomm_aug2004.pdf
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
        reference_coordinate_proxies=None,
        reference_dihedral_manager=None,
        ncs_dihedral_manager=None,
        den_manager=None,
        chirality_proxies=None,
        planarity_proxies=None,
        parallelity_proxies=None,
        ramachandran_manager=None,
        external_energy_function=None,
        plain_pairs_radius=None,
        max_reasonable_bond_distance=None,
        min_cubicle_edge=5,
        log=StringIO()):
    super(manager, self).__init__()
    if (site_symmetry_table is not None): assert crystal_symmetry is not None
    if (bond_params_table is not None and site_symmetry_table is not None):
      assert bond_params_table.size() == site_symmetry_table.indices().size()
    if (shell_sym_tables is not None and site_symmetry_table is not None):
      assert len(shell_sym_tables) > 0
      assert shell_sym_tables[0].size() == site_symmetry_table.indices().size()
    if (nonbonded_types is not None and site_symmetry_table is not None):
      assert nonbonded_types.size() == site_symmetry_table.indices().size(), "%d != %d" % (
          nonbonded_types.size(), site_symmetry_table.indices().size())
    if (nonbonded_types is not None) and (nonbonded_charges is not None):
      assert nonbonded_charges.size() == nonbonded_types.size(), "%d != %d" % (
          nonbonded_charges.size(), nonbonded_types.size())
    adopt_init_args(self, locals(), exclude=["log"])
    self.reset_internals()

# implement explicit pickling for the log object since StringIO doesn't support pickling
  # def __getstate__(self):
  #   state = self.__dict__.copy()
  #   if hasattr(self.log, "getvalue"):
  #     state[ "log" ] = self.log.getvalue()
  #   else:
  #     state[ "log" ] = ""
  #   return state

  # def __setstate__(self, state):
  #   state[ "log" ] = StringIO( state[ "log" ] )
  #   self.__dict__.update( state )

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

  def replace_site_symmetry(self, new_site_symmetry_table,
                                  special_position_settings,
                                  sites_cart):
    assert self.site_symmetry_table is not None
    self.site_symmetry_table = new_site_symmetry_table
    # Recompute shell_sym_tables
    asu_mappings = special_position_settings.asu_mappings(
      buffer_thickness = self.max_reasonable_bond_distance*3 if \
      self.max_reasonable_bond_distance is not None \
      else 150) # arbitrary big number if not defined
    asu_mappings.process_sites_cart(
      original_sites      = sites_cart,
      site_symmetry_table = self.site_symmetry_table)
    bond_asu_table = crystal.pair_asu_table(asu_mappings = asu_mappings)
    # Add all previously defined bonds
    bond_asu_table.add_pair_sym_table(self.shell_sym_tables[0])
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table = bond_asu_table,
      max_shell      = 3)
    self.shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]
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
      for proxy in self.get_dihedral_proxies():
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
      for proxy in self.get_dihedral_proxies():
        if (ec(edge=proxy.i_seqs[1:3]) in ["hinge", "loop"]):
          result.append(proxy)
      if (result.size() == 0):
        return None
      return result
    reduced_dihedral_proxies = get()
    #
    if not include_den_restraints:
      result = manager(
        crystal_symmetry=self.crystal_symmetry,
        site_symmetry_table=self.site_symmetry_table,
        bond_params_table=self.bond_params_table,
        shell_sym_tables=reduced_shell_sym_tables,
        angle_proxies=reduced_angle_proxies,
        dihedral_proxies=reduced_dihedral_proxies,
        ncs_dihedral_manager=self.ncs_dihedral_manager)
    else:
      result = manager(
        crystal_symmetry=self.crystal_symmetry,
        site_symmetry_table=self.site_symmetry_table,
        bond_params_table=self.bond_params_table,
        shell_sym_tables=reduced_shell_sym_tables,
        angle_proxies=reduced_angle_proxies,
        dihedral_proxies=reduced_dihedral_proxies,
        ncs_dihedral_manager=self.ncs_dihedral_manager,
        den_manager=self.den_manager,
        ramachandran_manager=self.ramachandran_manager)
    result.set_source(source = self.get_source())
    return result

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
      for i_seq in range(n_additional_sites):
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
    if (self.nonbonded_charges is not None):
      if (nonbonded_charges is None):
        nonbonded_charges = flex.int(n_additional_sites, 0)
      assert (nonbonded_charges.size() == n_additional_sites)
      nonbonded_charges = self.nonbonded_charges.concatenate(
        nonbonded_charges)
    result = manager(
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
      reference_coordinate_proxies=self.reference_coordinate_proxies,
      reference_dihedral_manager=self.reference_dihedral_manager,
      ramachandran_manager=self.ramachandran_manager,
      ncs_dihedral_manager=self.ncs_dihedral_manager,
      den_manager=self.den_manager,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      parallelity_proxies=self.parallelity_proxies,
      plain_pairs_radius=self.plain_pairs_radius,
      max_reasonable_bond_distance=self.max_reasonable_bond_distance)
    result.set_source(source = self.get_source())
    return result

  def select(self, selection=None, iselection=None):
    assert [selection, iselection].count(None) == 1

    n_seqs = dict_with_default_0()
    if (iselection is None):
      iselection = selection.iselection()
      n_seqs[selection.size()] += 1

    selected_stuff = [None]*7
    for i, stuff in enumerate([self.model_indices, self.conformer_indices,
        self.sym_excl_indices, self.donor_acceptor_excl_groups,
        self.nonbonded_types,
        self.nonbonded_charges]):
      if stuff is not None:
        selected_stuff[i] = stuff.select(iselection)
        n_seqs[stuff.size()] += 1
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
      selected_shell_sym_tables = [shell_sym_table.proxy_select(
          iselection)
        for shell_sym_table in self.shell_sym_tables]
      if (len(self.shell_sym_tables) > 0):
        n_seqs[self.shell_sym_tables[0].size()] += 1
    selected_nonbonded_types = None

    def get_n_seq():
      if (len(n_seqs) == 0):
        raise RuntimeError("Cannot determine n_seq.")
      if (len(n_seqs) != 1):
        raise RuntimeError("Selection size mismatches: %s." % str(n_seqs))
      return list(n_seqs.keys())[0]  # FIXME this might break py2/3 compat unless there is only 1 key

    n_seq = get_n_seq()

    selected_proxies = [None]*11
    for i, proxies in enumerate([self.angle_proxies, self.dihedral_proxies,
        self.reference_coordinate_proxies, self.reference_dihedral_manager,
        self.ncs_dihedral_manager, self.den_manager, self.chirality_proxies,
        self.planarity_proxies, self.parallelity_proxies,
        self.ramachandran_manager]):
      if proxies is not None:
        selected_proxies[i] = proxies.proxy_select(n_seq, iselection)

    result = manager(
      crystal_symmetry=self.crystal_symmetry,
      site_symmetry_table=selected_site_symmetry_table,
      model_indices=selected_stuff[0],
      conformer_indices=selected_stuff[1],
      sym_excl_indices=selected_stuff[2],
      donor_acceptor_excl_groups=selected_stuff[3],
      bond_params_table=selected_bond_params_table,
      shell_sym_tables=selected_shell_sym_tables,
      nonbonded_params=self.nonbonded_params,
      nonbonded_types=selected_stuff[4],
      nonbonded_charges=selected_stuff[5],
      nonbonded_function=self.nonbonded_function,
      nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
      nonbonded_buffer=self.nonbonded_buffer,
      angle_proxies=selected_proxies[0],
      dihedral_proxies=selected_proxies[1],
      reference_coordinate_proxies=selected_proxies[2],
      reference_dihedral_manager=selected_proxies[3],
      ramachandran_manager=selected_proxies[9],
      ncs_dihedral_manager=selected_proxies[4],
      den_manager=selected_proxies[5],
      chirality_proxies=selected_proxies[6],
      planarity_proxies=selected_proxies[7],
      parallelity_proxies=selected_proxies[8],
      plain_pairs_radius=self.plain_pairs_radius,
      max_reasonable_bond_distance=self.max_reasonable_bond_distance)
    result.set_source(source = self.get_source())
    return result

  def shift_sites_cart(self, shift):
    # just a stub from RSR. Proper treatment can be implemented here.
    self.remove_reference_coordinate_restraints_in_place()

  def discard_symmetry(self, new_unit_cell):
    assert self.site_symmetry_table is not None #XXX lazy
    assert self.shell_sym_tables is not None #XXX lazy
    result = manager(
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
      reference_coordinate_proxies=self.reference_coordinate_proxies,
      reference_dihedral_manager=self.reference_dihedral_manager,
      ramachandran_manager=self.ramachandran_manager,
      ncs_dihedral_manager=self.ncs_dihedral_manager,
      den_manager=self.den_manager,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      parallelity_proxies=self.parallelity_proxies,
      plain_pairs_radius=self.plain_pairs_radius,
      max_reasonable_bond_distance=self.max_reasonable_bond_distance)
    result.set_source(source = self.get_source())
    return result

  def add_angles_in_place(self, additional_angle_proxies):
    self.angle_proxies.extend(additional_angle_proxies)

  def remove_angles_in_place(self, selection):
    self.angle_proxies = self.angle_proxies.proxy_remove(
      selection=selection)

  #=================================================================
  # User-supplied restraints (geometry_restraints.edits scope)
  #=================================================================
  def get_user_supplied_restraints(self):
    bonds_simpe = None
    bonds_asu = None
    angles = None
    planarity = None
    parallelity = None
    pair_proxies = self.pair_proxies()
    specific_origin_id = origin_ids.get_origin_id('edits')
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        bonds_simpe = pair_proxies.bond_proxies.simple.proxy_select(
          origin_id=specific_origin_id)
        bonds_asu = pair_proxies.bond_proxies.asu.proxy_select(
          origin_id=specific_origin_id)
    angles = self.angle_proxies.proxy_select(
      origin_id=specific_origin_id)
    dihedrals = self.dihedral_proxies.proxy_select(
      origin_id=specific_origin_id)
    planarity = self.planarity_proxies.proxy_select(
      origin_id=specific_origin_id)
    parallelity = self.parallelity_proxies.proxy_select(
      origin_id=specific_origin_id)
    return bonds_simpe, bonds_asu, angles, dihedrals, planarity, parallelity

  def remove_user_supplied_restraints_in_place(self):
    """
    Problem is in removing bond restraints, which is not yet implemented at all
    and rather complicated.
    """
    raise NotImplementedError

  def get_bond_proxies_without_user_supplied(self, sites_cart=None):
    pair_proxies = self.pair_proxies(sites_cart=sites_cart)
    simple = None
    asu = None
    specific_origin_id = origin_ids.get_origin_id('edits')
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        if pair_proxies.bond_proxies.simple is not None:
          simple = pair_proxies.bond_proxies.simple.get_proxies_without_origin_id(
            origin_id=specific_origin_id)
        if pair_proxies.bond_proxies.asu is not None:
          asu = pair_proxies.bond_proxies.asu.get_proxies_without_origin_id(
            origin_id=specific_origin_id)
    return simple, asu

  def get_angle_proxies_without_user_supplied(self):
    specific_origin_id = origin_ids.get_origin_id('edits')
    if self.angle_proxies is not None:
      return self.angle_proxies.proxy_remove(origin_id=specific_origin_id)
    return None

  def get_dihedral_proxies_without_user_supplied(self):
    specific_origin_id = origin_ids.get_origin_id('edits')
    if self.dihedral_proxies is not None:
      return self.dihedral_proxies.proxy_remove(origin_id=specific_origin_id)
    return None

  def get_planarity_proxies_without_user_supplied(self):
    specific_origin_id = origin_ids.get_origin_id('edits')
    if self.planarity_proxies is not None:
      return self.planarity_proxies.proxy_remove(origin_id=specific_origin_id)
    return None

  def get_parallelity_proxies_without_user_supplied(self):
    specific_origin_id = origin_ids.get_origin_id('edits')
    if self.parallelity_proxies is not None:
      return self.parallelity_proxies.proxy_remove(origin_id=specific_origin_id)
    return None

  #=================================================================
  # Reference coordinate proxies methods
  #=================================================================
  def get_reference_coordinate_proxies(self):
    return self.reference_coordinate_proxies

  def adopt_reference_coordinate_restraints_in_place(self,
      reference_coordinate_proxies):
    self.reference_coordinate_proxies = reference_coordinate_proxies

  def remove_reference_coordinate_restraints_in_place(self,
      selection=None):
    if (selection is not None):
      self.reference_coordinate_proxies = \
        self.reference_coordinate_proxies.proxy_remove(selection=selection)
    else :
      self.reference_coordinate_proxies = None

  def get_n_reference_coordinate_proxies(self):
    if self.reference_coordinate_proxies is not None:
      return self.reference_coordinate_proxies.size()
    else:
      return 0

  def append_reference_coordinate_restraints_in_place(self,
      reference_coordinate_proxies):
    if self.reference_coordinate_proxies is not None:
      self.reference_coordinate_proxies.extend(
          reference_coordinate_proxies)
    else:
      self.reference_coordinate_proxies = reference_coordinate_proxies

  def add_reference_coordinate_restraints_in_place(self,
      all_chain_proxies=None,
      pdb_hierarchy=None,
      selection=None,
      exclude_outliers=True,
      sigma=0.2,
      limit=1.0,
      top_out=False,
      n_atoms_in_target_model=None,
      reference_is_average_alt_confs=None,
      ):
    assert [all_chain_proxies, pdb_hierarchy].count(None) == 1
    assert [all_chain_proxies, n_atoms_in_target_model].count(None) == 1
    if all_chain_proxies is None:
      assert isinstance(selection, flex.size_t)
    from mmtbx.geometry_restraints.reference import add_coordinate_restraints, \
        exclude_outliers_from_reference_restraints_selection

    if isinstance(selection, flex.size_t):
      isel = selection
      bsel = flex.bool(n_atoms_in_target_model, isel)
      sites_cart=pdb_hierarchy.atoms().extract_xyz()
      assert not reference_is_average_alt_confs
    else:
      # should be deleted if all_chain_proxies won't be used
      if reference_is_average_alt_confs:
        hierarchy = all_chain_proxies.pdb_hierarchy.as_forward_compatible_hierarchy()
        asel = hierarchy.average_alt_confs(pinch_limit=limit)
        sites_cart = hierarchy.atoms().extract_xyz()
      else:
        sites_cart = all_chain_proxies.pdb_hierarchy.atoms().extract_xyz()
      new_selection = flex.bool(sites_cart.size(), True)
      if selection is not None:
        new_selection = all_chain_proxies.selection(selection)
        if (new_selection.size() == 0):
          raise Sorry(("No atoms selected for harmonic restraints (input "+
            "selection string: %s)") % selection)
      # print >> self.log, "*** Restraining %d atoms to initial coordinates ***" % \
      #     new_selection.size()
      isel = new_selection.iselection()
    if exclude_outliers:
      new_selection = exclude_outliers_from_reference_restraints_selection(
          pdb_hierarchy=pdb_hierarchy if pdb_hierarchy is not None else all_chain_proxies.pdb_hierarchy,
          restraints_selection=new_selection if all_chain_proxies is not None else bsel)
      isel = new_selection.iselection()
    if reference_is_average_alt_confs:
      isel = asel
    proxies = add_coordinate_restraints(
        sites_cart=sites_cart.select(isel) if n_atoms_in_target_model is None else sites_cart,
        selection=isel,
        sigma=sigma,
        limit=limit,
        top_out_potential=top_out)
    if proxies.size() > 0:
      self.reference_coordinate_proxies = proxies

  #=================================================================
  # Torsion (dihedral) restraints on chi angles (side-chain rotamers)
  # proxies methods
  #=================================================================
  def add_chi_torsion_restraints_in_place(self,
      pdb_hierarchy,
      sites_cart,
      selection=None,
      sigma=2.5,
      limit=15.0,
      chi_angles_only=False,
      top_out_potential=False):
    from mmtbx.geometry_restraints.reference import generate_torsion_restraints
    specific_origin_id = origin_ids.get_origin_id('chi angles')
    chi_torsions = generate_torsion_restraints(
        pdb_hierarchy=pdb_hierarchy,
        sites_cart=sites_cart,
        selection=selection,
        sigma=sigma,
        limit=limit,
        chi_angles_only=chi_angles_only,
        top_out_potential=top_out_potential,
        origin_id=specific_origin_id)
    self.add_dihedrals_in_place(chi_torsions)

  def remove_chi_torsion_restraints_in_place(self, selection=None):
    specific_origin_id = origin_ids.get_origin_id('chi angles')
    if self.dihedral_proxies is not None:
      if selection is None:
        self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
          origin_id=specific_origin_id)
      else:
        chi_proxies = self.dihedral_proxies.proxy_select(
          origin_id=specific_origin_id)
        should_remain_proxies = chi_proxies.proxy_remove(selection=selection)
        self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
          origin_id=specific_origin_id)
        self.add_dihedrals_in_place(should_remain_proxies)

  def get_chi_torsion_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('chi angles')
    if self.dihedral_proxies is not None:
      return self.dihedral_proxies.proxy_select(origin_id=specific_origin_id)

  def get_n_chi_torsion_proixes(self):
    specific_origin_id = origin_ids.get_origin_id('chi angles')
    if self.dihedral_proxies is not None:
      return self.dihedral_proxies.proxy_select(origin_id=specific_origin_id).size()


  #=================================================================
  # Dihedral proxies methods
  #=================================================================
  def add_dihedrals_in_place(self,
      additional_dihedral_proxies,
      check_for_duplicates=True):
    if self.dihedral_proxies is not None:
      if check_for_duplicates:
        # first remove duplicates
        n_atoms = 0
        selection = None
        sel = []
        for p in additional_dihedral_proxies:
          for i_seq in p.i_seqs:
            sel.append(i_seq)
        selection = flex.size_t(sorted(list(set(sel))))
        self.remove_dihedrals_in_place(selection)
      self.dihedral_proxies.extend(additional_dihedral_proxies)
    else:
      self.dihedral_proxies = additional_dihedral_proxies

  def remove_dihedrals_in_place(self, selection):
    if self.dihedral_proxies is not None:
      self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
        selection=selection)

  def get_dihedral_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    if self.dihedral_proxies is not None:
      return self.dihedral_proxies.proxy_select(origin_id=specific_origin_id)
    return None

  #=================================================================
  # C-beta dihedral proxies methods
  #=================================================================
  def get_c_beta_torsion_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('C-beta')
    if self.dihedral_proxies is not None:
      return self.dihedral_proxies.proxy_select(origin_id=specific_origin_id)
    return None

  def get_n_c_beta_torsion_proxies(self):
    tp = self.get_c_beta_torsion_proxies()
    if tp is not None:
      return tp.size()
    return 0

  def remove_c_beta_torsion_restraints_in_place(self, selection=None):
    specific_origin_id = origin_ids.get_origin_id('C-beta')
    if selection is None:
      self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
        origin_id=specific_origin_id)
    else:
      self.dihedral_proxies = self.dihedral_proxies.proxy_select(
        origin_id=specific_origin_id).proxy_remove(selection=selection)

  #=================================================================
  # Reference dihedral proxies methods
  #=================================================================
  def adopt_reference_dihedral_manager(self, manager):
    self.reference_dihedral_manager = manager

  def remove_reference_dihedral_manager(self):
    self.reference_dihedral_manager = None

  def remove_reference_dihedrals_in_place(self, selection):
    if self.reference_dihedral_manager is not None:
      self.reference_dihedral_manager.proxy_remove(selection=selection)

  def remove_ncs_dihedrals_in_place(self):
    if self.ncs_dihedral_manager is not None:
      self.ncs_dihedral_manager.remove_ncs_dihedrals_in_place()

  def update_dihedral_ncs_restraints(self, model, log):
    if self.ncs_dihedral_manager is not None:
      self.ncs_dihedral_manager.update_dihedral_ncs_restraints(
          model=model,
          log=log)

  def get_n_reference_dihedral_proxies(self):
    if self.reference_dihedral_manager is not None:
      return self.reference_dihedral_manager.get_n_proxies()

  def sync_reference_dihedral_with_ncs(self, log):
    if (self.reference_dihedral_manager is not None and
        self.ncs_dihedral_manager is not None):
      self.reference_dihedral_manager.remove_restraints_with_ncs_matches(
          ncs_dihedral_proxies = self.ncs_dihedral_manager.ncs_dihedral_proxies,
          ncs_match_hash       = self.ncs_dihedral_manager.ncs_match_hash,
          log=log)

  #=================================================================
  # DEN manager/proxies methods
  #=================================================================
  def create_den_manager(self, den_params, pdb_hierarchy, log):
    from mmtbx.den import den_restraints
    if den_params.reference_file is not None:
      pdb_io_ref = iotbx.pdb.input(den_params.reference_file)
      pdb_hierarchy_ref = pdb_io_ref.construct_hierarchy()
      pdb_hierarchy_ref.atoms().reset_i_seq()
    else: #restrain model to starting coordinates
      pdb_hierarchy_ref = None
    den_manager = den_restraints(
      pdb_hierarchy=pdb_hierarchy,
      pdb_hierarchy_ref=pdb_hierarchy_ref,
      params=den_params,
      log=log)
    self.adopt_den_manager(den_manager)

  def adopt_den_manager(self, den_manager):
    self.den_manager = den_manager

  def get_n_den_proxies(self):
    if self.den_manager is not None:
      return self.den_manager.get_n_proixes()
    return 0

  def remove_chiralities_in_place(self, selection):
    self.chirality_proxies = self.chirality_proxies.proxy_remove(
      selection=selection)

  def add_planarities_in_place(self, additional_planarity_proxies):
    self.planarity_proxies.extend(additional_planarity_proxies)

  def remove_planarities_in_place(self, selection):
    self.planarity_proxies = self.planarity_proxies.proxy_remove(
      selection=selection)

  def add_parallelities_in_place(self, additional_parallelity_proxies):
    self.parallelity_proxies.extend(additional_parallelity_proxies)

  def remove_parallelities_in_place(self, selection):
    self.parallelity_proxies = self.parallelity_proxies.proxy_remove(
      selection=selection)

  #=================================================================
  # Ramachandran manager/proxies methods
  #=================================================================
  def set_ramachandran_restraints(self, manager):
    self.ramachandran_manager=manager

  def update_ramachandran_restraints_phi_psi_targets(self, hierarchy):
    if(self.ramachandran_manager is not None):
      self.ramachandran_manager.update_phi_psi_targets(hierarchy = hierarchy)

  def remove_ramachandran_in_place(self):
    self.ramachandran_manager = None

  def get_n_ramachandran_proxies(self):
    if self.ramachandran_manager is not None:
      return self.ramachandran_manager.get_n_proxies()
    return 0

  #=================================================================
  # Secondary structure manager/proxies methods
  #=================================================================
  def set_secondary_structure_restraints(self, ss_manager, hierarchy, log):
    (hb_proxies, hb_angle_proxies, planarity_proxies,
    parallelity_proxies) = ss_manager.create_all_new_restraints(
        pdb_hierarchy=hierarchy,
        grm=self,
        log=log)
    # print ('Total proxies:', len(hb_proxies))
    # for p in hb_proxies:
    #   print (p.i_seqs[0], p.i_seqs[1], hierarchy.atoms()[p.i_seqs[0]].id_str(), '<-->', hierarchy.atoms()[p.i_seqs[1]].id_str())
    #   print ("  ", hierarchy.atoms()[p.i_seqs[0]].distance(hierarchy.atoms()[p.i_seqs[1]]))
    self.add_new_hbond_restraints_in_place(
        proxies=hb_proxies,
        sites_cart=hierarchy.atoms().extract_xyz(),
        max_distance_between_connecting_atoms=10)
    self.add_angles_in_place(hb_angle_proxies)
    self.add_planarities_in_place(planarity_proxies)
    self.add_parallelities_in_place(parallelity_proxies)

  def remove_secondary_structure_restraints(self, sites_cart):
    bp_origin_id = origin_ids.get_origin_id('hydrogen bonds')
    bonds_to_remove = []
    for proxy_type in self.get_all_bond_proxies():
      for proxy in proxy_type:
        if proxy.origin_id == bp_origin_id:
          bonds_to_remove.append(proxy.i_seqs)
    self.remove_bond_restraints_in_place(
        bonded_pairs=bonds_to_remove,
        sites_cart=sites_cart)
    if self.angle_proxies is not None:
      self.angle_proxies = self.angle_proxies.proxy_remove(
          origin_id=origin_ids.get_origin_id('hydrogen bonds'))
    if self.planarity_proxies is not None:
      self.planarity_proxies = self.planarity_proxies.proxy_remove(
          origin_id=origin_ids.get_origin_id('basepair planarity'))
    if self.parallelity_proxies is not None:
      self.parallelity_proxies = self.parallelity_proxies.proxy_remove(
          origin_id=origin_ids.get_origin_id('basepair parallelity'))
    if self.parallelity_proxies is not None:
      self.parallelity_proxies = self.parallelity_proxies.proxy_remove(
          origin_id=origin_ids.get_origin_id('basepair stacking'))
    return None

  def add_parallelity_proxies_for_side_chain(self, hierarchy, log):
    parallelity_proxies = create_side_chain_restraints( hierarchy=hierarchy,
                                                        grm=self,
                                                        log=log)
    self.add_parallelities_in_place(parallelity_proxies)

  def set_external_energy_function(self, energy_function):
    self.external_energy_function = energy_function

  def _get_n_bond_proxies_origin(self, origin_id):
    pair_proxies = self.pair_proxies()
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        return len(pair_proxies.bond_proxies.simple.proxy_select(origin_id=origin_id))+\
            len(pair_proxies.bond_proxies.asu.proxy_select(origin_id=origin_id))
    return 0

  def get_n_bond_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    return self._get_n_bond_proxies_origin(origin_id=specific_origin_id)

  def get_covalent_bond_proxies(self, sites_cart=None):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    pair_proxies = self.pair_proxies(sites_cart=sites_cart)
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        return (pair_proxies.bond_proxies.simple.proxy_select(
          origin_id=specific_origin_id),
        pair_proxies.bond_proxies.asu.proxy_select(
          origin_id=specific_origin_id))

  def get_all_bond_proxies(self, sites_cart=None):
    pair_proxies = self.pair_proxies(sites_cart=sites_cart)
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        return (pair_proxies.bond_proxies.simple,
                pair_proxies.bond_proxies.asu)

  def get_covalent_angle_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    return self.angle_proxies.proxy_select(origin_id=specific_origin_id)

  def get_all_angle_proxies(self):
    return self.angle_proxies

  def get_n_hbond_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('hydrogen bonds')
    return self._get_n_bond_proxies_origin(origin_id=specific_origin_id)

  def get_n_angle_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    return len(self.angle_proxies.proxy_select(origin_id=specific_origin_id))

  def get_n_hangle_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('hydrogen bonds')
    return len(self.angle_proxies.proxy_select(origin_id=specific_origin_id))

  def get_n_stacking_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('basepair stacking')
    return len(self.parallelity_proxies.proxy_select(
      origin_id=specific_origin_id))

  def get_n_parallelity_bp_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('basepair parallelity')
    return len(self.parallelity_proxies.proxy_select(
      origin_id=specific_origin_id))

  def get_n_planarity_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('covalent geometry')
    return len(self.planarity_proxies.proxy_select(
      origin_id=specific_origin_id))

  def get_n_planarity_bp_proxies(self):
    specific_origin_id = origin_ids.get_origin_id('basepair planarity')
    return len(self.planarity_proxies.proxy_select(
      origin_id=specific_origin_id))

  def get_hbond_proxies_iseqs(self):
    specific_origin_id = origin_ids.get_origin_id('hydrogen bonds')
    result = []
    try:
      pair_proxies = self.pair_proxies()
    except AssertionError as e:
      # This is in case when somebody tries to access hbonds when
      # grm is not ready. See e.g. phenix/refinement/runtime.py
      # class refinement_status, grm.get_hbond_proxies_iseqs()
      if "pair_proxies not defined already." in e.args:
        return result
      else:
        raise
    if pair_proxies is not None:
      if pair_proxies.bond_proxies is not None:
        simple_p = pair_proxies.bond_proxies.simple.proxy_select(
          origin_id=specific_origin_id)
        asu_p = pair_proxies.bond_proxies.asu.proxy_select(
          origin_id=specific_origin_id)
        for p in simple_p:
          result.append((p.i_seqs[0], p.i_seqs[1]))
        for p in asu_p:
          result.append((p.i_seq, p.j_seq))
    return result

  def new_included_bonded_atoms(self, proxies, sites_cart,
      site_symmetry_table, nonbonded_types, nonbonded_charges,
      max_distance_between_connecting_atoms=5,
      skip_max_proxy_distance_calculation=False):

    """ Produce new geometry_restraints_manager object that will
    include new atoms each with exactly one bond to the existing atom.
    proxies - list of bond_proxy objects. Symmetry operation will be determined
    automatically, so proxy.rt_mx_ji is ignored.
    Essentially this function wraps self.new_including_isolated_sites() and
    self.add_new_bond_restraints_in_place().
    sites_cart should contain all coordinates, for old and new atoms. new
      coordinates should follow old ones and the order should be consistent
      with i_seqs mentioned in proxy objects.
    site_symmetry_table - table only for new atoms
    nonbonded_types - only for new atoms
    nonbonded_charges - only for new atoms."""

    number_of_new_atoms = len(proxies)

    assert site_symmetry_table.indices().size() == number_of_new_atoms
    assert len(nonbonded_types) == number_of_new_atoms
    assert len(nonbonded_charges) == number_of_new_atoms

    if (self.model_indices is None):
      model_indices = None
    else:
      model_indices = flex.size_t(number_of_new_atoms, 0)
    if (self.conformer_indices is None):
      conformer_indices = None
    else:
      conformer_indices = flex.size_t(number_of_new_atoms, 0)
    if (self.sym_excl_indices is None):
      sym_excl_indices = None
    else:
      sym_excl_indices = flex.size_t(number_of_new_atoms, 0)
    if (self.donor_acceptor_excl_groups is None):
      donor_acceptor_excl_groups = None
    else:
      donor_acceptor_excl_groups = flex.size_t(number_of_new_atoms, 0)
    new_grm = self.new_including_isolated_sites(
        n_additional_sites =number_of_new_atoms,
        model_indices=model_indices,
        conformer_indices=conformer_indices,
        sym_excl_indices=sym_excl_indices,
        donor_acceptor_excl_groups=donor_acceptor_excl_groups,
        site_symmetry_table=site_symmetry_table,
        nonbonded_types=nonbonded_types,
        nonbonded_charges=nonbonded_charges)
    sites_frac = self.crystal_symmetry.unit_cell().\
        fractionalize(sites_cart=sites_cart)
    new_grm.update_plain_pair_sym_table(sites_frac)
    new_grm.add_new_bond_restraints_in_place(proxies, sites_cart,
      max_distance_between_connecting_atoms=max_distance_between_connecting_atoms,
      skip_max_proxy_distance_calculation=skip_max_proxy_distance_calculation)
    new_grm.set_source(source = self.get_source())
    return new_grm

  def add_new_hbond_restraints_in_place(self, proxies, sites_cart,
      max_distance_between_connecting_atoms=5,
      skip_max_proxy_distance_calculation=False):
    self.add_new_bond_restraints_in_place(proxies, sites_cart,
        max_distance_between_connecting_atoms,
        skip_max_proxy_distance_calculation)

  def add_new_bond_restraints_in_place(self, proxies, sites_cart,
      max_distance_between_connecting_atoms=5,
      skip_max_proxy_distance_calculation=False):
    """ Add new bond restraints for list of proxies to this
    geometry restraints manager, _in_place_! Returns nothing.
    proxies - list of bond_proxy objects. The symmetry operation for the
    paired atoms is determined here, therefore the proxy.rt_mx_ji may be
    anything."""
    import time
    if len(proxies) == 0:
      return
    rt_mx_ji_options = [[] for x in proxies]
    # Get current max bond distance, copied from pair_proxies()
    t0 = time.time()
    bonded_distance_cutoff = max_distance_between_connecting_atoms
    sites_frac = self.crystal_symmetry.unit_cell().\
        fractionalize(sites_cart=sites_cart)
    existing_max_bonded_distance = 0
    for shell_sym_table in self.shell_sym_tables:
      existing_max_bonded_distance = flex.max_default(
              values=crystal.get_distances(
                  pair_sym_table=shell_sym_table,
                  orthogonalization_matrix=self.crystal_symmetry.unit_cell() \
                      .orthogonalization_matrix(),
                  sites_frac=sites_frac),
              default=0)
    t1 = time.time()
    max_p_distance=0
    if not skip_max_proxy_distance_calculation:
      for p in proxies:
        distance_model = geometry_restraints.bond(
            [sites_cart[p.i_seqs[0]],sites_cart[p.i_seqs[1]]],
            distance_ideal=0,
            weight=1).distance_model
        if distance_model > max_p_distance:
          max_p_distance = distance_model
      bonded_distance_cutoff = max(bonded_distance_cutoff,
          max_p_distance)
    bonded_distance_cutoff = min(bonded_distance_cutoff,
        max_distance_between_connecting_atoms)+0.1
    bonded_distance_cutoff = max(bonded_distance_cutoff,
        existing_max_bonded_distance)
    t2 = time.time()
    # print("bonded_distance_cutoff", bonded_distance_cutoff)
    # make asu mappings
    all_asu_mappings = self.crystal_symmetry.special_position_settings().\
        asu_mappings(buffer_thickness=bonded_distance_cutoff)
    all_asu_mappings.process_sites_cart(
      original_sites=sites_cart,
      site_symmetry_table=self.site_symmetry_table)
    # Add all previously defined bonds
    all_bonds_asu_table = crystal.pair_asu_table(asu_mappings=all_asu_mappings)
    all_bonds_asu_table.add_pair_sym_table(self.shell_sym_tables[0])

    t3 = time.time()
    proxies_i_seqs = {}
    proxies_iselection = []
    for np, p in enumerate(proxies):
      # print('p iseqs', p.i_seqs)
      proxies_i_seqs[p.i_seqs] = np
      proxies_i_seqs[(p.i_seqs[1], p.i_seqs[0])] = np
      for i in list(p.i_seqs):
        if i not in proxies_iselection:
          proxies_iselection.append(i)
    # Not sure whether we want to sort it...
    # proxies_iselection = flex.size_t(sorted(proxies_iselection))
    # print('proxies_i_seqs', proxies_i_seqs)
    proxies_iselection = flex.size_t(proxies_iselection)
    t4 = time.time()

    # Generate pairs for connecting atoms only - should be much faster then
    # doing the same for all atoms in geometry_restraints_manager.
    # Note, that this will generate not only useful pairs but all pairs
    # within bonded_distance_cutoff cutoff. Therefore
    # later we will filter them.
    conn_asu_mappings = self.crystal_symmetry.special_position_settings().\
      asu_mappings(buffer_thickness=bonded_distance_cutoff)
    connecting_sites_cart = sites_cart.select(proxies_iselection)
    conn_site_symmetry_table = self.site_symmetry_table.select(
        proxies_iselection)
    conn_asu_mappings.process_sites_cart(
        original_sites=connecting_sites_cart,
        site_symmetry_table=conn_site_symmetry_table)
    conn_pair_asu_table = crystal.pair_asu_table(
        asu_mappings=conn_asu_mappings)
    conn_pair_asu_table.add_all_pairs(
        distance_cutoff=bonded_distance_cutoff)
    pair_generator = crystal.neighbors_fast_pair_generator(
        conn_asu_mappings,
        distance_cutoff=bonded_distance_cutoff)

    t5 = time.time()
    # r_a connects i_seqs in pair_generator with original i_seqs
    r_a = list(reindexing_array(len(sites_cart), proxies_iselection.as_int()))
    r_a_index = [0]*(len(sites_cart)+1)
    for index, value in enumerate(r_a):
      r_a_index[value] = index
    n_added_proxies = 0
    for pair in pair_generator:
      # !!! Warning !!! pairs here may come reverted: (32,24) and (24,32) and
      # they may result in different rt_mx_ji. Therefore in the next for cycle
      # we need to make sure we are picking correct rt_mx_ji
      pair_in_origninal_indeces = (r_a_index[pair.i_seq], r_a_index[pair.j_seq])
      n_proxy = None
      # Is this pair should be restrained?
      try:
        n_proxy = proxies_i_seqs.get(pair_in_origninal_indeces, None)
      except ValueError:
        # there is no proxy for this pair, so we will not make a bond for it
        continue
      if n_proxy is None:
        continue
      # Sanity check (not necessary because of 'continue' in previous line)
      # print "  n_proxy", n_proxy
      if n_proxy is not None:
        # print("Adding to options,", n_proxy, rt_mx_ji_options)
        #Trying to find rt_mx_ji for connecting atoms
        rt_mx_i = conn_asu_mappings.get_rt_mx_i(pair)
        rt_mx_j = conn_asu_mappings.get_rt_mx_j(pair)
        rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
        rt_mx_ji_options[n_proxy].append(rt_mx_ji)
        # print("  pair:",  pair.i_seq, pair.j_seq, "n_proxy ", n_proxy,rt_mx_ji)
    # print ("rt_mx_ji_options:", rt_mx_ji_options)
    # for i, opts in enumerate(rt_mx_ji_options):
    #   print("  ",i,' ', end='')
    #   for o in opts:
    #     print(o,",", end='')
    #   print()
    # print()
    # print("STARTING TO PICK rt_mx_ji for proxies")
    for proxy, rt_mx_ji_option in zip(proxies, rt_mx_ji_options):
      # print("rt_mx_ji_options:",end='')
      # for op in rt_mx_ji_option:
      #   print (op,' | ' , end='')
      # print()
      # print("  choose rt_mx_ji")
      rt_mx_ji = None
      if len(rt_mx_ji_option) == 1:
        rt_mx_ji = rt_mx_ji_option[0]
      if len(rt_mx_ji_option) > 1:
        # search for unit mx or pick the first one
        for rmj in rt_mx_ji_option:
          if rmj.is_unit_mx():
            rt_mx_ji = rmj
        if rt_mx_ji is None:
          # no unit mx found and we have >1 options, need to pick correct one
          for rmj in rt_mx_ji_option:
            t_rt_mx_i = all_asu_mappings.get_rt_mx(proxy.i_seqs[0], 0)
            t_rt_mx_j = t_rt_mx_i.multiply(rmj)
            j_sym = all_asu_mappings.find_i_sym(proxy.i_seqs[1], t_rt_mx_j)
            # print('checking j_sym', j_sym, proxy.i_seqs[0], proxy.i_seqs[1])
            if j_sym >= 0:
              rt_mx_ji = rmj
      # Add new defined bond
      # print ("  chosen rt_mx_ji", rt_mx_ji)
      if rt_mx_ji is not None:
        # print("  Adding new bond:", proxy.i_seqs[0], proxy.i_seqs[1], rt_mx_ji)
        all_bonds_asu_table.add_pair(
          i_seq=proxy.i_seqs[0],
          j_seq=proxy.i_seqs[1],
          rt_mx_ji=rt_mx_ji)
        # Update with new bond
        self.bond_params_table.update(
            i_seq=proxy.i_seqs[0],
            j_seq=proxy.i_seqs[1],
            params=geometry_restraints.bond_params(
                distance_ideal = proxy.distance_ideal,
                weight         = proxy.weight,
                slack          = proxy.slack,
                limit          = proxy.limit,
                top_out        = proxy.top_out,
                origin_id      = proxy.origin_id)
            )
        n_added_proxies += 1
    t6 = time.time()
    # print('n_added_proxies', n_added_proxies)
    # update self.shell_sym_tables with new bonds
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table=all_bonds_asu_table,
      max_shell=3)
    self.shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]
    self.reset_internals()
    # Run this function so that new pair_proxies are ready to go
    t61 = time.time()
    self.pair_proxies(sites_cart=sites_cart)
    t7 = time.time()
    # print "times in add_new_bond_restraints_in_place:"
    # print "t1: %f" % (t1-t0)
    # print "t2: %f" % (t2-t1)
    # print "t3: %f" % (t3-t2)
    # print "t4: %f" % (t4-t3)
    # print "t5: %f" % (t5-t4)
    # print "t6: %f" % (t6-t5)
    # print "t61: %f" % (t61-t6)
    # print "t7: %f" % (t7-t61)
    # STOP()

  def remove_bond_restraints_in_place(self, bonded_pairs, sites_cart):
    """ Removing bond restraints in place.
    Current limitations:
      * symmetry related cases (1' means symmetry mate of 1):
        -: 2' --- 1 --- 2
        -: 1' --- 1

    Args:
        bonded_pairs (list): list of bonded i_seqs to remove bonds, e.g.
        [(0,1), (5,15)] for bonds between atoms 0--1 and 5--15 to be removed.
    """
    def _show_bpt(bpt):
      print('bpt size:', bpt.size())
      for i_seq in range(bpt.size()):
        print('i_seq', i_seq)
        d = bpt[i_seq]
        for k in d.keys():
          print('  j_seq: %d, %f' % (k, d[k].distance_ideal))

    # print("1-2,"*30)
    # print('BEGIN shell_sym_tables[0]')
    # self.shell_sym_tables[0].show()
    filtered_bonded_pairs = [sorted(x) for x in bonded_pairs if self.is_bonded_atoms(x[0], x[1])]
    # print("Original BPT")
    # _show_bpt(self.bond_params_table)

    # Removing bonds from self.shell_sym_tables
    asu_mappings = self.crystal_symmetry.special_position_settings().\
        asu_mappings(buffer_thickness=5)
    asu_mappings.process_sites_cart(
        original_sites=sites_cart,
        site_symmetry_table=self.site_symmetry_table)

    pair_asu_table = crystal.pair_asu_table(asu_mappings=asu_mappings)
    for i, j in filtered_bonded_pairs:
      del self.shell_sym_tables[0][i][j]
    pair_asu_table.add_pair_sym_table(self.shell_sym_tables[0])
    shell_asu_tables = crystal.coordination_sequences.shell_asu_tables(
      pair_asu_table=pair_asu_table,
      max_shell=3)
    self.shell_sym_tables = [shell_asu_table.extract_pair_sym_table()
      for shell_asu_table in shell_asu_tables]

    self.reset_internals()
    self.pair_proxies(sites_cart=sites_cart)

  def is_bonded_atoms(self, i_seq, j_seq):
    i_s = i_seq
    j_s = j_seq
    if i_seq > j_seq:
      i_s = j_seq
      j_s = i_seq
    return j_s in self.shell_sym_tables[0][i_s].keys()

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
     reference_coordinate_proxies,
     reference_dihedral_manager,
     ncs_dihedral_manager,
     den_manager,
     chirality_proxies,
     planarity_proxies,
     parallelity_proxies,
     ramachandran_manager) = [None]*13
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
    if flags.reference_coordinate:
      reference_coordinate_proxies = self.reference_coordinate_proxies
    if (flags.reference_dihedral):
      reference_dihedral_manager = self.reference_dihedral_manager
    if (flags.ncs_dihedral): ncs_dihedral_manager = self.ncs_dihedral_manager
    if flags.den_restraints: den_manager = self.den_manager
    if (flags.chirality): chirality_proxies = self.chirality_proxies
    if (flags.planarity): planarity_proxies = self.planarity_proxies
    if (flags.parallelity): parallelity_proxies = self.parallelity_proxies
    if flags.ramachandran_restraints:
      ramachandran_manager = self.ramachandran_manager
    return geometry_restraints.energies.energies(
      sites_cart=sites_cart,
      bond_proxies=bond_proxies,
      nonbonded_proxies=nonbonded_proxies,
      nonbonded_function=nonbonded_function,
      angle_proxies=angle_proxies,
      dihedral_proxies=dihedral_proxies,
      reference_coordinate_proxies=reference_coordinate_proxies,
      reference_dihedral_manager=reference_dihedral_manager,
      ncs_dihedral_manager=ncs_dihedral_manager,
      den_manager=den_manager,
      chirality_proxies=chirality_proxies,
      planarity_proxies=planarity_proxies,
      parallelity_proxies=parallelity_proxies,
      ramachandran_manager = ramachandran_manager,
      external_energy_function=external_energy_function,
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

  def update_atom_nonbonded_type(self,
        i_seq,
        nonbonded_type,
        charge=0):
    if (self.nonbonded_types is not None):
      self.nonbonded_types[i_seq] = nonbonded_type
    if (self.nonbonded_charges is not None):
      self.nonbonded_charges[i_seq] = charge

  def write_geo_file(self,
      sites_cart=None,
      site_labels=None,
      file_name=None,
      file_descriptor=sys.stdout,
      header="# Geometry restraints\n"):
    outf_descriptor = None
    if file_name is None:
      outf_descriptor = file_descriptor
    else:
      outf_descriptor = open(file_name, "w")
    print(header, file=outf_descriptor)
    self.show_sorted(
      sites_cart=sites_cart,
      site_labels=site_labels,
      f=outf_descriptor)
    if file_name is not None:
      outf_descriptor.close()

  def show_sorted(self,
        flags=None,
        sites_cart=None,
        site_labels=None,
        f=None):
    from cctbx.geometry_restraints.auto_linking_types import covalent_headers
    from cctbx.geometry_restraints.auto_linking_types import internal_labels
    covalent_key = 'covalent geometry'
    default_origin_id = origin_ids.get_origin_id(covalent_key)
    if (f is None): f = sys.stdout
    pair_proxies = self.pair_proxies(flags=flags, sites_cart=sites_cart)
    if (sites_cart is None):
      sites_cart = self._sites_cart_used_for_pair_proxies

    if pair_proxies.bond_proxies is not None:
      #
      # write covalent bonds
      #
      label='%s | %s |' % (covalent_headers[0], covalent_key)
      tempbuffer = StringIO()
      pair_proxies.bond_proxies.show_sorted(
          by_value="residual",
          sites_cart=sites_cart,
          site_labels=site_labels,
          f=tempbuffer,
          origin_id=default_origin_id)
      if tempbuffer.getvalue().find(': 0')==-1:
        print(label, tempbuffer.getvalue()[5:], file=f)
      #
      # write bonds with other origin_id
      #
      label=covalent_headers[0]
      for key in origin_ids.get_bond_origin_id_labels():
        origin_id=origin_ids.get_origin_id(key)
        if origin_id==default_origin_id: continue
        t_label=origin_ids.get_geo_file_header(key)
        tempbuffer = StringIO()
        pair_proxies.bond_proxies.show_sorted(
            by_value="residual",
            sites_cart=sites_cart,
            site_labels=site_labels,
            f=tempbuffer,
            prefix="",
            origin_id=origin_id)
        if tempbuffer.getvalue().find(': 0')==-1:
          print('%s | %s |' % (label, t_label), tempbuffer.getvalue()[5:], file=f)
    #
    # write of the other internals for each origin_id
    #
    for i, (proxies, keys, start) in enumerate([
      (
       self.angle_proxies, # self.get_all_angle_proxies(),
       origin_ids.get_angle_origin_id_labels(),
       11),
      (
       self.dihedral_proxies, # self.get_dihedral_proxies(),
       origin_ids.get_dihedral_origin_id_labels(),
       15),
      (
       self.chirality_proxies,
       origin_ids.get_chiral_origin_id_labels(),
       10),
      (
       self.planarity_proxies,
       origin_ids.get_plane_origin_id_labels(),
       10),
      (
       self.parallelity_proxies,
       origin_ids.get_parallelity_origin_id_labels(),
       12),
      ]):
      p_label=covalent_headers[i+1]
      internals=internal_labels[i+1]
      if (proxies is not None):
        #
        # "covalent" original ids
        #
        if p_label not in ['Parallelity']: # not default origin for parallelity
          tempbuffer = StringIO()
          proxies.show_sorted(
            by_value="residual",
            sites_cart=sites_cart,
            site_labels=site_labels,
            f=tempbuffer,
            origin_id=default_origin_id)
          print('%s | %s | %s' % (tempbuffer.getvalue()[:start-1],
                                  covalent_key,
                                  tempbuffer.getvalue()[start:]),
                file=f)
        #
        # other origin ids
        #
        for key in keys:
          origin_id=origin_ids.get_origin_id(key)
          if origin_id==default_origin_id: continue
          label=origin_ids.get_geo_file_header(key, internals=internals)
          if label is None: continue
          tempbuffer = StringIO()
          proxies.show_sorted(
              by_value="residual",
              sites_cart=sites_cart,
              site_labels=site_labels,
              f=tempbuffer,
              prefix="",
              origin_id=origin_id)
          if len(tempbuffer.getvalue()) and tempbuffer.getvalue().find(': 0')==-1:
            print('%s | %s | %s' % (p_label, label, tempbuffer.getvalue()[start:]),
                  file=f)
    #
    # not parsed for geo viewer
    #
    for p_label, proxies in [
        ("Reference torsion angle", self.reference_dihedral_manager),
        ("NCS torsion angle",       self.ncs_dihedral_manager),
        ("",                        self.ramachandran_manager),
        ('Reference coordinate',    self.reference_coordinate_proxies),
        ]:
      if proxies is not None:
        if hasattr(proxies, 'show_sorted'):
          proxies.show_sorted(
              by_value="residual",
              sites_cart=sites_cart,
              site_labels=site_labels,
              proxy_label=p_label,
              f=f)
        else:
          def dist2(xyz1, xyz2):
            d2=0
            for i in range(3):
              d2+=(xyz1[i]-xyz2[i])**2
            return d2
          print('Harmonic | Reference coordinate | restraints :%d' % (len(proxies)),
                file=f)
          for proxy in proxies:
            i_seq=proxy.i_seqs[0]
            d2=dist2(proxy.ref_sites, sites_cart[i_seq])
            print('bond %s, d2=%6.3fA, sigma=%7.3f, weight=%7.3f' % (
              site_labels[i_seq],
              d2,
              # '(%9.5f, %9.5f, %9.5f)' % tuple(proxy.ref_sites),
              1/proxy.weight**0.5,
              proxy.weight,
              ),
              file=f,
              )
        print(file=f)
    #
    # Here should be showing DEN manager...
    #
    if (pair_proxies.nonbonded_proxies is not None):
      tempbuffer = StringIO()
      pair_proxies.nonbonded_proxies.show_sorted(
        by_value="delta",
        sites_cart=sites_cart, site_labels=site_labels,
        f=tempbuffer,
        suppress_model_minus_vdw_greater_than=None)
      start=10
      print('%s| unspecified | %s' % (tempbuffer.getvalue()[:start],
                                      tempbuffer.getvalue()[start:]
                                     ),
            file=f)

# This should be in model class?
#  def nb_overlaps_info(
#    self,
#    sites_cart,
#    hd_sel,
#    macro_mol_sel=None,
#    site_labels=None):
#    """ non-bonded overlaps information """
#    from cctbx.geometry_restraints.nonbonded_overlaps import info
#    if not macro_mol_sel:
#      from cctbx.geometry_restraints.nonbonded_overlaps import get_macro_mol_sel
#      macro_mol_sel = get_macro_mol_sel(pdb_processed_file=self)
#
#    return info(
#      geometry_restraints_manager=self,
#      macro_molecule_selection=macro_mol_sel,
#      sites_cart=sites_cart,
#      hd_sel=hd_sel,
#      site_labels=site_labels).result

  def _bond_generator(self):
    simple, asu = self.get_all_bond_proxies()
    simple = simple.get_proxies_without_origin_id(0)
    asu = asu.get_proxies_without_origin_id(0)
    for bond in list(simple)+list(asu):
      yield bond

  def get_struct_conn_mmcif(self, hierarchy):
    atoms = hierarchy.atoms()
    def _atom_info(atom):
      return [  # auth
                atom.parent().resname.strip(),
                atom.parent().parent().parent().id,
                atom.parent().parent().resseq.strip(),
                atom.name.strip(),
                # label
                atom.parent().resname.strip(),
                hierarchy.get_label_asym_id(atom.parent().parent()),
                hierarchy.get_label_seq_id(atom.parent()),
                atom.name.strip(),
                '.', # role
               ]
      return res
    def _atom_info_grouped(bond):
      row = []
      if hasattr(bond, 'i_seqs'):
        i1 = bond.i_seqs[0]
        i2 = bond.i_seqs[1]
      else:
        i1 = bond.i_seq
        i2 = bond.j_seq
      row += _atom_info(atoms[i1])
      row += _atom_info(atoms[i2])
      #row.append('1_555') # symmetry!
      return row
    from cctbx.geometry_restraints.auto_linking_types import origin_ids
    struct_conn_loop = iotbx.cif.model.loop(header=(
      '_struct_conn.id',
      '_struct_conn.conn_type_id',
      '_struct_conn.ptnr1_auth_comp_id',
      '_struct_conn.ptnr1_auth_asym_id',
      '_struct_conn.ptnr1_auth_seq_id',
      '_struct_conn.ptnr1_auth_atom_id',
      '_struct_conn.ptnr1_label_comp_id',
      '_struct_conn.ptnr1_label_asym_id',
      '_struct_conn.ptnr1_label_seq_id',
      '_struct_conn.ptnr1_label_atom_id',
      '_struct_conn.ptnr1_role',
      #'_struct_conn.ptnr1_symmetry',
      '_struct_conn.ptnr2_auth_comp_id',
      '_struct_conn.ptnr2_auth_asym_id',
      '_struct_conn.ptnr2_auth_seq_id',
      '_struct_conn.ptnr2_auth_atom_id',
      '_struct_conn.ptnr2_label_comp_id',
      '_struct_conn.ptnr2_label_asym_id',
      '_struct_conn.ptnr2_label_seq_id',
      '_struct_conn.ptnr2_label_atom_id',
      '_struct_conn.ptnr2_role',
      #'_struct_conn.ptnr2_symmetry',
      '_struct_conn.details',
      ))
    for i, bond in enumerate(self._bond_generator()):
      row = ['C%05d' % (i+1)]
      origin_id_info = origin_ids[0].get(bond.origin_id, None)
      assert origin_id_info
      if origin_id_info[0]=='SS BOND': row.append('disulf')
      elif origin_id_info[0]=='metal coordination': row.append('metalc')
      elif origin_id_info[0]=='hydrogen bonds': row.append('hydrog')
      else: row.append('covale')
      row += _atom_info_grouped(bond)
      if len(origin_id_info)>2 and origin_id_info[2]:
        row.append(origin_id_info[2])
      elif origin_id_info[1]: row.append(origin_id_info[1])
      else: row.append('.') # details
      struct_conn_loop.add_row(row)
    return struct_conn_loop

  def get_cif_link_entries(self, mon_lib_srv):
    from cctbx.geometry_restraints.auto_linking_types import origin_ids
    links = iotbx.cif.model.cif()
    # not sure why this is the case but these are missing from origin_ids
    done = ['POST-BETA-TRANS',
            'PRE-BETA-TRANS',
            'SSRAD',
            'rna2p',
            'rna3p',
            ]
    for i, bond in enumerate(self._bond_generator()):
      row = ['C%05d' % (i+1)]
      origin_id_info = origin_ids[0].get(bond.origin_id, None)
      assert origin_id_info
      key = origin_id_info[0].replace('link_', '')
      link_key = 'link_%s' % origin_id_info[0]
      if origin_id_info[0]=='SS BOND':
        links['link_SS'] = mon_lib_srv.link_link_id_dict['SS'].as_cif_block()
      elif origin_id_info[0]=='User supplied cif_link':
        tlinks = []
        for origin_id in origin_ids:
          for oi, item in origin_id.items():
            tlinks.append(item[0])
          break
        for tkey, item in sorted(mon_lib_srv.link_link_id_dict.items()):
          if tkey in done: continue
          done.append(tkey)
          if tkey not in tlinks and 'link_%s' % tkey not in tlinks:
            links['link_%s' % tkey] = item.as_cif_block()
      elif key in mon_lib_srv.link_link_id_dict:
        links['link_%s' % key] = mon_lib_srv.link_link_id_dict[key].as_cif_block()
      elif origin_id_info[0] in ['hydrogen bonds',
                                 'edits',
                                 'metal coordination',
                                 'glycosidic custom',
                                 'Misc. bond',
                                 ]:
        # not writing cif_link for various reasons like programatic links
        pass
      else:
        print(origin_id_info[0])
        assert 0
    return links

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
  new_grm = manager(
    crystal_symmetry=box.crystal_symmetry(),
    site_symmetry_table=asu_mappings.site_symmetry_table(),
    bond_params_table=bond_params_table,
    shell_sym_tables=shell_sym_tables,
    nonbonded_params=nonbonded_params,
    nonbonded_types=nonbonded_types,
    nonbonded_function=geometry_restraints.prolsq_repulsion_function(),
    max_reasonable_bond_distance=max_reasonable_bond_distance)
  return box.sites_cart, new_grm

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

from __future__ import absolute_import, division, print_function

class Base_geometry(object):
  """
  Base class for geometry restraints manager. This will go inside
  restraints manager defined in mmtbx/restraints.py
  The reason for it is for external packages, like AMBER to be able
  to inherit from it and be placed instead of classical geometry.
  Therefore all functions defined for classical geometry are defined
  here with hopefully reasonable default actions.
  """
  def __init__(self):
    self._source = None

  def set_source(self, source):
    assert self._source is None or self._source == source
    self._source = source

  def get_source(self):
    return self._source

  def reset_internals(self):
    pass
  def replace_site_symmetry(self, new_site_symmetry_table):
    pass
  def simple_edge_list(self, omit_slack_greater_than=0):
    return None
  def rigid_clusters_due_to_dihedrals_and_planes(self,
      constrain_dihedrals_with_sigma_less_than):
    return None
  def construct_tardy_tree(self,
      sites=None,
      sites_cart=None,
      selection=None,
      omit_bonds_with_slack_greater_than=0,
      constrain_dihedrals_with_sigma_less_than=10,
      near_singular_hinges_angular_tolerance_deg=5):
    return None
  def reduce_for_tardy(self,
        tardy_tree,
        omit_bonds_with_slack_greater_than=0,
        include_den_restraints=False):
    return self
  def sites_cart_used_for_pair_proxies(self):
    return None
  def new_including_isolated_sites(self,
        n_additional_sites,
        model_indices=None,
        conformer_indices=None,
        sym_excl_indices=None,
        donor_acceptor_excl_groups=None,
        site_symmetry_table=None,
        nonbonded_types=None,
        nonbonded_charges=None):
    return self

  def select(self, selection=None, iselection=None):
    """
    This is cruicial function and should be defined in child.
    Given selection return manager only with selected atoms and restraints
    """
    raise NotImplementedError
  def shift_sites_cart(self, shift):
    """
    Sometimes coordinates of the model are shifted (e.g. in boxing with map).
    In case any element of restraints cares, do what is appropriate here.
    """
    raise NotImplementedError

  def discard_symmetry(self, new_unit_cell):
    return self
  def add_angles_in_place(self, additional_angle_proxies):
    pass
  def remove_angles_in_place(self, selection):
    pass
  def get_user_supplied_restraints(self):
    return None, None, None, None, None
  def remove_user_supplied_restraints_in_place(self):
    raise NotImplementedError
  def get_bond_proxies_without_user_supplied(self, sites_cart=None):
    return None, None
  def get_angle_proxies_without_user_supplied(self):
    return None
  def get_planarity_proxies_without_user_supplied(self):
    return None
  def get_parallelity_proxies_without_user_supplied(self):
    return None

  #=================================================================
  # Reference coordinate proxies methods
  #=================================================================
  def get_reference_coordinate_proxies(self):
    return None
  def adopt_reference_coordinate_restraints_in_place(self,
      reference_coordinate_proxies):
    pass
  def remove_reference_coordinate_restraints_in_place(self,
      selection=None):
    pass
  def get_n_reference_coordinate_proxies(self):
    return 0

  def append_reference_coordinate_restraints_in_place(self,
      reference_coordinate_proxies):
    pass
  def add_reference_coordinate_restraints_in_place(self,
      all_chain_proxies=None,
      pdb_hierarchy=None,
      selection=None,
      exclude_outliers=True,
      sigma=0.2,
      limit=1.0,
      top_out=False):
    pass

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
    pass
  def remove_chi_torsion_restraints_in_place(self, selection=None):
    pass
  def get_chi_torsion_proxies(self):
    return None

  def get_n_chi_torsion_proixes(self):
    return None
  #=================================================================
  # Dihedral proxies methods
  #=================================================================
  def add_dihedrals_in_place(self,
      additional_dihedral_proxies,
      check_for_duplicates=True):
    pass
  def remove_dihedrals_in_place(self, selection):
    pass
  def get_dihedral_proxies(self):
    return None

  #=================================================================
  # C-beta dihedral proxies methods
  #=================================================================
  def get_c_beta_torsion_proxies(self):
    return None
  def get_n_c_beta_torsion_proxies(self):
    return 0

  def remove_c_beta_torsion_restraints_in_place(self, selection=None):
    pass
  #=================================================================
  # Reference dihedral proxies methods
  #=================================================================
  def adopt_reference_dihedral_manager(self, manager):
    pass
  def remove_reference_dihedrals_in_place(self, selection):
    pass
  def remove_ncs_dihedrals_in_place(self):
    pass
  def update_dihedral_ncs_restraints(self, model, log):
    pass
  def get_n_reference_dihedral_proxies(self):
    return None

  def sync_reference_dihedral_with_ncs(self, log):
    pass
  #=================================================================
  # DEN manager/proxies methods
  #=================================================================
  def create_den_manager(self, den_params, pdb_hierarchy, log):
    pass
  def adopt_den_manager(self, den_manager):
    pass
  def get_n_den_proxies(self):
    return 0
  def remove_chiralities_in_place(self, selection):
    pass
  def add_planarities_in_place(self, additional_planarity_proxies):
    pass
  def remove_planarities_in_place(self, selection):
    pass
  def add_parallelities_in_place(self, additional_parallelity_proxies):
    pass
  def remove_parallelities_in_place(self, selection):
    pass
  #=================================================================
  # Ramachandran manager/proxies methods
  #=================================================================
  def set_ramachandran_restraints(self, manager):
    pass
  def update_ramachandran_restraints_phi_psi_targets(self, sites_cart):
    pass
  def remove_ramachandran_in_place(self):
    pass
  def get_n_ramachandran_proxies(self):
    return 0
  #=================================================================
  # Secondary structure manager/proxies methods
  #=================================================================
  def set_secondary_structure_restraints(self, ss_manager, hierarchy, log):
    pass
  def remove_secondary_structure_restraints(self):
    # Not implemented. The problem here is to remove hbond restraints, which
    # requires modification of pair_proxies and as complicated as addition
    # of bond restraint.
    raise NotImplementedError
  def set_external_energy_function(self, energy_function):
    pass
  def _get_n_bond_proxies_origin(self, origin_id):
    return 0
  def get_n_bond_proxies(self):
    return 0
  def get_covalent_bond_proxies(self, sites_cart=None):
    return None
  def get_all_bond_proxies(self, sites_cart=None):
    return None
  def get_covalent_angle_proxies(self):
    return None
  def get_all_angle_proxies(self):
    return None
  def get_n_hbond_proxies(self):
    return None
  def get_n_angle_proxies(self):
    return 0
  def get_n_hangle_proxies(self):
    return 0
  def get_n_stacking_proxies(self):
    return 0
  def get_n_parallelity_bp_proxies(self):
    return 0
  def get_n_planarity_proxies(self):
    return 0
  def get_n_planarity_bp_proxies(self):
    return 0
  def get_hbond_proxies_iseqs(self):
    return None

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
    return self

  def add_new_hbond_restraints_in_place(self, proxies, sites_cart,
      max_distance_between_connecting_atoms=5,
      skip_max_proxy_distance_calculation=False):
    pass

  def add_new_bond_restraints_in_place(self, proxies, sites_cart,
      max_distance_between_connecting_atoms=5,
      skip_max_proxy_distance_calculation=False):
    """ Add new bond restraints for list of proxies to this
    geometry restraints manager, _in_place_! Returns nothing.
    proxies - list of bond_proxy objects. The symmetry operation for the
    paired atoms is determined here, therefore the proxy.rt_mx_ji may be
    anything."""
    pass

  def is_bonded_atoms(self, i_seq, j_seq):
    return True

  def pair_proxies(self,
        sites_cart=None,
        flags=None,
        asu_is_inside_epsilon=None,
        bonded_distance_cutoff_epsilon=None,
        site_labels=None):
    return None

  def nonbonded_model_distances(self, sites_cart=None):
    return None
  def update_plain_pair_sym_table(self, sites_frac):
    pass

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
    """
    Crucial function to override in child. Maybe we should put some info
    here later.
    """
    raise NotImplementedError

  def harmonic_restraints(self, variables, type_indices, type_weights):
    """
    Not clear what it is
    """
    raise NotImplementedError

  def ta_harmonic_restraints(self, sites_cart, ta_harmonic_restraint_info, weight = 0.001, slack = 0.5):
    """
    Not clear what it is
    """
    raise NotImplementedError

  def update_atom_nonbonded_type(self,
        i_seq,
        nonbonded_type,
        charge=0):
    pass

  def write_geo_file(self,
      sites_cart=None,
      site_labels=None,
      file_name=None,
      file_descriptor=None,
      header="# Geometry restraints\n"):
    """
    Define this if you want the .geo file
    """
    pass

  def show_sorted(self,
        flags=None,
        sites_cart=None,
        site_labels=None,
        f=None):
    """
    Similar to geo file
    """
    pass

  def nb_overlaps_info(
    self,
    sites_cart,
    hd_sel,
    macro_mol_sel=None,
    site_labels=None):
    """ non-bonded overlaps information """
    return None

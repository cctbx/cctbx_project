from cctbx import restraints
from cctbx.restraints import energies
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args

class manager:

  def __init__(self,
        crystal_symmetry=None,
        site_symmetry_table=None,
        bond_params_table=None,
        shell_sym_tables=None,
        repulsion_params=None,
        repulsion_types=None,
        nonbonded_distance_cutoff=5,
        nonbonded_buffer=1,
        angle_proxies=None,
        dihedral_proxies=None,
        chirality_proxies=None,
        planarity_proxies=None):
    if (site_symmetry_table is not None): assert crystal_symmetry is not None
    if (bond_params_table is not None and site_symmetry_table is not None):
      assert bond_params_table.size() == site_symmetry_table.indices().size()
    if (shell_sym_tables is not None and site_symmetry_table is not None):
      assert len(shell_sym_tables) > 0
      assert shell_sym_tables[0].size() == site_symmetry_table.indices().size()
    if (repulsion_types is not None and site_symmetry_table is not None):
      assert repulsion_types.size() == site_symmetry_table.indices().size()
    adopt_init_args(self, locals())
    self._sites_cart_used_for_pair_proxies = None
    self._pair_proxies = None
    self.n_updates_pair_proxies = 0

  def pair_proxies(self,
        sites_cart=None,
        lock=00000,
        asu_is_inside_epsilon=None):
    if (self.repulsion_types is None):
      if (self._pair_proxies is None):
        self.n_updates_pair_proxies += 1
        self._pair_proxies = restraints.pair_proxies(
          bond_params_table=self.bond_params_table)
    elif (sites_cart is not None
          and (self._sites_cart_used_for_pair_proxies is None
          or (not lock
          and self._sites_cart_used_for_pair_proxies.max_distance(sites_cart)
              > self.nonbonded_buffer))):
      self.n_updates_pair_proxies += 1
      self._sites_cart_used_for_pair_proxies = sites_cart.deep_copy()
      if (self.crystal_symmetry is None):
        asu_mappings = \
          crystal.direct_space_asu.non_crystallographic_asu_mappings(
            sites_cart=sites_cart)
        bonded_distance_cutoff = 0
      else:
        unit_cell = self.crystal_symmetry.unit_cell()
        sites_frac = unit_cell.fractionalization_matrix() * sites_cart
        bonded_distance_cutoff = 0
        for shell_sym_table in self.shell_sym_tables:
          bonded_distance_cutoff = max(bonded_distance_cutoff,
            flex.max(crystal.get_distances(
              pair_sym_table=shell_sym_table,
              orthogonalization_matrix=unit_cell.orthogonalization_matrix(),
              sites_frac=sites_frac)))
        asu_mappings = crystal.symmetry.asu_mappings(self.crystal_symmetry,
          buffer_thickness=max(bonded_distance_cutoff,
                               self.nonbonded_distance_cutoff),
          is_inside_epsilon=asu_is_inside_epsilon)
        asu_mappings.process_sites_frac(
          original_sites=sites_frac,
          site_symmetry_table=self.site_symmetry_table)
      shell_asu_tables = [
        crystal.pair_asu_table(asu_mappings=asu_mappings)
          .add_pair_sym_table(sym_table=shell_sym_table)
            for shell_sym_table in self.shell_sym_tables]
      self._pair_proxies = restraints.pair_proxies(
        repulsion_params=self.repulsion_params,
        repulsion_types=self.repulsion_types,
        bond_params_table=self.bond_params_table,
        shell_asu_tables=shell_asu_tables,
        bonded_distance_cutoff=bonded_distance_cutoff,
        nonbonded_distance_cutoff=self.nonbonded_distance_cutoff,
        nonbonded_buffer=self.nonbonded_buffer)
    elif (self._pair_proxies is None):
      raise AssertionError("pair_proxies not defined already.")
    return self._pair_proxies

  def energies(self,
        sites_cart,
        compute_gradients=00000,
        disable_asu_cache=00000,
        lock_pair_proxies=00000):
    pair_proxies = self.pair_proxies(
      sites_cart=sites_cart,
      lock=lock_pair_proxies)
    return energies.energies(
      sites_cart=sites_cart,
      bond_proxies=pair_proxies.bond_proxies,
      repulsion_proxies=pair_proxies.repulsion_proxies,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      compute_gradients=compute_gradients,
      disable_asu_cache=disable_asu_cache)

from cctbx import geometry_restraints
import cctbx.geometry_restraints.flags
import cctbx.geometry_restraints.energies
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import adopt_init_args
import sys

class manager:

  def __init__(self,
        crystal_symmetry=None,
        model_indices=None,
        conformer_indices=None,
        site_symmetry_table=None,
        bond_params_table=None,
        shell_sym_tables=None,
        nonbonded_params=None,
        nonbonded_types=None,
        nonbonded_function=None,
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
    if (nonbonded_types is not None and site_symmetry_table is not None):
      assert nonbonded_types.size() == site_symmetry_table.indices().size()
    adopt_init_args(self, locals())
    self._sites_cart_used_for_pair_proxies = None
    self._flags_bond_used_for_pair_proxies = False
    self._flags_nonbonded_used_for_pair_proxies = False
    self._pair_proxies = None
    self.n_updates_pair_proxies = 0

  def pair_proxies(self,
        sites_cart=None,
        flags=None,
        lock=False,
        asu_is_inside_epsilon=None,
        bonded_distance_cutoff_epsilon=None):
    if (bonded_distance_cutoff_epsilon is None):
      bonded_distance_cutoff_epsilon = 1.e-6
    if (self.nonbonded_types is None):
      if (self._pair_proxies is None):
        self.n_updates_pair_proxies += 1
        self._pair_proxies = geometry_restraints.pair_proxies(
          flags=flags,
          bond_params_table=self.bond_params_table)
    elif (sites_cart is not None
          and (self._sites_cart_used_for_pair_proxies is None
               or flags is not None
                  and (self._flags_bond_used_for_pair_proxies
                          != flags.bond
                       or self._flags_nonbonded_used_for_pair_proxies
                             != flags.nonbonded)
          or (not lock
          and self._sites_cart_used_for_pair_proxies.max_distance(sites_cart)
              > self.nonbonded_buffer))):
      self.n_updates_pair_proxies += 1
      self._sites_cart_used_for_pair_proxies = sites_cart.deep_copy()
      if (flags is None):
        self._flags_bond_used_for_pair_proxies = True
        self._flags_nonbonded_used_for_pair_proxies = True
      else:
        self._flags_bond_used_for_pair_proxies = flags.bond
        self._flags_nonbonded_used_for_pair_proxies = flags.nonbonded
      bonded_distance_cutoff = 0
      if (self.crystal_symmetry is None):
        for shell_sym_table in self.shell_sym_tables:
          bonded_distance_cutoff = max(bonded_distance_cutoff,
            flex.max(crystal.get_distances(
              pair_sym_table=shell_sym_table,
              sites_cart=sites_cart)))
        bonded_distance_cutoff *= (1 + bonded_distance_cutoff_epsilon)
        asu_mappings = \
          crystal.direct_space_asu.non_crystallographic_asu_mappings(
            sites_cart=sites_cart)
      else:
        unit_cell = self.crystal_symmetry.unit_cell()
        sites_frac = unit_cell.fractionalization_matrix() * sites_cart
        for shell_sym_table in self.shell_sym_tables:
          bonded_distance_cutoff = max(bonded_distance_cutoff,
            flex.max(crystal.get_distances(
              pair_sym_table=shell_sym_table,
              orthogonalization_matrix=unit_cell.orthogonalization_matrix(),
              sites_frac=sites_frac)))
        bonded_distance_cutoff *= (1 + bonded_distance_cutoff_epsilon)
        asu_mappings = crystal.symmetry.asu_mappings(self.crystal_symmetry,
          buffer_thickness=max(
            bonded_distance_cutoff,
            self.nonbonded_distance_cutoff + self.nonbonded_buffer),
          is_inside_epsilon=asu_is_inside_epsilon)
        asu_mappings.process_sites_frac(
          original_sites=sites_frac,
          site_symmetry_table=self.site_symmetry_table)
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
        nonbonded_params=self.nonbonded_params,
        nonbonded_types=self.nonbonded_types,
        nonbonded_distance_cutoff_plus_buffer=
          self.nonbonded_distance_cutoff+self.nonbonded_buffer)
    elif (self._pair_proxies is None):
      raise AssertionError("pair_proxies not defined already.")
    return self._pair_proxies

  def energies(self,
        sites_cart,
        flags=None,
        compute_gradients=False,
        disable_asu_cache=False,
        lock_pair_proxies=False):
    if (flags is None):
      flags = geometry_restraints.flags.flags(default=True)
    pair_proxies = self.pair_proxies(
      flags=flags,
      sites_cart=sites_cart,
      lock=lock_pair_proxies)
    (bond_proxies,
     nonbonded_proxies,
     nonbonded_function,
     angle_proxies,
     dihedral_proxies,
     chirality_proxies,
     planarity_proxies) = [None]*7
    if (flags.bond):
      assert pair_proxies.bond_proxies is not None
      bond_proxies = pair_proxies.bond_proxies
    if (flags.nonbonded and self.nonbonded_types is not None):
      assert pair_proxies.nonbonded_proxies is not None
      nonbonded_proxies = pair_proxies.nonbonded_proxies
      nonbonded_function = self.nonbonded_function
    if (flags.angle):     angle_proxies = self.angle_proxies
    if (flags.dihedral):  dihedral_proxies = self.dihedral_proxies
    if (flags.chirality): chirality_proxies = self.chirality_proxies
    if (flags.planarity): planarity_proxies = self.planarity_proxies
    return geometry_restraints.energies.energies(
      sites_cart=sites_cart,
      bond_proxies=bond_proxies,
      nonbonded_proxies=nonbonded_proxies,
      nonbonded_function=nonbonded_function,
      angle_proxies=angle_proxies,
      dihedral_proxies=dihedral_proxies,
      chirality_proxies=chirality_proxies,
      planarity_proxies=planarity_proxies,
      compute_gradients=compute_gradients,
      disable_asu_cache=disable_asu_cache)

  def show_interactions(self, flags=None, sites_cart=None, i_seq=None, f=None):
    if (f is None): f = sys.stdout
    pair_proxies = self.pair_proxies(flags=flags, sites_cart=sites_cart)
    if (pair_proxies.bond_proxies is not None):
      for proxy in pair_proxies.bond_proxies.simple:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "bond simple:", proxy.i_seqs
      for proxy in pair_proxies.bond_proxies.asu:
        if (i_seq is None or i_seq in [proxy.i_seq, proxy.j_seq]):
          print >> f, "bond asu:", (proxy.i_seq, proxy.j_seq, proxy.j_sym)
    if (self.angle_proxies is not None):
      for proxy in self.angle_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "angle:", proxy.i_seqs
    if (self.dihedral_proxies is not None):
      for proxy in self.dihedral_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "dihedral:", proxy.i_seqs
    if (self.chirality_proxies is not None):
      for proxy in self.chirality_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "chirality:", proxy.i_seqs
    if (self.planarity_proxies is not None):
      for proxy in self.planarity_proxies:
        if (i_seq is None or i_seq in proxy.i_seqs):
          print >> f, "planarity:", tuple(proxy.i_seqs)

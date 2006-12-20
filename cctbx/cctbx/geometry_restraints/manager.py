from cctbx import geometry_restraints
import cctbx.geometry_restraints.flags
import cctbx.geometry_restraints.energies
from cctbx import crystal
from cctbx.array_family import flex
from scitbx.python_utils.misc import store
from libtbx import introspection
from libtbx import adopt_init_args
import sys

class manager(object):

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
        nonbonded_distance_cutoff=None,
        nonbonded_buffer=1,
        angle_proxies=None,
        dihedral_proxies=None,
        chirality_proxies=None,
        planarity_proxies=None,
        plain_pairs_radius=None,
        max_reasonable_bond_distance=None):
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
    self.plain_pair_sym_table = None
    self.nonbonded_distance_cutoff_was_determined_automatically = False
    self.adjusted_nonbonded_distance_cutoff = self.nonbonded_distance_cutoff
    self.effective_nonbonded_buffer = self.nonbonded_buffer
    self.n_updates_pair_proxies = 0

  def sites_cart_used_for_pair_proxies(self):
    return self._sites_cart_used_for_pair_proxies

  def new_including_isolated_sites(self,
        n_additional_sites,
        model_indices=None,
        conformer_indices=None,
        site_symmetry_table=None,
        nonbonded_types=None):
    assert n_additional_sites >= 0
    assert (model_indices is None) == (self.model_indices is None)
    assert (conformer_indices is None) == (self.conformer_indices is None)
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
      chirality_proxies=self.chirality_proxies,
      planarity_proxies=self.planarity_proxies,
      plain_pairs_radius=self.plain_pairs_radius)

  def select(self, selection):
    iselection = selection.iselection()
    selected_model_indices = None
    if (self.model_indices is not None):
      selected_model_indices = self.model_indices.select(
        iselection)
    selected_conformer_indices = None
    if (self.conformer_indices is not None):
      selected_conformer_indices = self.conformer_indices.select(
        iselection)
    selected_site_symmetry_table = None
    if (self.site_symmetry_table is not None):
      selected_site_symmetry_table = self.site_symmetry_table.select(
        iselection)
    selected_bond_params_table = None
    if (self.bond_params_table is not None):
      selected_bond_params_table = self.bond_params_table.proxy_select(
        iselection)
    selected_shell_sym_tables = None
    if (self.shell_sym_tables is not None):
      selected_shell_sym_tables = [shell_sym_table.proxy_select(iselection)
        for shell_sym_table in self.shell_sym_tables]
    selected_nonbonded_types = None
    if (self.nonbonded_types is not None):
      selected_nonbonded_types = self.nonbonded_types.select(
        iselection)
    selected_angle_proxies = None
    if (self.angle_proxies is not None):
      selected_angle_proxies = self.angle_proxies.proxy_select(
        selection.size(), iselection)
    selected_dihedral_proxies = None
    if (self.dihedral_proxies is not None):
      selected_dihedral_proxies = self.dihedral_proxies.proxy_select(
        selection.size(), iselection)
    selected_chirality_proxies = None
    if (self.chirality_proxies is not None):
      selected_chirality_proxies = self.chirality_proxies.proxy_select(
        selection.size(), iselection)
    selected_planarity_proxies = None
    if (self.planarity_proxies is not None):
      selected_planarity_proxies = self.planarity_proxies.proxy_select(
        selection.size(), iselection)
    return manager(
      crystal_symmetry=self.crystal_symmetry,
      model_indices=selected_model_indices,
      conformer_indices=selected_conformer_indices,
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
      chirality_proxies=selected_chirality_proxies,
      planarity_proxies=selected_planarity_proxies,
      plain_pairs_radius=self.plain_pairs_radius)

  def remove_angles_in_place(self, selection):
    self.angle_proxies = self.angle_proxies.proxy_remove(
      selection=selection)

  def remove_dihedrals_in_place(self, selection):
    self.dihedral_proxies = self.dihedral_proxies.proxy_remove(
      selection=selection)

  def remove_chiralities_in_place(self, selection):
    self.chirality_proxies = self.chirality_proxies.proxy_remove(
      selection=selection)

  def remove_planarities_in_place(self, selection):
    self.planarity_proxies = self.planarity_proxies.proxy_remove(
      selection=selection)

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
              > self.effective_nonbonded_buffer))):
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
        self.nonbonded_distance_cutoff \
          = self.nonbonded_params.find_max_vdw_distance(
              nonbonded_types=self.nonbonded_types)
        self.nonbonded_distance_cutoff_was_determined_automatically = True
        self.adjusted_nonbonded_distance_cutoff \
          = self.nonbonded_distance_cutoff
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
            if (self.max_reasonable_bond_distance is not None
                and   bonded_distance_cutoff
                    > self.max_reasonable_bond_distance):
              raise RuntimeError(
                "Bond distance > max_reasonable_bond_distance: %.6g > %.6g" % (
                  bonded_distance_cutoff, self.max_reasonable_bond_distance))
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
          nonbonded_params=self.nonbonded_params,
          nonbonded_types=self.nonbonded_types,
          nonbonded_distance_cutoff_plus_buffer
            =current_nonbonded_distance_cutoff_plus_buffer)
        introspection.virtual_memory_info().update_max()
        if (self._pair_proxies.nonbonded_proxies is None):
          break
        self.adjusted_nonbonded_distance_cutoff = \
          max(0, self._pair_proxies.nonbonded_proxies.max_vdw_distance)
        if (self.nonbonded_distance_cutoff
            < self._pair_proxies.nonbonded_proxies.max_vdw_distance):
          if (self.nonbonded_distance_cutoff_was_determined_automatically):
            raise AssertionError("Internal error.")
          raise AssertionError(
            "nonbonded_distance_cutoff=%.6g is too small:"
            " max_vdw_distance=%.6g" % (
              self.nonbonded_distance_cutoff,
              self._pair_proxies.nonbonded_proxies.max_vdw_distance))
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
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False,
        lock_pair_proxies=False,
        normalization=False):
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
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization)

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
          sorted_asu_proxies=pair_proxies.nonbonded_proxies,
          function=self.nonbonded_function)
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

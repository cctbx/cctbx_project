from libtbx import adopt_init_args
import cctbx.adp_restraints.flags
import cctbx.adp_restraints.energies
import cctbx.geometry_restraints.energies
import cctbx.geometry_restraints.flags
import cctbx.geometry_restraints
import sys

class manager(object):

  def __init__(self,
        bond_params_table=None,
        angle_proxies=None,
        dihedral_proxies=None,
        chirality_proxies=None,
        planarity_proxies=None,
        bond_similarity_proxies=None,
        adp_similarity_proxies=None,
        rigid_bond_proxies=None,
        isotropic_adp_proxies=None):
    adopt_init_args(self, locals())
    self._sites_cart_used_for_pair_proxies = None
    self._flags_bond_used_for_pair_proxies = False
    self._pair_proxies = None
    self.n_updates_pair_proxies = 0

  def pair_proxies(self, flags=None):
    if (self._pair_proxies is None):
      self.n_updates_pair_proxies += 1
      self._pair_proxies = cctbx.geometry_restraints.pair_proxies(
        flags=flags,
        bond_params_table=self.bond_params_table)
    return self._pair_proxies

  def show_sorted(self,
        flags=None,
        unit_cell=None,
        sites_cart=None,
        site_labels=None,
        u_cart=None,
        u_iso=None,
        use_u_aniso=None,
        f=None,
        prefix="",
        max_items=None):
    if (f is None): f = sys.stdout
    pair_proxies = self.pair_proxies(flags=flags)
    if (sites_cart is None):
      sites_cart = self._sites_cart_used_for_pair_proxies
    if (pair_proxies.bond_proxies is not None):
      pair_proxies.bond_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.angle_proxies is not None):
      self.angle_proxies.show_sorted(
        by_value="residual",
        unit_cell=unit_cell, sites_cart=sites_cart, site_labels=site_labels,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.dihedral_proxies is not None):
      self.dihedral_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels,
        unit_cell=unit_cell, f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.chirality_proxies is not None):
      self.chirality_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.planarity_proxies is not None):
      self.planarity_proxies.show_sorted(
        by_value="residual",
        unit_cell=unit_cell, sites_cart=sites_cart, site_labels=site_labels,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.bond_similarity_proxies is not None):
      self.bond_similarity_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels,
        unit_cell=unit_cell, f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.adp_similarity_proxies is not None):
      self.adp_similarity_proxies.show_sorted(
        by_value="residual", site_labels=site_labels,
        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.rigid_bond_proxies is not None):
      self.rigid_bond_proxies.show_sorted(
        by_value="residual",
        sites_cart=sites_cart, site_labels=site_labels, u_cart=u_cart,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.isotropic_adp_proxies is not None):
      self.isotropic_adp_proxies.show_sorted(
        by_value="residual",
        site_labels=site_labels, u_cart=u_cart,
        f=f, prefix=prefix, max_items=max_items)
      print >> f

  def energies_sites(self,
        sites_cart,
        unit_cell=None,
        flags=None,
        compute_gradients=False,
        gradients=None,
        disable_asu_cache=False,
        normalization=False):
    if (flags is None):
      flags = cctbx.geometry_restraints.flags.flags(default=True)
    pair_proxies = self.pair_proxies(flags=flags)
    (bond_proxies,
     angle_proxies,
     dihedral_proxies,
     chirality_proxies,
     planarity_proxies,
     bond_similarity_proxies) = [None]*6
    if (flags.bond):
      assert pair_proxies.bond_proxies is not None
      bond_proxies = pair_proxies.bond_proxies
    if (flags.angle):     angle_proxies = self.angle_proxies
    if (flags.dihedral):  dihedral_proxies = self.dihedral_proxies
    if (flags.chirality): chirality_proxies = self.chirality_proxies
    if (flags.planarity): planarity_proxies = self.planarity_proxies
    if (flags.bond_similarity):
      bond_similarity_proxies = self.bond_similarity_proxies
    return cctbx.geometry_restraints.energies.energies(
      sites_cart=sites_cart,
      unit_cell=unit_cell,
      bond_proxies=bond_proxies,
      angle_proxies=angle_proxies,
      dihedral_proxies=dihedral_proxies,
      chirality_proxies=chirality_proxies,
      planarity_proxies=planarity_proxies,
      bond_similarity_proxies=bond_similarity_proxies,
      compute_gradients=compute_gradients,
      gradients=gradients,
      disable_asu_cache=disable_asu_cache,
      normalization=normalization)

  def energies_adps(self,
        sites_cart,
        u_cart,
        u_iso=None,
        use_u_aniso=None,
        flags=None,
        compute_gradients=False,
        gradients_aniso_cart=None,
        gradients_iso=None,
        normalization=False):
    if (flags is None):
      flags = cctbx.adp_restraints.flags.flags(default=True)
    (adp_similarity_proxies,
     rigid_bond_proxies,
     isotropic_adp_proxies) = [None]*3
    if (flags.adp_similarity): adp_similarity_proxies = self.adp_similarity_proxies
    if (flags.rigid_bond): rigid_bond_proxies = self.rigid_bond_proxies
    if (flags.isotropic_adp): isotropic_adp_proxies = self.isotropic_adp_proxies
    return cctbx.adp_restraints.energies.energies(
      sites_cart=sites_cart,
      u_cart=u_cart,
      u_iso=u_iso,
      use_u_aniso=use_u_aniso,
      adp_similarity_proxies=adp_similarity_proxies,
      rigid_bond_proxies=rigid_bond_proxies,
      isotropic_adp_proxies=isotropic_adp_proxies,
      compute_gradients=compute_gradients,
      gradients_aniso_cart=gradients_aniso_cart,
      gradients_iso=gradients_iso,
      normalization=normalization)

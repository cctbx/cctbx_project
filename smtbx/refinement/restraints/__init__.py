import scitbx.sparse # import dependency

import boost.python
boost.python.import_ext("smtbx_refinement_restraints_ext")

from smtbx_refinement_restraints_ext import *

from cctbx.xray import parameter_map
from cctbx.adp_restraints import adp_restraint_params

from libtbx import adopt_optional_init_args

import sys

class manager(object):
  bond_proxies=None
  angle_proxies=None
  dihedral_proxies=None
  chirality_proxies=None
  planarity_proxies=None
  bond_similarity_proxies=None
  adp_similarity_proxies=None
  rigid_bond_proxies=None
  isotropic_adp_proxies=None
  #new
  fixed_u_eq_adp_proxies=None
  adp_u_eq_similarity_proxies=None
  adp_volume_similarity_proxies=None

  def __init__(self, **kwds):
    adopt_optional_init_args(self, kwds)

  def show_sorted(self, xray_structure,
                  f=None,
                  prefix="",
                  max_items=None):
    unit_cell = xray_structure.unit_cell()
    sites_cart = xray_structure.sites_cart()
    u_cart = xray_structure.scatterers().extract_u_cart(unit_cell)
    u_iso = xray_structure.scatterers().extract_u_iso()
    use_u_aniso = xray_structure.use_u_aniso()
    site_labels = xray_structure.scatterers().extract_labels()
    if (f is None): f = sys.stdout
    if (self.bond_proxies is not None):
      self.bond_proxies.show_sorted(
        by_value="residual", unit_cell=unit_cell,
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
    if (self.adp_u_eq_similarity_proxies is not None):
      self.adp_u_eq_similarity_proxies.show_sorted(
        by_value="residual", site_labels=site_labels,
        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.adp_volume_similarity_proxies is not None):
      self.adp_volume_similarity_proxies.show_sorted(
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
        by_value="residual", site_labels=site_labels,
        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso,
        f=f, prefix=prefix, max_items=max_items)
      print >> f
    if (self.fixed_u_eq_adp_proxies is not None):
      self.fixed_u_eq_adp_proxies.show_sorted(
        by_value="residual", site_labels=site_labels,
        u_cart=u_cart, u_iso=u_iso, use_u_aniso=use_u_aniso,
        f=f, prefix=prefix, max_items=max_items)
      print >> f

  def add_to_cif_block(self, cif_block, xray_structure):
    import iotbx.cif.restraints
    iotbx.cif.restraints.add_to_cif_block(
      cif_block, xray_structure,
      bond_proxies=self.bond_proxies,
      angle_proxies=self.angle_proxies,
      dihedral_proxies=self.dihedral_proxies,
      bond_similarity_proxies=self.bond_similarity_proxies,
      rigid_bond_proxies=self.rigid_bond_proxies,
      adp_similarity_proxies=self.adp_similarity_proxies,
      isotropic_adp_proxies=self.isotropic_adp_proxies)

  def build_linearised_eqns(self, xray_structure, parameter_map):
    n_restraints = 0
    n_params = parameter_map.n_parameters
    geometry_proxies = [proxies for proxies in (
      self.bond_proxies, self.angle_proxies, self.dihedral_proxies)
                        if proxies is not None]
    # count restraints, i.e. number of rows for restraint matrix
    n_restraints = sum([proxies.size() for proxies in geometry_proxies])
    if self.bond_similarity_proxies is not None:
      for proxy in self.bond_similarity_proxies:
        n_restraints += proxy.i_seqs.size()
      geometry_proxies.append(self.bond_similarity_proxies)
    adp_proxies = []
    if self.adp_similarity_proxies is not None:
      adp_proxies.append(self.adp_similarity_proxies)
      n_restraints += 6 * self.adp_similarity_proxies.size()
    if self.adp_u_eq_similarity_proxies is not None:
      adp_proxies.append(self.adp_u_eq_similarity_proxies)
      n_restraints += self.adp_u_eq_similarity_proxies.size()
    if self.adp_volume_similarity_proxies is not None:
      adp_proxies.append(self.adp_volume_similarity_proxies)
      n_restraints += self.adp_volume_similarity_proxies.size()
    if self.isotropic_adp_proxies is not None:
      adp_proxies.append(self.isotropic_adp_proxies)
      n_restraints += 6 * self.isotropic_adp_proxies.size()
    if self.fixed_u_eq_adp_proxies is not None:
      adp_proxies.append(self.fixed_u_eq_adp_proxies)
      n_restraints += self.fixed_u_eq_adp_proxies.size()
    if self.rigid_bond_proxies is not None:
      adp_proxies.append(self.rigid_bond_proxies)
      n_restraints += self.rigid_bond_proxies.size()
    # construct restraints matrix
    linearised_eqns = linearised_eqns_of_restraint(
      n_restraints, n_params)
    for proxies in geometry_proxies:
      linearise_restraints(
        xray_structure.unit_cell(), xray_structure.sites_cart(),
        parameter_map, proxies, linearised_eqns)
    u_cart = xray_structure.scatterers().extract_u_cart(
      xray_structure.unit_cell())
    params = adp_restraint_params(
      sites_cart=xray_structure.sites_cart(),
      u_cart=u_cart,
      u_iso=xray_structure.scatterers().extract_u_iso(),
      use_u_aniso=xray_structure.use_u_aniso())
    for proxies in adp_proxies:
      linearise_restraints(
        xray_structure.unit_cell(), params,
        parameter_map, proxies, linearised_eqns)
    return linearised_eqns

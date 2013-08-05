from __future__ import division

from libtbx import adopt_init_args

# catch-all class for handling any higher-level restraints (such as
# Ramachandran, rotamer, H-bonds, etc.)

class manager (object) :
  def __init__ (self,
                ramachandran_proxies=None,
                ramachandran_lookup=None,
                hydrogen_bond_proxies=None,
                hydrogen_bond_params=None,
                reference_manager=None,
                den_manager=None,
                ncs_manager=None,
                c_beta_dihedral_proxies=None,
                flags=None) :
    adopt_init_args(self, locals())
    if self.flags is None:
      import mmtbx.geometry_restraints.flags
      self.flags = mmtbx.geometry_restraints.flags.flags(default=True)
    assert (ramachandran_proxies is None) or (ramachandran_lookup is not None)
    if (self.hydrogen_bond_params is None) :
      from mmtbx.geometry_restraints import hbond
      self.hydrogen_bond_params = hbond.master_phil.fetch().extract()
    if (self.reference_manager is None) :
      from mmtbx.geometry_restraints import reference
      self.reference_manager = reference.manager()

  def get_n_proxies (self) :
    n_proxies = 0
    if (self.ramachandran_proxies is not None) :
      n_proxies += len(self.ramachandran_proxies)
    if (self.hydrogen_bond_proxies is not None) :
      if isinstance(self.hydrogen_bond_proxies, list) :
        n_proxies += len(self.hydrogen_bond_proxies)
      else :
        n_proxies += self.hydrogen_bond_proxies.size()
    if (self.reference_manager is not None):
      if (self.reference_manager.reference_coordinate_proxies is not None):
        n_proxies += \
          len(self.reference_manager.reference_coordinate_proxies)
      if (self.reference_manager.reference_torsion_proxies is not None):
        n_proxies += \
          len(self.reference_manager.reference_torsion_proxies)
    if (self.den_manager is not None) :
      n_proxies += len(self.den_manager.den_proxies)
    if (self.c_beta_dihedral_proxies is not None) :
      n_proxies += len(self.c_beta_dihedral_proxies)
    return n_proxies

  def get_n_hbonds (self) :
    if (self.hydrogen_bond_proxies is not None) :
      return len(self.hydrogen_bond_proxies)
    return 0

  def restraints_residual_sum (self,
                               sites_cart,
                               gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    target = 0
    if (self.ramachandran_proxies is not None and
        self.flags.ramachandran) :
      target += self.ramachandran_lookup.restraints_residual_sum(
        sites_cart=sites_cart,
        proxies=self.ramachandran_proxies,
        gradient_array=gradient_array)
    if (self.hydrogen_bond_proxies is not None and
        self.flags.hydrogen_bond) :
      from mmtbx.geometry_restraints import hbond
      lj_potential = self.hydrogen_bond_params.lennard_jones.potential
      target += hbond.target_and_gradients(
        proxies=self.hydrogen_bond_proxies,
        sites_cart=sites_cart,
        gradient_array=gradient_array,
        falloff_distance=self.hydrogen_bond_params.falloff_distance,
        lennard_jones_potential=lj_potential)
    if (self.reference_manager is not None and
        self.flags.reference) :
      if (self.reference_manager.reference_coordinate_proxies is not None or
          self.reference_manager.reference_torsion_proxies is not None):
        target += self.reference_manager.target_and_gradients(
          sites_cart=sites_cart,
          gradient_array=gradient_array)
    if (self.den_manager is not None and
        self.flags.den) :
      #print "DEN target is in geneneric manager"
      den_target = self.den_manager.target_and_gradients(
        sites_cart=sites_cart,
        gradient_array=gradient_array)
      #print "DEN target: %.1f" % den_target
      target += den_target
    if (self.c_beta_dihedral_proxies is not None and
        self.flags.c_beta) :
      from mmtbx.geometry_restraints import c_beta
      c_beta_target = c_beta.target_and_gradients(
        sites_cart=sites_cart,
        c_beta_dihedral_proxies=self.c_beta_dihedral_proxies,
        gradient_array=gradient_array)
      target += c_beta_target
    return target

  def hbonds_as_simple_bonds (self) :
    if (self.hydrogen_bond_proxies is not None) :
      from mmtbx.geometry_restraints import hbond
      return hbond.get_simple_bonds(self.hydrogen_bond_proxies)
    return []

  def rotamers (self) :
    return None #self.rotamer_manager

  def update_hydrogen_bonds (self,
                             pdb_hierarchy,
                             xray_structure,
                             params=None,
                             log=None) :
    from mmtbx.geometry_restraints import hbond
    if (params is None) :
      params = self.hydrogen_bond_params
    self.hydrogen_bond_proxies = hbond.find_implicit_hydrogen_bonds(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=xray_structure,
      params=self.hydrogen_bond_params,
      log=log).proxies

  def select (self,
              n_seq,
              iselection) :
    ramachandran_proxies = hydrogen_bond_proxies = den_manager = None
    c_beta_dihedral_proxies = ncs_manager = None
    if (self.ramachandran_proxies is not None) :
      ramachandran_proxies = self.ramachandran_proxies.proxy_select(
        n_seq, iselection)
    if (self.hydrogen_bond_proxies is not None) :
      hydrogen_bond_proxies = self.hydrogen_bond_proxies.proxy_select(
        n_seq, iselection)
    if (self.den_manager is not None) :
      den_manager = self.den_manager.select(n_seq, iselection)
    if (self.c_beta_dihedral_proxies is not None) :
      c_beta_dihedral_proxies = self.c_beta_dihedral_proxies.proxy_select(
        n_seq, iselection)
    if (self.ncs_manager is not None) :
      ncs_manager = self.ncs_manager.select(n_seq, iselection)
    return manager(
      ramachandran_proxies=ramachandran_proxies,
      ramachandran_lookup=self.ramachandran_lookup,
      hydrogen_bond_proxies=hydrogen_bond_proxies,
      hydrogen_bond_params=self.hydrogen_bond_params,
      den_manager=den_manager,
      ncs_manager=ncs_manager,
      c_beta_dihedral_proxies=c_beta_dihedral_proxies,
      flags=self.flags)

  def add_c_beta_torsion_restraints(self,
                                    pdb_hierarchy,
                                    selection=None,
                                    sigma=2.5):
    from mmtbx.geometry_restraints import c_beta
    self.c_beta_dihedral_proxies = \
      c_beta.get_c_beta_torsion_proxies(
        pdb_hierarchy=pdb_hierarchy,
        selection=selection,
        sigma=2.5)

  def remove_c_beta_torsion_restraints(self, selection):
    self.c_beta_dihedral_proxies = \
      self.c_beta_dihedral_proxies.proxy_remove(selection=selection)

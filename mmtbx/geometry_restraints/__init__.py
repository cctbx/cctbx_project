
from libtbx import adopt_init_args

# XXX catch-all class for handling any higher-level restraints (such as
# Ramachandran, rotamer, H-bonds, etc.)

class manager (object) :
  def __init__ (self,
                ramachandran_proxies=None,
                ramachandran_lookup=None,
                hydrogen_bond_proxies=None,
                hydrogen_bond_params=None,
                rotamer_manager=None,
                den_manager=None,
                flags=None) :
    adopt_init_args(self, locals())
    if self.flags is None:
      import mmtbx.geometry_restraints.flags
      self.flags = mmtbx.geometry_restraints.flags.flags(default=True)
    assert (ramachandran_proxies is None) or (ramachandran_lookup is not None)
    if (self.hydrogen_bond_params is None) :
      from mmtbx.geometry_restraints import hbond
      self.hydrogen_bond_params = hbond.master_phil.fetch().extract()

  def get_n_proxies (self) :
    n_proxies = 0
    if (self.ramachandran_proxies is not None) :
      n_proxies += len(self.ramachandran_proxies)
    if (self.hydrogen_bond_proxies is not None) :
      if isinstance(self.hydrogen_bond_proxies, list) :
        n_proxies += len(self.hydrogen_bond_proxies)
      else :
        n_proxies += self.hydrogen_bond_proxies.size()
    if (self.den_manager is not None) :
      n_proxies += len(self.den_manager.den_proxies)
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
    if (self.den_manager is not None and
        self.flags.den) :
      #print "DEN target is in geneneric manager"
      den_target = self.den_manager.target_and_gradients(
        sites_cart=sites_cart,
        gradient_array=gradient_array)
      #print "DEN target: %.1f" % den_target
      target += den_target
    return target

  def hbonds_as_simple_bonds (self) :
    if (self.hydrogen_bond_proxies is not None) :
      from mmtbx.geometry_restraints import hbond
      return hbond.get_simple_bonds(self.hydrogen_bond_proxies)
    return []

  def rotamers (self) :
    return self.rotamer_manager

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

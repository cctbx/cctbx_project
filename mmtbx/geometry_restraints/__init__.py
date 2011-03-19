
from libtbx import adopt_init_args

# XXX catch-all class for handling any higher-level restraints (such as
# Ramachandran, rotamer, H-bonds, etc.)

class manager (object) :
  def __init__ (self,
                ramachandran_proxies=None,
                ramachandran_lookup=None,
                hydrogen_bond_proxies=None,
                hydrogen_bond_params=None) :
    adopt_init_args(self, locals())
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
    return n_proxies

  def restraints_residual_sum (self,
                               sites_cart,
                               gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    target = 0
    if (self.ramachandran_proxies is not None) :
      target += self.ramachandran_lookup.restraints_residual_sum(
        sites_cart=sites_cart,
        proxies=self.ramachandran_proxies,
        gradient_array=gradient_array)
    if (self.hydrogen_bond_proxies is not None) :
      from mmtbx.geometry_restraints import hbond
      lj_potential = self.hydrogen_bond_params.lennard_jones.potential
      target += hbond.target_and_gradients(
        proxies=self.hydrogen_bond_proxies,
        sites_cart=sites_cart,
        gradient_array=gradient_array,
        falloff_distance=self.hydrogen_bond_params.falloff_distance,
        lennard_jones_potential=lj_potential)
    return target

  def hbonds_as_simple_bonds (self) :
    if (self.hydrogen_bond_proxies is not None) :
      from mmtbx.geometry_restraints import hbond
      return hbond.get_simple_bonds(self.hydrogen_bond_proxies)
    return []

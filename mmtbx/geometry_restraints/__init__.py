from __future__ import division

from libtbx import adopt_init_args
from scitbx.array_family import flex


# catch-all class for handling any higher-level restraints (such as
# Ramachandran, rotamer, H-bonds, etc.)

class manager (object) :
  def __init__ (self,
                reference_manager=None,
                flags=None) :
    adopt_init_args(self, locals())
    if self.flags is None:
      import mmtbx.geometry_restraints.flags
      self.flags = mmtbx.geometry_restraints.flags.flags(default=True)
    if (self.reference_manager is None) :
      from mmtbx.geometry_restraints import reference
      self.reference_manager = reference.manager()

  def get_n_proxies(self):
    return self.get_n_reference_torsion_proxies()

  def get_n_reference_torsion_proxies(self):
    if self.reference_manager is not None:
      if self.reference_manager.reference_torsion_proxies is not None:
        return len(self.reference_manager.reference_torsion_proxies)
    return 0

  def restraints_residual_sum (self,
                               sites_cart,
                               gradient_array=None) :
    if (gradient_array is None) :
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    target = 0
    if (self.reference_manager is not None and
        self.flags.reference) :
      if (self.reference_manager.reference_torsion_proxies is not None):
        target += self.reference_manager.target_and_gradients(
          sites_cart=sites_cart,
          gradient_array=gradient_array)
    return target

  def rotamers (self) :
    return None #self.rotamer_manager

  def select (self,
              n_seq,
              iselection) :
    if self.reference_manager is not None:
      reference_manager = self.reference_manager.select(n_seq, iselection)
    return manager(
      reference_manager=reference_manager,
      flags=self.flags)

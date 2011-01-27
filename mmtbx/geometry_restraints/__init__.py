
from libtbx import adopt_init_args

# XXX catch-all class for handling any higher-level restraints (such as
# Ramachandran, rotamer, H-bonds, etc.)

class manager (object) :
  def __init__ (self,
                ramachandran_proxies=None,
                ramachandran_lookup=None) :
    adopt_init_args(self, locals())
    assert (ramachandran_proxies is None) or (ramachandran_lookup is not None)

  def get_proxies (self) :
    proxies = []
    if (self.ramachandran_proxies is not None) :
      proxies.extend(self.ramachandran_proxies)
    return proxies

  def restraints_residual_sum (self,
                               proxies, # XXX ignored, but passed anyway
                               sites_cart,
                               gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    target = 0
    if (self.ramachandran_proxies is not None) :
      target += self.ramachandran_lookup.restraints_residual_sum(
        sites_cart=sites_cart,
        proxies=proxies,
        gradient_array=gradient_array)
    return target


class manager (object) :
  def __init__ (self) :
    self._proxy_sets = []
    self._handlers = []

  def add_proxies (self, proxies, handler_method) :
    assert hasattr(handler_method, "__call__")
    self._proxy_sets.append(proxies)
    self._handlers.append(handler_method)

  def target_and_gradients (self,
                            sites_cart,
                            gradient_array=None) :
    if (gradient_array is None) :
      from scitbx.array_family import flex
      gradient_array = flex.vec3_double(sites_cart.size(), (0.0,0.0,0.0))
    sum = 0
    for (proxies, handler) in zip(self._proxy_sets, self._handlers) :
      sum += handler(proxies=proxies,
                     sites_cart=sites_cart,
                     gradient_array=gradient_array)
    return sum

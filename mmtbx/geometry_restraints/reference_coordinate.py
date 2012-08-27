from __future__ import division
from cctbx.array_family import flex

class build_proxies(object):
  def __init__(self,
               sites_cart,
               selection=None,
               sigma=0.5):
    import boost.python
    self.ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
    self.reference_coordinate_proxies = \
      self.ext.shared_reference_coordinate_proxy()
    if (selection is not None):
      if (isinstance(selection, flex.bool)):
        selection = selection.iselection()
    if selection is None:
      selection = flex.bool(
        len(sites_cart),
        True).iselection()
    assert len(sites_cart) == len(selection)
    weight = 1.0 / (sigma**2)
    for k, i_seq in enumerate(selection):
      i_seqs = [i_seq]
      ref_sites = sites_cart[k]
      proxy = self.ext.reference_coordinate_proxy(
                i_seqs=i_seqs,
                ref_sites=ref_sites,
                weight=weight)
      self.reference_coordinate_proxies.append(proxy)

def target_and_gradients(proxies,
                         sites_cart,
                         gradient_array):
  import boost.python
  ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
  return ext.reference_coordinate_residual_sum(
    sites_cart,
    proxies,
    gradient_array)

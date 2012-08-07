from cctbx.array_family import flex

class build_proxies(object):
  def __init__(self,
               sites_cart,
               pdb_hierarchy,
               selection=None,
               sites_cart_reference=None,
               c_alpha_only=False,
               sigma=0.5):
    import boost.python
    self.ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
    self.reference_coordinate_proxies = \
      self.ext.shared_reference_coordinate_proxy()
    if c_alpha_only:
      ca_iselection = self.get_c_alpha_iselection(pdb_hierarchy)
      selection = flex.bool(
        len(sites_cart),
        ca_iselection)
    elif selection is None:
      selection = flex.bool(
        len(sites_cart),
        True)
    if sites_cart_reference is None:
      sites_cart_reference = sites_cart
    assert len(sites_cart) == len(sites_cart_reference)
    weight = 1.0 / (sigma**2)
    for i_seq, sel in enumerate(selection):
      if sel:
        i_seqs = tuple([i_seq])
        ref_sites = sites_cart_reference[i_seq]
        proxy = self.ext.reference_coordinate_proxy(
                  i_seqs=i_seqs,
                  ref_sites=ref_sites,
                  weight=weight)
        self.reference_coordinate_proxies.append(proxy)

  def get_c_alpha_iselection(self, pdb_hierarchy):
    ca_iselection = flex.size_t()
    for atom in pdb_hierarchy.atoms():
      if atom.name == " CA ":
        ca_iselection.append(atom.i_seq)
    return ca_iselection

def target_and_gradients(proxies,
                         sites_cart,
                         gradient_array):
  import boost.python
  ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
  return ext.reference_coordinate_residual_sum(
    sites_cart,
    proxies,
    gradient_array)

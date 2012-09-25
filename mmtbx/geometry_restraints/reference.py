
from __future__ import division
from cctbx.array_family import flex
from cctbx import geometry_restraints
import boost.python
ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")

class manager(object):
  def __init__(self):
    self.reference_coordinate_proxies = None
    self.reference_torsion_proxies = None

  def add_coordinate_restraints(
        self,
        sites_cart,
        selection=None,
        sigma=0.5):
    self.add_reference_coordinate_proxies(
      sites_cart=sites_cart,
      selection=selection,
      sigma=sigma)

  def add_torsion_restraints(
        self,
        pdb_hierarchy,
        sites_cart,
        selection=None,
        sigma=2.5):
    self.add_reference_torsion_proxies(
      pdb_hierarchy=pdb_hierarchy,
      sites_cart=sites_cart,
      selection=selection,
      sigma=sigma)

  def add_reference_coordinate_proxies(
        self,
        sites_cart,
        selection=None,
        sigma=0.5):
    import boost.python
    self.ext = boost.python.import_ext("mmtbx_reference_coordinate_ext")
    if self.reference_coordinate_proxies is None:
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

  def add_reference_torsion_proxies(
        self,
        pdb_hierarchy,
        sites_cart,
        selection=None,
        sigma=2.5):
    from mmtbx.torsion_restraints.reference_model import build_torsion_proxies
    local_reference_torsion_proxies = \
        build_torsion_proxies(
          pdb_hierarchy=pdb_hierarchy,
          sites_cart=sites_cart,
          selection=selection,
          sigma=sigma)
    if self.reference_torsion_proxies is None:
      self.reference_torsion_proxies = local_reference_torsion_proxies
    else:
      for dp in local_reference_torsion_proxies:
        self.reference_torsion_proxies.append(dp)

  def remove_coordinate_restraints(self, selection):
    self.reference_coordinate_proxies = \
      self.reference_coordinate_proxies.proxy_remove(selection=selection)

  def remove_torsion_restraints(self, selection):
    self.reference_torsion_proxies = \
      self.reference_torsion_proxies.proxy_remove(selection=selection)

  def target_and_gradients(self,
                           sites_cart,
                           gradient_array):
    target = 0.0
    if self.reference_coordinate_proxies is not None:
      target += ext.reference_coordinate_residual_sum(
        sites_cart,
        self.reference_coordinate_proxies,
        gradient_array)
    if self.reference_torsion_proxies is not None:
      target += geometry_restraints.dihedral_residual_sum(
                  sites_cart=sites_cart,
                  proxies=self.reference_torsion_proxies,
                  gradient_array=gradient_array)
    return target

from __future__ import division

from cctbx import geometry_restraints
from scitbx import matrix
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
#from time import time

class manager(object):
  def __init__(self,
      hierarchy,
      geometry_restraints,
      crystal_symmetry     = None,
      xray_structure       = None,
      idealize             = True):
    cs = crystal_symmetry
    # maybe not necessary to allow reading in xrs - think about it later
    assert cs is not None or xray_structure is not None
    if xray_structure is not None:
      assert hierarchy.atoms_size() == xray_structure.scatterers().size()
    if xray_structure is None:
      self.xray_structure = hierarchy.extract_xray_structure(
        crystal_symmetry=cs)
      self.cs = cs
    else:
      self.xray_structure = xray_structure
      self.cs = xray_structure.crystal_symmetry()
    self.hierarchy = hierarchy
    sites_cart = hierarchy.atoms().extract_xyz()

    self.connectivity_obj = connectivity.determine_H_neighbors(
      geometry_restraints           = geometry_restraints,
      pdb_hierarchy       = self.hierarchy)
    self.h_connectivity = self.connectivity_obj.h_connectivity
    self.h_parameterization = parameterization.get_h_parameterization(
      h_connectivity   = self.h_connectivity,
      sites_cart     = sites_cart,
      idealize       = idealize)

  def idealize_hydrogens(self, hierarchy, xray_structure=None):
    """ Doing idealization in place, maybe it is better to return a copy, but
    not sure. """
    # some safeguarding is necessary, like
    # I'm sure you'll come up with something more to make sure this hierarchy
    # is consistent with what was used to initialize this object.
    assert hierarchy.atoms_size() == self.hierarchy.atoms_size()
    if xray_structure is None:
      xray_structure = hierarchy.extract_xray_structure(
        crystal_symmetry=self.cs)
    sites_cart = hierarchy.atoms().extract_xyz()

    for ih in self.h_parameterization.keys():
      hp = self.h_parameterization[ih]
      rh = matrix.col(sites_cart[ih])
      rh_calc = None
      if (hp.htype in ['unk','unk_ideal']):
        rh_calc = rh
      else:
        rh_calc = parameterization.compute_H_position(
          ih         = ih,
          sites_cart = sites_cart,
          hp         = hp)
        if (rh_calc is None):
          rh_calc = rh
      sites_cart[ih] = rh_calc
      #if h_obj.rH_gen is not None:
        # XXX This is bug guard. For some reason in my model I'm getting None
        # somewhere in the middle --> DL: should work now but keep in mind
    xray_structure.set_sites_cart(sites_cart)
    hierarchy.adopt_xray_structure(xray_structure)


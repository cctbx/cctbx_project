from __future__ import division
from cctbx import geometry_restraints
from scitbx import matrix
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
#from time import time

class manager(object):
  def __init__(self,
      pdb_hierarchy,
      geometry_restraints,
      xray_structure         = None,
      use_ideal_bonds_angles = True):
    # maybe not necessary to allow reading in xrs - think about it later
    # TODO: more safeguarding
#    assert cs is not None or xray_structure is not None
    self.cs = geometry_restraints.crystal_symmetry
    if xray_structure is not None:
      assert pdb_hierarchy.atoms_size() == xray_structure.scatterers().size()
    if xray_structure is None:
      self.xray_structure = pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.cs)
    else:
      self.xray_structure = xray_structure
    self.pdb_hierarchy = pdb_hierarchy
    sites_cart = pdb_hierarchy.atoms().extract_xyz()

    self.connectivity_obj = connectivity.determine_H_neighbors(
      geometry_restraints = geometry_restraints,
      pdb_hierarchy       = self.pdb_hierarchy)
    self.h_connectivity = self.connectivity_obj.h_connectivity
    self.h_parameterization = parameterization.get_h_parameterization(
      h_connectivity         = self.h_connectivity,
      sites_cart             = sites_cart,
      use_ideal_bonds_angles = use_ideal_bonds_angles)

  def idealize_hydrogens_inplace(self, pdb_hierarchy=None, xray_structure=None):
    """ Doing idealization in place, maybe it is better to return a copy, but
    not sure. """
    # some safeguarding is necessary, like
    # I'm sure you'll come up with something more to make sure this hierarchy
    # is consistent with what was used to initialize this object.
    if pdb_hierarchy is not None:
      assert pdb_hierarchy.atoms_size() == self.pdb_hierarchy.atoms_size()
    if xray_structure is None:
      xray_structure = pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.cs)
    #sites_cart = pdb_hierarchy.atoms().extract_xyz()
    sites_cart = xray_structure.sites_cart()

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
    if pdb_hierarchy is not None:
      pdb_hierarchy.adopt_xray_structure(xray_structure)

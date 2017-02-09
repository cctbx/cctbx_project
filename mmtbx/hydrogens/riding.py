from __future__ import division
from cctbx import geometry_restraints
from scitbx import matrix
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
from mmtbx.hydrogens import modify_gradients
#from mmtbx_hydrogens_ext import *

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

    self.connectivity_manager = connectivity.determine_connectivity(
      pdb_hierarchy       = self.pdb_hierarchy,
      geometry_restraints = geometry_restraints)
    self.h_connectivity = self.connectivity_manager.h_connectivity
    self.double_H = self.connectivity_manager.double_H
    self.parameterization_manager = parameterization.manager(
      h_connectivity         = self.h_connectivity,
      sites_cart             = sites_cart,
      use_ideal_bonds_angles = use_ideal_bonds_angles)
    self.h_parameterization = \
      self.parameterization_manager.determine_parameterization()

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

    for hp in self.h_parameterization:
      if (hp == None): continue
      ih = hp.ih
      #hp = self.h_parameterization[ih]
      rh = matrix.col(sites_cart[ih])
      rh_calc = None
      if (hp.htype in ['unk','unk_ideal']):
        rh_calc = rh
      else:
        rh_calc = parameterization.compute_H_position(
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

  #def idealize_hydrogens_inplace_cpp(self, pdb_hierarchy=None, xray_structure=None):
  #  if pdb_hierarchy is not None:
  #    assert pdb_hierarchy.atoms_size() == self.pdb_hierarchy.atoms_size()
  #  if xray_structure is None:
  #    xray_structure = pdb_hierarchy.extract_xray_structure(
  #      crystal_symmetry=self.cs)
  #  sites_cart = xray_structure.sites_cart()
#
  #  for ih in self.h_parameterization:
  #    hp = self.h_parameterization[ih]
  #    rc = riding_coefficients(
  #      htype=hp.htype, a0=hp.a0, a1=hp.a1, a2=hp.a2, a=hp.a, b=hp.b,
  #      h=hp.h, n=hp.n, disth=hp.dist_h)
  #    rh = matrix.col(sites_cart[ih])
  #    rh_calc = None
  #    if (hp.htype in ['unk','unk_ideal']):
  #      rh_calc = rh
  #    else:
  #      rh_calc = compute_H_position(
  #        riding_coefficients = rc,
  #        sites_cart          = sites_cart,
  #        ih                  = ih)
  #      if (rh_calc is None):
  #        rh_calc = rh
  #    sites_cart[ih] = rh_calc
  #  xray_structure.set_sites_cart(sites_cart)
  #  if pdb_hierarchy is not None:
  #    pdb_hierarchy.adopt_xray_structure(xray_structure)

  def gradients_reduced(self, grads, sites_cart, hd_selection):
    modify_gradients.modify_gradients(
      sites_cart          = sites_cart,
      h_parameterization  = self.h_parameterization,
      grads               = grads)
    grads = grads.select(~hd_selection)
    return grads

  def diagnostics(self, sites_cart, threshold):
    diagnostics = parameterization.diagnostics(
      h_parameterization = self.h_parameterization,
      h_connectivity     = self.h_connectivity,
      sites_cart         = sites_cart,
      threshold          = threshold)
    return diagnostics


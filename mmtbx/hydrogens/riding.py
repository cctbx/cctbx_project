from __future__ import division
import time
from libtbx import group_args
from cctbx import geometry_restraints
from scitbx import matrix
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
from mmtbx.hydrogens import modify_gradients
from mmtbx_hydrogens_ext import *

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

    self.parameterization_cpp = []
#    t0=time.time()
    for hp in self.h_parameterization:
      if (hp is not None):
        self.parameterization_cpp.append(hp)
#    print "initialize cpp para", time.time()-t0


  def idealize_hydrogens_inplace(
      self,
      pdb_hierarchy  = None,
      xray_structure = None):
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

  def idealize_hydrogens_inplace_cpp(self, pdb_hierarchy=None, xray_structure=None):
    if pdb_hierarchy is not None:
      assert pdb_hierarchy.atoms_size() == self.pdb_hierarchy.atoms_size()
    if xray_structure is None:
      xray_structure = pdb_hierarchy.extract_xray_structure(
        crystal_symmetry=self.cs)
    sites_cart = xray_structure.sites_cart()
    sites_cart_new = apply_new_H_positions(
      sites_cart = sites_cart,
      parameterization = self.parameterization_cpp)
    xray_structure.set_sites_cart(sites_cart_new)
    if pdb_hierarchy is not None:
      pdb_hierarchy.adopt_xray_structure(xray_structure)

  def gradients_reduced(self, grads, sites_cart, hd_selection):
    modify_gradients.modify_gradients(
      sites_cart          = sites_cart,
      h_parameterization  = self.parameterization_cpp,
      grads               = grads)
    grads = grads.select(~hd_selection)
    return grads

  def gradients_reduced_cpp(self, gradients, sites_cart, hd_selection):
    new_gradients = modify_gradients_cpp(
      gradients         = gradients,
      sites_cart        = sites_cart,
      parameterization  = self.parameterization_cpp)
    new_gradients = new_gradients.select(~hd_selection)
    return new_gradients

  def diagnostics(self, sites_cart, threshold):
    h_parameterization = self.h_parameterization
    h_connectivity     = self.h_connectivity
    h_distances = {}
    long_distance_list, list_h, type_list = [], [], []
    for hp in h_parameterization:
      if (hp == None): continue
      ih = hp.ih
      list_h.append(ih)
      rh = matrix.col(sites_cart[ih])
      rh_calc = parameterization.compute_H_position(
        sites_cart = sites_cart,
        hp         = hp)
      if (rh_calc is not None):
        h_distance = (rh_calc - rh).length()
        h_distances[ih] = h_distance
        type_list.append(hp.htype)
      if (h_distance > threshold):
        long_distance_list.append(ih)
    #set_temp = set(list_h)
    #set_temp = set(list(h_parameterization.keys()))
    #slipped = [x for x in h_connectivity if x not in set_temp]
    return group_args(
      h_distances        = h_distances,
      unk_list           = self.parameterization_manager.unk_list,
      unk_ideal_list     = self.parameterization_manager.unk_ideal_list,
      long_distance_list = long_distance_list,
#      slipped            = slipped,
      type_list          = type_list,
      number_h_para      = len(list_h),
      threshold          = threshold)

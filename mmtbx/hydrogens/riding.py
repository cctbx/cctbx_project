from __future__ import division
from libtbx import group_args
from cctbx import geometry_restraints
from scitbx import matrix
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
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
    self.connectivity_slipped = self.connectivity_manager.connectivity_slipped
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

  def idealize_hydrogens_inplace_cpp(
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
    sites_cart = xray_structure.sites_cart()
    sites_cart_new = apply_new_H_positions(
      sites_cart = sites_cart,
      parameterization = self.parameterization_cpp)
      # XXX Think about bug guard if calculated position is None
    xray_structure.set_sites_cart(sites_cart_new)
    if pdb_hierarchy is not None:
      pdb_hierarchy.adopt_xray_structure(xray_structure)

  def gradients_reduced_cpp(self, gradients, sites_cart, hd_selection):
    new_gradients = modify_gradients_cpp(
      gradients         = gradients,
      sites_cart        = sites_cart,
      parameterization  = self.parameterization_cpp)
    new_gradients = new_gradients.select(~hd_selection)
    return new_gradients

  def diagnostics(self, sites_cart, threshold):
    h_distances = {}
    unk_list = self.parameterization_manager.unk_list
    unk_ideal_list = self.parameterization_manager.unk_ideal_list
    long_distance_list, list_h, type_list = [], [], []
    for rc in self.parameterization_cpp:
      ih = rc.ih
      list_h.append(ih)
      rh = matrix.col(sites_cart[ih])
      rh_calc = matrix.col(compute_h_position(
        riding_coefficients = rc,
        sites_cart          = sites_cart))
      # add safeguard? What if no coordinates?
      h_distance = (rh_calc - rh).length()
      h_distances[ih] = h_distance
      type_list.append(rc.htype)
      if (h_distance > threshold):
        long_distance_list.append(ih)
    all_h_in_para = list_h + unk_list + unk_ideal_list
    list_h_connect = []
    for item in self.h_connectivity:
      if (item):
        list_h_connect.append(item.ih)
    #set_temp = set(list(h_parameterization.keys()))
    slipped = [x for x in list_h_connect if x not in set(all_h_in_para)]
    return group_args(
      h_distances          = h_distances,
      unk_list             = unk_list,
      unk_ideal_list       = unk_ideal_list,
      long_distance_list   = long_distance_list,
      slipped              = slipped,
      type_list            = type_list,
      number_h_para        = len(list_h),
      n_connect            = len(list_h_connect),
      threshold            = threshold,
      double_H             = self.double_H,
      connectivity_slipped = self.connectivity_slipped)

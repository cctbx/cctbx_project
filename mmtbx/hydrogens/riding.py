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
      use_ideal_bonds_angles = True):
    self.pdb_hierarchy = pdb_hierarchy
    connectivity_manager = connectivity.determine_connectivity(
      pdb_hierarchy       = self.pdb_hierarchy,
      geometry_restraints = geometry_restraints)
    h_connectivity = connectivity_manager.h_connectivity

    diagnostics_connectivity = connectivity_manager.get_diagnostics()
    self.h_in_connectivity = diagnostics_connectivity.h_in_connectivity
    self.double_H = diagnostics_connectivity.double_H
    self.connectivity_slipped = diagnostics_connectivity.connectivity_slipped

    self.parameterization_manager = parameterization.manager(
      h_connectivity         = h_connectivity,
      sites_cart             = self.pdb_hierarchy.atoms().extract_xyz(),
      use_ideal_bonds_angles = use_ideal_bonds_angles)
    self.h_parameterization = self.parameterization_manager.h_parameterization
    self.parameterization_cpp = []
    for hp in self.h_parameterization:
      if (hp is not None):
        self.parameterization_cpp.append(hp)
    self.hd_selection = self.pdb_hierarchy.atom_selection_cache().\
      selection("element H or element D")
    self.not_hd_selection = ~self.hd_selection

  def refinable_parameters_init(self):
    return flex.double(self.n_parameters(), 0)

  def n_parameters(self):
    return self.not_hd_selection.count(True)*3

  def idealize(self, sites_cart=None, pdb_hierarchy=None, xray_structure=None):
    assert [sites_cart, pdb_hierarchy, xray_structure].count(None) in [2,3]
    if(xray_structure is not None):
      sites_cart = xray_structure.sites_cart()
      apply_new_H_positions(
        sites_cart = sites_cart,
        parameterization = self.parameterization_cpp)
      xray_structure.set_sites_cart(sites_cart)
    elif(pdb_hierarchy is not None):
      sites_cart = pdb_hierarchy.atoms().extract_xyz()
      apply_new_H_positions(
        sites_cart = sites_cart,
        parameterization = self.parameterization_cpp)
      pdb_hierarchy.atoms().set_xyz(sites_cart)
    elif(sites_cart is not None):
      apply_new_H_positions(
        sites_cart = sites_cart,
        parameterization = self.parameterization_cpp)
    else:
      sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
      apply_new_H_positions(
        sites_cart = sites_cart,
        parameterization = self.parameterization_cpp)
      self.pdb_hierarchy.atoms().set_xyz(sites_cart)

  def gradients_reduced_cpp(self, gradients, sites_cart, hd_selection):
    new_gradients = modify_gradients_cpp(
      gradients         = gradients,
      sites_cart        = sites_cart,
      parameterization  = self.parameterization_cpp)
    new_gradients = new_gradients.select(~hd_selection)
    return new_gradients

  def print_parameterization_info(self, log):
    unk_list = self.parameterization_manager.unk_list
    double_H = self.double_H
    atoms = self.pdb_hierarchy.atoms()
    if unk_list or double_H:
      number = len(unk_list) + len(double_H)
      m = """  The following atoms are not used in the riding H procedure. This typically
  happens when
  - A neighboring atom (H or non H) is missing.
  - Restraints involving the H atom are incomplete (this occurs most likely when
    custom restraints are supplied).
  - An H atom has less than 3 non-H covalently bound partners.
    (This will be addressed in the future)
  It is not worrysome if a few atoms are listed here; but it is always helpful
  to check the environment of these H atoms, as it might hint to missing atoms
  or restraints.\n"""
      print >> log, m
      print >> log, '  Number of H atoms not used in the riding H procedure: ', \
        number, '\n'
      for ih in unk_list:
        atom = atoms[ih]
        labels = atom.fetch_labels()
        if (labels.resname in ['DOD', 'WAT', 'HOH']):
          continue
        print >> log, '\t', atom.id_str()
      for ih in double_H.keys():
        atom = atoms[ih]
        print >> log, '\t', atom.id_str()

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
#    list_h_connect = []
#    for item in self.h_connectivity:
#      if (item):
#        list_h_connect.append(item.ih)
    slipped = [x for x in self.h_in_connectivity if x not in set(all_h_in_para)]
    return group_args(
      h_distances          = h_distances,
      unk_list             = unk_list,
      unk_ideal_list       = unk_ideal_list,
      long_distance_list   = long_distance_list,
      slipped              = slipped,
      type_list            = type_list,
      number_h_para        = len(list_h),
#      n_connect            = len(list_h_connect),
      threshold            = threshold,
      double_H             = self.double_H,
      connectivity_slipped = self.connectivity_slipped)

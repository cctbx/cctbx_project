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
    parameterization_manager = parameterization.manager(
      h_connectivity         = h_connectivity,
      sites_cart             = self.pdb_hierarchy.atoms().extract_xyz(),
      use_ideal_bonds_angles = use_ideal_bonds_angles)
    self.h_parameterization = parameterization_manager.h_parameterization
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
    list_h = []
    for rc in self.parameterization_cpp:
      ih = rc.ih
      list_h.append(ih)
    atoms = self.pdb_hierarchy.atoms()

    all_H_model = \
      list(self.pdb_hierarchy.select(self.hd_selection).atoms().extract_i_seq())
    unk_list = [x for x in all_H_model if x not in list_h]

    if unk_list:
      number = len(unk_list)
      m = """  The following atoms are not used in the riding H procedure. This typically
  happens when
  - A neighboring atom (H or non H) is missing.
  - Restraints involving the H atom are incomplete (this occurs most likely when
    custom restraints are supplied).
  - An H atom has less than 3 non-H covalently bound partners.
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

  def diagnostics(self, sites_cart, threshold):
    h_distances = {}
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

    return group_args(
      h_distances          = h_distances,
      long_distance_list   = long_distance_list,
      type_list            = type_list)


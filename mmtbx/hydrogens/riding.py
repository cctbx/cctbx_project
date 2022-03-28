from __future__ import absolute_import, division, print_function
from libtbx import group_args
from scitbx import matrix
from cctbx.array_family import flex
from scitbx_array_family_flex_ext import reindexing_array
from mmtbx.hydrogens import connectivity
from mmtbx.hydrogens import parameterization
from mmtbx_hydrogens_ext import *

class manager(object):
  def __init__(self,
      pdb_hierarchy,
      geometry_restraints,
      use_ideal_bonds_angles = True,
      process_manager        = True,
      use_ideal_dihedral = False,
      ignore_h_with_dof  = False):
    self.pdb_hierarchy = pdb_hierarchy
    self.geometry_restraints = geometry_restraints
    self.use_ideal_bonds_angles = use_ideal_bonds_angles
    self.ignore_h_with_dof = ignore_h_with_dof
    #
    self.hd_selection = self.pdb_hierarchy.atom_selection_cache().\
      selection("element H or element D")
    self.not_hd_selection = ~self.hd_selection
    self.h_parameterization = []
    self.parameterization_cpp = []
    if process_manager is True:
      connectivity_manager = connectivity.determine_connectivity(
        pdb_hierarchy       = self.pdb_hierarchy,
        geometry_restraints = geometry_restraints)
      h_connectivity = connectivity_manager.h_connectivity
      diagnostics_connectivity = connectivity_manager.get_diagnostics()
      parameterization_manager = parameterization.manager(
        h_connectivity         = h_connectivity,
        sites_cart             = self.pdb_hierarchy.atoms().extract_xyz(),
        use_ideal_bonds_angles = use_ideal_bonds_angles,
        site_labels         = [atom.id_str().replace('pdb=','').replace('"','')
                                  for atom in pdb_hierarchy.atoms()],
        use_ideal_dihedral = use_ideal_dihedral,
        ignore_h_with_dof  = self.ignore_h_with_dof)
      self.h_parameterization = parameterization_manager.h_parameterization
      self.parameterization_cpp = self.get_parameterization_cpp(
        h_parameterization = self.h_parameterization)

  def copy_h_parameterization(self, h_parameterization):
    new_h_parameterization = []
    for rc in h_parameterization:
      if rc is not None:
        rc_copy = riding_coefficients(rc)
        new_h_parameterization.append(rc_copy)
      else:
        new_h_parameterization.append(None)
    return new_h_parameterization

  def update(self, pdb_hierarchy, geometry_restraints, n_new_atoms):
    new_h_parameterization = self.copy_h_parameterization(
      h_parameterization = self.h_parameterization)
    new_h_parameterization + [None]*n_new_atoms
    new_manager = manager(
      pdb_hierarchy = pdb_hierarchy,
      geometry_restraints = geometry_restraints,
      use_ideal_bonds_angles = self.use_ideal_bonds_angles,
      process_manager = False)
    new_manager.h_parameterization = new_h_parameterization
    new_manager.parameterization_cpp = \
      self.get_parameterization_cpp(h_parameterization = new_h_parameterization)
    return new_manager

  def deep_copy(self):
    size = self.pdb_hierarchy.atoms_size()
    return self.select(selection=flex.bool(size, True))

  def select(self, selection):
    h_parameterization = self.copy_h_parameterization(self.h_parameterization)
    new_manager = manager(
      pdb_hierarchy = self.pdb_hierarchy.deep_copy(),
      geometry_restraints = self.geometry_restraints,
      use_ideal_bonds_angles = self.use_ideal_bonds_angles,
      process_manager = False)
    # Properties from current manager
    n_atoms = new_manager.pdb_hierarchy.atoms_size()
    iselection_original = new_manager.pdb_hierarchy.atoms().extract_i_seq()
    # Properties for the new (selected) manager
    #new_hierachy = new_manager.pdb_hierarchy.deep_copy()
    new_hierachy = new_manager.pdb_hierarchy.select(selection)
    #new_hierachy.select(selection)

    new_geometry_restraints = new_manager.geometry_restraints.select(selection)
    pdb_atoms = new_hierachy.atoms()
    pdb_atoms.reset_i_seq()
    hd_selection_new = new_hierachy.atom_selection_cache().\
          selection("element H or element D")
    new_manager.pdb_hierarchy = new_hierachy

    new_manager.hd_selection = hd_selection_new
    new_manager.not_hd_selection = ~hd_selection_new
    new_manager.geometry_restraints = new_geometry_restraints
    if selection.count(True) == n_atoms:
      new_h_parameterization = h_parameterization
    else:
      new_h_parameterization = self.reindex_h_parameterization(
        selection           = selection,
        iselection_original = iselection_original,
        n_atoms             = n_atoms,
        h_para              = h_parameterization)
    new_manager.h_parameterization = new_h_parameterization
    new_manager.parameterization_cpp = \
      self.get_parameterization_cpp(h_parameterization = new_h_parameterization)
    return new_manager

  def reindex_h_parameterization(
      self, selection, iselection_original, n_atoms, h_para):
    iselection = selection.iselection().as_int()
    r_a = list(reindexing_array(n_atoms, iselection))
    reindexing_dict = {}
    for index, i_seq in enumerate(iselection_original):
      reindexing_dict[iselection_original[index]] = r_a[index]
    # Change h_parameterization (contains i_seq --> have to be updated)
    new_h_parameterization = []
    for index, rc in enumerate(h_para):
      key = iselection_original[index]
      # No entry for non-selected atoms (H or non-H)
      if reindexing_dict[key] == n_atoms: continue
      # Non-H atoms (if included in selection) have entry None
      if rc is None:
        new_h_parameterization.append(None)
        continue
      # For other entries: 2 possibilities
      # a) update all i_seqs according to reindexing dictionary
      # b) if any neighbors of H is not in selection --> change this
      #    entry to None (because if neighbor is missing, H cannot be built)
      if reindexing_dict[rc.ih] != n_atoms:
        rc.ih = reindexing_dict[rc.ih]
      else:
        new_h_parameterization.append(None)
        continue
      if reindexing_dict[rc.a0] != n_atoms:
        rc.a0 = reindexing_dict[rc.a0]
      else:
        new_h_parameterization.append(None)
        continue
      if reindexing_dict[rc.a1] != n_atoms:
        rc.a1 = reindexing_dict[rc.a1]
      else:
        new_h_parameterization.append(None)
        continue
      if reindexing_dict[rc.a2] != n_atoms:
        rc.a2 = reindexing_dict[rc.a2]
      else:
        new_h_parameterization.append(None)
        continue
      # a3 only exists for htype "3neigbs"
      if rc.a3 != -1:
        if reindexing_dict[rc.a3] != n_atoms:
          rc.a3 = reindexing_dict[rc.a3]
        else:
          new_h_parameterization.append(None)
          continue
      new_h_parameterization.append(rc)
    return new_h_parameterization

  def get_parameterization_cpp(self, h_parameterization):
    parameterization_cpp = []
    for hp in h_parameterization:
      if (hp is not None):
        parameterization_cpp.append(hp)
    return parameterization_cpp

  def refinable_parameters_init(self):
    return flex.double(self.n_parameters(), 0)

  def n_parameters(self):
    return self.not_hd_selection.count(True)*3

  def idealize_riding_h_positions(self,
               sites_cart = None,
               pdb_hierarchy = None,
               xray_structure = None,
               selection_bool = None):
    assert [sites_cart, pdb_hierarchy, xray_structure].count(None) in [2,3]
    # prevent idealization if refinement flag for xyz is False
    if selection_bool:
      para_cpp = list()
      for p in self.parameterization_cpp:
        if selection_bool[p.ih]:
          para_cpp.append(p)
      self.parameterization_cpp = para_cpp
    #
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
      print(m, file=log)
      print('  Number of H atoms not used in the riding H procedure: ', \
        number, '\n', file=log)
      for ih in unk_list:
        atom = atoms[ih]
        labels = atom.fetch_labels()
        if (labels.resname in ['DOD', 'WAT', 'HOH']):
          continue
        print('\t', atom.id_str(), file=log)

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


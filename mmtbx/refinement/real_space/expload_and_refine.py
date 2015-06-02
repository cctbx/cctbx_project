from __future__ import division
import mmtbx.refinement.real_space.individual_sites
from cctbx import maptbx
from libtbx import adopt_init_args

class scorer(object):
  def __init__(self, pdb_hierarchy, unit_cell, map_data):
    adopt_init_args(self, locals())
    self.sites_cart = self.pdb_hierarchy.atoms().extract_xyz()
    self.target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = self.sites_cart)

  def update(self, sites_cart):
    target = maptbx.real_space_target_simple(
      unit_cell   = self.unit_cell,
      density_map = self.map_data,
      sites_cart  = sites_cart)
    if(target > self.target):
      self.target = target
      self.sites_cart = sites_cart
    print self.target, target # XXX for debugging

class run(object):
  def __init__(
        self,
        xray_structure,
        pdb_hierarchy,
        map_data,
        restraints_manager,
        number_of_macro_cycles=5,
        target_bond_rmsd=0.02,
        target_angle_rmsd=2.0,
        number_of_trials=20,
        xyz_shake=10.0,
        states=None):
    self.xray_structure = xray_structure
    self.pdb_hierarchy = pdb_hierarchy
    self.states = states
    # Refine without exploding - get refined starting model
    ro=None
    for macro_cycle in range(1, number_of_macro_cycles):
      scale = float(number_of_macro_cycles)/macro_cycle
      ro = mmtbx.refinement.real_space.individual_sites.easy(
        map_data                    = map_data,
        xray_structure              = self.xray_structure,
        pdb_hierarchy               = self.pdb_hierarchy,
        geometry_restraints_manager = restraints_manager,
        rms_bonds_limit             = target_bond_rmsd*scale,
        rms_angles_limit            = target_angle_rmsd*scale,
        max_iterations              = 50,
        selection                   = None,
        states_accumulator          = self.states, # this is good for demo, but makes it slower
        log                         = None)
      self.xray_structure = ro.xray_structure
    if ro:
      weight = ro.w
    else:
      weight=None
    # Expload and refine
    sc = scorer(
      pdb_hierarchy = self.pdb_hierarchy,
      unit_cell     = self.xray_structure.unit_cell(),
      map_data      = map_data)
    for trial in xrange(number_of_trials):
      self.xray_structure.shake_sites_in_place(
        rms_difference = None,
        mean_distance  = xyz_shake)
      ro = mmtbx.refinement.real_space.individual_sites.easy(
        map_data                    = map_data,
        xray_structure              = self.xray_structure,
        pdb_hierarchy               = self.pdb_hierarchy,
        geometry_restraints_manager = restraints_manager,
        rms_bonds_limit             = target_bond_rmsd,
        rms_angles_limit            = target_angle_rmsd,
        max_iterations              = 250,
        selection                   = None,
        w                           = weight,
        log                         = None)
      self.xray_structure = ro.xray_structure
      sc.update(sites_cart = self.xray_structure.sites_cart())
      self.xray_structure = self.xray_structure.replace_sites_cart(
        new_sites=sc.sites_cart)
      self.states.add(sites_cart = sc.sites_cart)
    print "LOOK:",maptbx.real_space_target_simple( # XXX for debugging
        unit_cell   = self.xray_structure.unit_cell(),
        density_map = map_data,
        sites_cart  = self.xray_structure.sites_cart())
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

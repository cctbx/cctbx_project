from __future__ import division
from cctbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.refinement.real_space
from mmtbx.refinement.real_space import individual_sites
import math
from cctbx import maptbx
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager

import boost.python
ext = boost.python.import_ext("mmtbx_rotamer_fit_ext")

def flatten(l):
  if l is None: return None
  return sum(([x] if not (isinstance(x, list) or isinstance(x, flex.size_t))
    else flatten(x) for x in l), [])

class run(object):
  def __init__(self,
               residue,
               unit_cell,
               target_map,
               mon_lib_srv,
               rotamer_manager,
               sin_cos_table,
               backbone_sample=True):
    adopt_init_args(self, locals())
    # Initial state
    rotamer_start = rotamer_manager.rotamer(residue=self.residue)
    sites_cart_start=self.residue.atoms().extract_xyz()
    target_start = self.get_target_value(sites_cart=sites_cart_start)
    # Actual calculations
    self.chi_angles = self.rotamer_manager.get_chi_angles(
      resname=self.residue.resname)
    co = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = True)
    if(backbone_sample):
      self.fit_c_beta(c_beta_rotation_cluster = co.clusters[0])
    self.fit_side_chain(clusters = co.clusters[1:])
    # Final state
    rotamer_final = rotamer_manager.rotamer(residue=self.residue)
    target_final = self.get_target_value(
      sites_cart=self.residue.atoms().extract_xyz())
    # Sanity and consistency check
    #XXX not always hold due to approx fast math assert rotamer_final != "OUTLIER"
    # Potentially this will keep an OUTLIER if no better fit found
    if(target_start > target_final):
      self.residue.atoms().set_xyz(sites_cart_start)

  def get_target_value(self, sites_cart, selection=None):
    if(selection is None):
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = self.target_map,
        sites_cart  = sites_cart)
    else:
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = self.target_map,
        sites_cart  = sites_cart,
        selection   = selection)

  def fit_side_chain(self, clusters):
    rotamer_iterator = \
      mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = self.residue)
    if(rotamer_iterator is None): return
    selection = flex.size_t(flatten(clusters[0].vector))
    start_target_value = self.get_target_value(
      sites_cart = self.residue.atoms().extract_xyz(),
      selection  = selection)
    sites_cart_first_rotamer = list(rotamer_iterator)[0][1]
    self.residue.atoms().set_xyz(sites_cart_first_rotamer)
    axes = []
    atr = []
    for i, angle in enumerate(self.chi_angles[0]):
      cl = clusters[i]
      axes.append(flex.size_t(cl.axis))
      atr.append(flex.size_t(cl.atoms_to_rotate))
    ro = ext.fit(
      target_value             = start_target_value,
      axes                     = axes,
      rotatable_points_indices = atr,
      angles_array             = self.chi_angles,
      density_map              = self.target_map,
      all_points               = self.residue.atoms().extract_xyz(),
      unit_cell                = self.unit_cell,
      selection                = selection,
      sin_table                = self.sin_cos_table.sin_table,
      cos_table                = self.sin_cos_table.cos_table,
      step                     = self.sin_cos_table.step,
      n                        = self.sin_cos_table.n)
    sites_cart_result = ro.result()
    if(sites_cart_result.size()>0):
      self.residue.atoms().set_xyz(sites_cart_result)

  def fit_c_beta(self, c_beta_rotation_cluster):
    selection = flex.size_t(c_beta_rotation_cluster.selection)
    sites_cart = self.residue.atoms().extract_xyz()
    start_target_value = self.get_target_value(
      sites_cart = sites_cart,
      selection  = selection)
    ro = ext.fit(
      target_value             = start_target_value,
      axes                     = [c_beta_rotation_cluster.axis],
      rotatable_points_indices = [c_beta_rotation_cluster.atoms_to_rotate],
      angles_array             = [[i*math.pi/180] for i in range(-20,21,1)],
      density_map              = self.target_map,
      all_points               = sites_cart,
      unit_cell                = self.unit_cell,
      selection                = selection,
      sin_table                = self.sin_cos_table.sin_table,
      cos_table                = self.sin_cos_table.cos_table,
      step                     = self.sin_cos_table.step,
      n                        = self.sin_cos_table.n)
    sites_cart_result = ro.result()
    if(sites_cart_result.size()>0):
      self.residue.atoms().set_xyz(sites_cart_result)

class run_with_minimization(object):
  def __init__(self,
               target_map,
               residue,
               xray_structure,
               mon_lib_srv,
               rotamer_manager,
               # This is cctbx.geometry_restraints.manager.manager
               geometry_restraints_manager,
               real_space_gradients_delta,
               selection_radius = 5,
               rms_bonds_limit = 0.03, # XXX probably needs to be much lower
               rms_angles_limit = 3.0, # XXX
               backbone_sample_angle=None,
               allow_modified_residues=False):
    adopt_init_args(self, locals())
    # load rotamer manager
    self.rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load()
    # pre-compute sin and cos tables
    self.sin_cos_table = scitbx.math.sin_cos_table(n=10000)
    self.backbone_atom_names = ["N", "CA", "O", "CB", "C"]
    self.residue_iselection = self.residue.atoms().extract_i_seq()
    assert (not self.residue_iselection.all_eq(0))
    self.residue_selection = flex.bool(
      xray_structure.scatterers().size(), self.residue_iselection)
    self.residue_backbone_selection = flex.size_t()
    for atom in self.residue.atoms():
      if(atom.name.strip() in self.backbone_atom_names):
        self.residue_backbone_selection.append(atom.i_seq)
    self.residue_backbone_selection = flex.bool(
      xray_structure.scatterers().size(), self.residue_backbone_selection)
    self.target_map_work = target_map
    self.fit_backbone()
    negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
      xray_structure          = self.xray_structure,
      selection_within_radius = self.selection_radius,
      iselection              = self.residue.atoms().extract_i_seq())
    self.target_map_work = mmtbx.refinement.real_space.\
      negate_map_around_selected_atoms_except_selected_atoms(
        xray_structure   = self.xray_structure,
        map_data         = target_map,
        negate_selection = negate_selection,
        atom_radius      = 1.5)
    self.fit_rotamers()

  def fit_backbone(self):
    # move in place (pure geometry regularizaition of residue in question)
    self.real_space_refine(optimize_weight=False, start_trial_weight_value=0,
      use_selection_real_space=False)
    # fit n-c-o-ca-cb only (ignore side chain!). XXX BAD: amino-acid specific!
    self.grid_sample_around_c_n_axis()
    # fine-tune
    self.real_space_refine(optimize_weight=True, start_trial_weight_value=50,
      use_selection_real_space=True)

  def fit_rotamers(self):
    sps = self.xray_structure.special_position_settings()
    mmtbx.refinement.real_space.fit_residue.run(
      target_map      = self.target_map_work,
      mon_lib_srv     = self.mon_lib_srv,
      unit_cell       = self.xray_structure.unit_cell(),
      residue         = self.residue,
      sin_cos_table   = self.sin_cos_table,
      rotamer_manager = self.rotamer_manager)
    sites_cart_poor = self.xray_structure.sites_cart()
    sites_cart_poor.set_selected(self.residue_iselection,
      self.residue.atoms().extract_xyz())
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_poor)

  def grid_sample_around_c_n_axis(self):
    sps = self.xray_structure.special_position_settings()
    scorer = mmtbx.refinement.real_space.score(
      target_map = self.target_map_work,
      residue    = self.residue,
      unit_cell  = self.xray_structure.unit_cell())
    def get_cluster(self):
      axis=[]
      atoms_to_rotate=[]
      use_in_target_selection = flex.size_t()
      counter = 0
      for atom in self.residue.atoms():
        if(atom.name.strip() in ["N", "C"]):
          axis.append(counter)
        else:
          atoms_to_rotate.append(counter)
        if(atom.name.strip() in self.backbone_atom_names):
          use_in_target_selection.append(counter)
        counter += 1
      return mmtbx.refinement.real_space.cluster(
        axis            = axis,
        atoms_to_rotate = atoms_to_rotate,
        selection       = use_in_target_selection)
    cl = get_cluster(self)
    residue_sites_cart = self.residue.atoms().extract_xyz()
    scorer.reset_with(
      sites_cart = residue_sites_cart,
      selection  = cl.selection)
    angle_start = 0
    angle_end = 360
    if (self.backbone_sample_angle is not None) :
      assert (self.backbone_sample_angle > 0)
      angle_start = - self.backbone_sample_angle
      angle_end = self.backbone_sample_angle
    mmtbx.refinement.real_space.torsion_search(
      clusters   = [cl],
      sites_cart = residue_sites_cart,
      scorer     = scorer,
      start      = 0,
      stop       = 360,
      step       = 1)
    self.residue.atoms().set_xyz(new_xyz=scorer.sites_cart)
    selection = self.residue.atoms().extract_i_seq()
    sites_cart_poor = self.xray_structure.sites_cart()
    sites_cart_poor.set_selected(selection, scorer.sites_cart)
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_poor)

  def real_space_refine(self, optimize_weight, start_trial_weight_value,
        use_selection_real_space):
    brm = individual_sites.box_refinement_manager(
      xray_structure              = self.xray_structure,
      target_map                  = self.target_map_work,
      geometry_restraints_manager = self.geometry_restraints_manager,
      real_space_gradients_delta  = 1./4,
      max_iterations              = 500)
    brm.refine(
      selection                = self.residue_selection,
      optimize_weight          = optimize_weight,
      start_trial_weight_value = start_trial_weight_value,
      selection_buffer_radius  = self.selection_radius,
      box_cushion              = 2,
      rms_bonds_limit          = self.rms_bonds_limit,
      rms_angles_limit         = self.rms_angles_limit)
    self.xray_structure = brm.xray_structure
    self.residue.atoms().set_xyz(brm.sites_cart.select(self.residue_iselection))

def get_rotamer_iterator(mon_lib_srv, residue):
  rotamer_iterator = mon_lib_srv.rotamer_iterator(
    fine_sampling = True,
    comp_id=residue.resname,
    atom_names=residue.atoms().extract_name(),
    sites_cart=residue.atoms().extract_xyz())
  if (rotamer_iterator is None) :
    return None
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

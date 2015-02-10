from __future__ import division
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.lbfgs
from libtbx import adopt_init_args
import iotbx.pdb
import mmtbx.refinement.real_space
import cctbx.geometry_restraints.flags
from mmtbx.refinement.real_space import individual_sites
import sys

class run(object):
  def __init__(self,
               target_map,
               residue,
               xray_structure,
               mon_lib_srv,
               rotamer_manager,
               # FIXME this actually needs to be an mmtbx.restraints.manager
               # object, not cctbx.geometry_restraints.manager.manager - the
               # keyword should be changed to reflect this
               geometry_restraints_manager,
               real_space_gradients_delta,
               selection_radius = 5,
               rms_bonds_limit = 0.03, # XXX probably needs to be much lower
               rms_angles_limit = 3.0, # XXX
               backbone_sample_angle=None,
               allow_modified_residues=False):
    adopt_init_args(self, locals())
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
    mmtbx.refinement.real_space.fit_residue.manager(
      target_map                = self.target_map_work,
      mon_lib_srv               = self.mon_lib_srv,
      special_position_settings = sps,
      residue                   = self.residue,
      use_clash_filter          = True,
      sites_cart_all            = self.xray_structure.sites_cart(),
      rotamer_manager           = self.rotamer_manager,
      allow_modified_residues   = self.allow_modified_residues)
    sites_cart_poor = self.xray_structure.sites_cart()
    sites_cart_poor.set_selected(self.residue_iselection,
      self.residue.atoms().extract_xyz())
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_poor)

  def grid_sample_around_c_n_axis(self):
    sps = self.xray_structure.special_position_settings()
    scorer = mmtbx.refinement.real_space.score(
      target_map                = self.target_map_work,
      sites_cart_all            = self.xray_structure.sites_cart(),
      residue                   = self.residue,
      special_position_settings = sps)
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
        selection       = use_in_target_selection,
        start = 0,
        stop = 360,
        step=1)
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
      scorer     = scorer)
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
      geometry_restraints_manager = self.geometry_restraints_manager.geometry,
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

class manager(object):
  def __init__(self,
               target_map,
               residue,
               special_position_settings,
               mon_lib_srv,
               rotamer_manager,
               sites_cart_all=None,
               use_clash_filter = False,
               use_slope=True,
               debug=False,
               use_torsion_search=True,
               use_rotamer_iterator=True,
               torsion_search_backrub_start =-25,
               torsion_search_backrub_stop  = 25,
               torsion_search_backrub_step  = 5,
               torsion_search_all_start = 0,
               torsion_search_all_stop  = 360,
               torsion_search_all_step  = 1,
               torsion_search_local_start = -50,
               torsion_search_local_stop  = 50,
               torsion_search_local_step  = 5,
               allow_modified_residues    =  False,
               backbone_sample            = True,
               log                        = None):
    adopt_init_args(self, locals())
    if(self.log is None): self.log = sys.stdout
    if(not allow_modified_residues):
      get_class = iotbx.pdb.common_residue_names_get_class
      assert get_class(residue.resname) == "common_amino_acid", residue.resname
    if(self.use_torsion_search):
      start = self.torsion_search_all_start
      stop  = self.torsion_search_all_stop
      step  = self.torsion_search_all_step
    else:
      start = self.torsion_search_local_start
      stop  = self.torsion_search_local_stop
      step  = self.torsion_search_local_step
    self.clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = self.backbone_sample,
      torsion_search_backbone_start  = torsion_search_backrub_start,
      torsion_search_backbone_stop   = torsion_search_backrub_stop,
      torsion_search_backbone_step   = torsion_search_backrub_step,
      torsion_search_sidechain_start = start,
      torsion_search_sidechain_stop  = stop,
      torsion_search_sidechain_step  = step).clusters
    self.vector_selections = self.clusters[0].vector
    self.fit(debug=debug)

  def get_scorer(self, vector=None, use_binary=False, use_clash_filter=False):
    return mmtbx.refinement.real_space.score(
      target_map                = self.target_map,
      sites_cart_all            = self.sites_cart_all,
      special_position_settings = self.special_position_settings,
      residue                   = self.residue,
      vector                    = vector,
      use_binary                = use_binary,
      use_clash_filter          = use_clash_filter,
      rotamer_eval              = self.rotamer_manager)

  def fit(self, debug=True):
    if(self.use_rotamer_iterator):
      rotamer_iterator = get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = self.residue)
      if(rotamer_iterator is not None):
        score_rotamers = self.get_scorer(
          vector           = self.vector_selections,
          use_clash_filter = self.use_clash_filter,
          use_binary       = True)
        cntr=0
        for rotamer, rotamer_sites_cart in rotamer_iterator:
          if(debug): print rotamer.id,"-"*50
          if(cntr==0):
            score_rotamers.reset_with(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 0, rotamer.id
          else:
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 0, rotamer.id
          cntr+=1
          if(self.use_torsion_search):
            score_rotamer = self.get_scorer()
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            mmtbx.refinement.real_space.torsion_search(
              clusters   = self.clusters,
              sites_cart = rotamer_sites_cart.deep_copy(),
              scorer     = score_rotamer)
            if(debug): print score_rotamers.target, 1
            score_rotamers.update(
              sites_cart = score_rotamer.sites_cart.deep_copy(), tmp=rotamer.id)
            if(debug): print score_rotamers.target, 2
          else:
            if(debug): print score_rotamers.target, 3
            score_rotamers.update(sites_cart = rotamer_sites_cart.deep_copy())
            if(debug): print score_rotamers.target, 4
          if(debug): print rotamer.id, score_rotamers.target, score_rotamers.tmp
        if(debug): print "final:", score_rotamers.target, score_rotamers.tmp
        if(score_rotamers.sites_cart is not None):
          self.residue.atoms().set_xyz(new_xyz = score_rotamers.sites_cart)
          if(self.sites_cart_all is not None):
            self.sites_cart_all = self.sites_cart_all.set_selected(
              self.residue.atoms().extract_i_seq(), score_rotamers.sites_cart)
    elif(self.use_torsion_search):
      score_residue = self.get_scorer()
      mmtbx.refinement.real_space.torsion_search(
        clusters   = self.clusters,
        sites_cart = self.residue.atoms().extract_xyz(),
        scorer     = score_residue)
      self.residue.atoms().set_xyz(new_xyz=score_residue.sites_cart)

  def rigid_body_refine(self, max_iterations = 250):
    import mmtbx.geometry_restraints
    import cctbx.geometry_restraints.manager
    import mmtbx.refinement.real_space.rigid_body
    reference_sites = flex.vec3_double()
    reference_selection = flex.size_t()
    cntr = 0
    for atom in self.residue.atoms():
      if(atom.name.strip().upper() in ["N", "C"]):
        reference_sites.append(atom.xyz)
        reference_selection.append(cntr)
        cntr += 1
    generic_restraints_manager = mmtbx.geometry_restraints.manager()
    restraints_manager = cctbx.geometry_restraints.manager.manager(
      generic_restraints_manager = generic_restraints_manager)
    restraints_manager.generic_restraints_manager.reference_manager.\
      add_coordinate_restraints(
        sites_cart = reference_sites,
        selection  = reference_selection,
        sigma      = 0.5)
    flags = cctbx.geometry_restraints.flags.flags(generic_restraints=True)
    lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    minimized = mmtbx.refinement.real_space.rigid_body.refine(
      residue                     = self.residue,
      density_map                 = self.target_map,
      geometry_restraints_manager = restraints_manager,
      real_space_target_weight    = 1,
      real_space_gradients_delta  = 2.0*0.25,
      lbfgs_termination_params    = lbfgs_termination_params,
      unit_cell                   = self.unit_cell,
      cctbx_geometry_restraints_flags = flags)
    return minimized.sites_cart_residue

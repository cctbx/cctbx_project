from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.refinement.real_space
from mmtbx.refinement.real_space import individual_sites
import math
from cctbx import maptbx
import scitbx.math
import mmtbx.idealized_aa_residues.rotamer_manager
import iotbx.pdb

import boost.python
ext = boost.python.import_ext("mmtbx_rotamer_fit_ext")

def flatten(l):
  if l is None: return None
  return sum(([x] if not (isinstance(x, list) or isinstance(x, flex.size_t))
    else flatten(x) for x in l), [])

class run(object):
  def __init__(self,
               residue,
               mon_lib_srv,
               rotamer_manager,
               sin_cos_table,
               vdw_radii=None,
               xyzrad_bumpers=None,
               target_map=None,
               target_map_for_cb=None,
               unit_cell=None,
               backbone_sample=True,
               accept_only_if_max_shift_is_smaller_than=None):
    adopt_init_args(self, locals())
    #
    if(target_map is None):
      assert not backbone_sample
    # Initial state
    rotamer_start = rotamer_manager.rotamer(residue=self.residue)
    sites_cart_start = self.residue.atoms().extract_xyz()
    self.co = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = True)
    if(target_map is not None):
      target_start = self.get_target_value(sites_cart = sites_cart_start)
    # Actual calculations
    self.chi_angles = self.rotamer_manager.get_chi_angles(
      resname = self.residue.resname)
    if(len(self.co.clusters)>0):
      if(backbone_sample):
        self.fit_c_beta(c_beta_rotation_cluster = self.co.clusters[0])
      self.fit_side_chain(clusters = self.co.clusters[1:])
      # Final state
      if(target_map is not None):
        target_final = self.get_target_value(
          sites_cart=self.residue.atoms().extract_xyz())
        if(target_start > target_final):
          self.residue.atoms().set_xyz(sites_cart_start)

  def get_target_value(self, sites_cart, selection=None, target_map=None):
    if(target_map is None): target_map = self.target_map
    if(selection is None):
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = target_map,
        sites_cart  = sites_cart)
    else:
      return maptbx.real_space_target_simple(
        unit_cell   = self.unit_cell,
        density_map = target_map,
        sites_cart  = sites_cart,
        selection   = selection)

  def fit_side_chain(self, clusters):
    rotamer_iterator = \
      mmtbx.refinement.real_space.fit_residue.get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = self.residue)
    if(rotamer_iterator is None): return
    #selection_rsr = flex.size_t(flatten(clusters[0].vector))
    selection_clash = self.co.clash_eval_selection
    selection_rsr   = self.co.rsr_eval_selection
    if(self.target_map is not None):
      start_target_value = self.get_target_value(
        sites_cart = self.residue.atoms().extract_xyz(),
        selection  = selection_rsr)
    sites_cart_start = self.residue.atoms().extract_xyz()
    sites_cart_first_rotamer = list(rotamer_iterator)[0][1]
    self.residue.atoms().set_xyz(sites_cart_first_rotamer)
    axes = []
    atr = []
    for i, angle in enumerate(self.chi_angles[0]):
      cl = clusters[i]
      axes.append(flex.size_t(cl.axis))
      atr.append(flex.size_t(cl.atoms_to_rotate))
    sites = self.residue.atoms().extract_xyz()
    if(self.target_map is not None and self.xyzrad_bumpers is not None):
      # Get vdW radii
      radii = flex.double()
      atom_names = []
      for a in self.residue.atoms():
        atom_names.append(a.name.strip())
      converter = iotbx.pdb.residue_name_plus_atom_names_interpreter(
        residue_name=self.residue.resname, atom_names = atom_names)
      mon_lib_names = converter.atom_name_interpretation.mon_lib_names()
      for n in mon_lib_names:
        try: radii.append(self.vdw_radii[n.strip()]-0.25)
        except KeyError: radii.append(1.5) # XXX U, Uranium, OXT are problems!
      #
      xyzrad_residue = ext.xyzrad(sites_cart = sites, radii = radii)
      #
      ro = ext.fit(
        target_value             = start_target_value,
        xyzrad_bumpers           = self.xyzrad_bumpers,
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        density_map              = self.target_map,
        all_points               = xyzrad_residue,
        unit_cell                = self.unit_cell,
        selection_clash          = selection_clash,
        selection_rsr            = selection_rsr,
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    elif(self.target_map is not None and self.xyzrad_bumpers is None):
      ro = ext.fit(
        target_value             = start_target_value,
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        density_map              = self.target_map,
        all_points               = sites,
        unit_cell                = self.unit_cell,
        selection                = selection_rsr,
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    else:
      ro = ext.fit(
        sites_cart_start         = sites_cart_start.deep_copy(),
        axes                     = axes,
        rotatable_points_indices = atr,
        angles_array             = self.chi_angles,
        all_points               = self.residue.atoms().extract_xyz(),
        sin_table                = self.sin_cos_table.sin_table,
        cos_table                = self.sin_cos_table.cos_table,
        step                     = self.sin_cos_table.step,
        n                        = self.sin_cos_table.n)
    sites_cart_result = ro.result()
    if(sites_cart_result.size()>0):
      dist = None
      if(self.accept_only_if_max_shift_is_smaller_than is not None):
        dist = flex.max(flex.sqrt((sites_cart_start - sites_cart_result).dot()))
      if(dist is None):
        self.residue.atoms().set_xyz(sites_cart_result)
      else:
        if(dist is not None and
           dist < self.accept_only_if_max_shift_is_smaller_than):
          self.residue.atoms().set_xyz(sites_cart_result)
        else:
          self.residue.atoms().set_xyz(sites_cart_start)
    else:
      self.residue.atoms().set_xyz(sites_cart_start)

  def fit_c_beta(self, c_beta_rotation_cluster):
    selection = flex.size_t(c_beta_rotation_cluster.selection)
    sites_cart = self.residue.atoms().extract_xyz()
    sites_cart_start = sites_cart.deep_copy() # XXX
    start_target_value = self.get_target_value(
      sites_cart = sites_cart,
      selection  = selection,
      target_map = self.target_map_for_cb)
    ro = ext.fit(
      target_value             = start_target_value+1.e-6,
      axes                     = [c_beta_rotation_cluster.axis],
      rotatable_points_indices = [c_beta_rotation_cluster.atoms_to_rotate],
      angles_array             = [[i*math.pi/180] for i in range(-20,21,1)],
      density_map              = self.target_map_for_cb,
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
    else:
      self.residue.atoms().set_xyz(sites_cart_start)

class run_with_minimization(object):
  def __init__(self,
               target_map,
               residue,
               vdw_radii,
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
    self.target_map_orig = target_map.deep_copy()
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
    self.real_space_refine(optimize_weight=False, start_trial_weight_value=0)
    # fit n-c-o-ca-cb only (ignore side chain!). XXX BAD: amino-acid specific!
    self.grid_sample_around_c_n_axis()
    # fine-tune
    self.real_space_refine(optimize_weight=True, start_trial_weight_value=50)

  def fit_rotamers(self):
    sps = self.xray_structure.special_position_settings()
    mmtbx.refinement.real_space.fit_residue.run(
      vdw_radii         = self.vdw_radii,
      target_map        = self.target_map_work,
      target_map_for_cb = self.target_map_orig,
      mon_lib_srv       = self.mon_lib_srv,
      unit_cell         = self.xray_structure.unit_cell(),
      residue           = self.residue,
      sin_cos_table     = self.sin_cos_table,
      rotamer_manager   = self.rotamer_manager)
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
    scorer.reset(
      sites_cart = residue_sites_cart,
      selection  = cl.selection)
    angle_start = 0
    angle_end = 360
    if (self.backbone_sample_angle is not None):
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

  def real_space_refine(self, optimize_weight, start_trial_weight_value):
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
  if (rotamer_iterator is None):
    return None
  if (rotamer_iterator.problem_message is not None):
    return None
  if (rotamer_iterator.rotamer_info is None):
    return None
  return rotamer_iterator

class tune_up(object):
  def __init__(self,
               target_map,
               residue,
               mon_lib_srv,
               rotamer_manager,
               unit_cell,
               torsion_search_start = -20,
               torsion_search_stop  = 20,
               torsion_search_step  = 2):
    adopt_init_args(self, locals())
    self.clusters = mmtbx.refinement.real_space.aa_residue_axes_and_clusters(
      residue         = self.residue,
      mon_lib_srv     = self.mon_lib_srv,
      backbone_sample = False).clusters
    score_residue = mmtbx.refinement.real_space.score3(
      unit_cell    = self.unit_cell,
      target_map   = self.target_map,
      residue      = self.residue,
      rotamer_eval = self.rotamer_manager)
    mmtbx.refinement.real_space.torsion_search(
      clusters   = self.clusters,
      sites_cart = self.residue.atoms().extract_xyz(),
      scorer     = score_residue,
      start      = self.torsion_search_start,
      stop       = self.torsion_search_stop,
      step       = self.torsion_search_step)
    self.residue.atoms().set_xyz(new_xyz=score_residue.sites_cart)

#
# These functions are not used anywhere. And not tested anymore.
# They are here as an example of correct backrub move, according to
# original paper https://doi.org/10.1016/j.str.2005.10.007
# Unfortunately, for proper backrub move we need previous and next residues,
# but current code is build under assumption that one residue is enough for
# rotamer fitting. One will have to reconsider this idea and do some changes
# to make it possible to do proper backrub move.
#
def _find_theta(ap1, ap2, cur_xyz, needed_xyz):
  from mmtbx.building.loop_closure.ccd import ccd_python
  f, s_home, r_norm, r_home = ccd_python._get_f_r_s(
      axis_point_1=ap1,
      axis_point_2=ap2,
      moving_coor=cur_xyz,
      fixed_coor=needed_xyz)
  b = list(2*r_norm*(f.dot(r_home)))[0]
  c = list(2*r_norm*(f.dot(s_home)))[0]
  znam = math.sqrt(b*b+c*c)
  sin_alpha = c/znam
  cos_alpha = b/znam
  alpha = math.atan2(sin_alpha, cos_alpha)
  return math.degrees(alpha)

def backrub_move(
    prev_res,
    cur_res,
    next_res,
    angle,
    move_oxygens=False,
    accept_worse_rama=False,
    rotamer_manager=None,
    rama_manager=None):
  import boost.python
  ext = boost.python.import_ext("mmtbx_validation_ramachandran_ext")
  from mmtbx_validation_ramachandran_ext import rama_eval
  from scitbx.matrix import rotate_point_around_axis
  from mmtbx.conformation_dependent_library.multi_residue_class import ThreeProteinResidues, \
      RestraintsRegistry

  if abs(angle) < 1e-4:
    return
  if prev_res is None or next_res is None:
    return
  saved_res = [{},{},{}]
  for i, r in enumerate([prev_res, cur_res, next_res]):
    for a in r.atoms():
      saved_res[i][a.name.strip()] = a.xyz
  if rotamer_manager is None:
    rotamer_manager = RotamerEval()
  prev_ca = prev_res.find_atom_by(name=" CA ")
  cur_ca = cur_res.find_atom_by(name=" CA ")
  next_ca = next_res.find_atom_by(name=" CA ")
  if prev_ca is None or next_ca is None or cur_ca is None:
    return
  atoms_to_move = []
  atoms_to_move.append(prev_res.find_atom_by(name=" C  "))
  atoms_to_move.append(prev_res.find_atom_by(name=" O  "))
  for atom in cur_res.atoms():
    atoms_to_move.append(atom)
  atoms_to_move.append(next_res.find_atom_by(name=" N  "))
  for atom in atoms_to_move:
    assert atom is not None
    new_xyz = rotate_point_around_axis(
        axis_point_1 = prev_ca.xyz,
        axis_point_2 = next_ca.xyz,
        point        = atom.xyz,
        angle        = angle,
        deg          = True)
    atom.xyz = new_xyz
  if move_oxygens:
    registry = RestraintsRegistry()
    if rama_manager is None:
      rama_manager = rama_eval()
    tpr = ThreeProteinResidues(geometry=None, registry=registry)
    tpr.append(prev_res)
    tpr.append(cur_res)
    tpr.append(next_res)
    phi_psi_angles = tpr.get_phi_psi_angles()
    rama_key = tpr.get_ramalyze_key()
    ev_before = rama_manager.evaluate_angles(rama_key, phi_psi_angles[0], phi_psi_angles[1])
    theta1 = _find_theta(
        ap1 = prev_ca.xyz,
        ap2 = cur_ca.xyz,
        cur_xyz = prev_res.find_atom_by(name=" O  ").xyz,
        needed_xyz = saved_res[0]["O"])
    theta2 = _find_theta(
        ap1 = cur_ca.xyz,
        ap2 = next_ca.xyz,
        cur_xyz = cur_res.find_atom_by(name=" O  ").xyz,
        needed_xyz = saved_res[1]["O"])
    for a in [prev_res.find_atom_by(name=" C  "),
        prev_res.find_atom_by(name=" O  "),
        cur_res.find_atom_by(name=" C  ")]:
      new_xyz = rotate_point_around_axis(
              axis_point_1 = prev_ca.xyz,
              axis_point_2 = cur_ca.xyz,
              point        = a.xyz,
              angle        = theta1,
              deg          = True)
      a.xyz = new_xyz
    for a in [cur_res.find_atom_by(name=" C  "),
        cur_res.find_atom_by(name=" O  "),
        next_res.find_atom_by(name=" N  ")]:
      new_xyz = rotate_point_around_axis(
              axis_point_1 = cur_ca.xyz,
              axis_point_2 = next_ca.xyz,
              point        = a.xyz,
              angle        = theta2,
              deg          = True)
      a.xyz = new_xyz
    phi_psi_angles = tpr.get_phi_psi_angles()
    rama_key = tpr.get_ramalyze_key()
    ev_after = rama_manager.evaluate_angles(rama_key, phi_psi_angles[0], phi_psi_angles[1])
    if ev_before > ev_after and not accept_worse_rama:
      for a in [prev_res.find_atom_by(name=" C  "),
          prev_res.find_atom_by(name=" O  "),
          cur_res.find_atom_by(name=" C  ")]:
        new_xyz = rotate_point_around_axis(
                axis_point_1 = prev_ca.xyz,
                axis_point_2 = cur_ca.xyz,
                point        = a.xyz,
                angle        = -theta1,
                deg          = True)
        a.xyz = new_xyz
      for a in [cur_res.find_atom_by(name=" C  "),
          cur_res.find_atom_by(name=" O  "),
          next_res.find_atom_by(name=" N  ")]:
        new_xyz = rotate_point_around_axis(
                axis_point_1 = cur_ca.xyz,
                axis_point_2 = next_ca.xyz,
                point        = a.xyz,
                angle        = -theta2,
                deg          = True)
        a.xyz = new_xyz

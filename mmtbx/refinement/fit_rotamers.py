from mmtbx.refinement import print_statistics
from cctbx.array_family import flex
from libtbx import adopt_init_args
from mmtbx.command_line import lockit
from libtbx import adopt_init_args
from cctbx import maptbx
import scitbx.lbfgs
import iotbx.pdb
import mmtbx.monomer_library
import mmtbx.model
import mmtbx.refinement.real_space
from mmtbx import map_tools
from cctbx import maptbx
import math, sys
import mmtbx.monomer_library.server
import iotbx.pdb.atom_name_interpretation
import scitbx.graph.tardy_tree
import iotbx.phil
import mmtbx.monomer_library
import mmtbx.monomer_library.rotamer_utils

torsion_search_params_str = """\
torsion_search
  .style = box auto_align
{
  min_angle_between_solutions = 5
    .type = float
    .short_caption = Min. angle between solutions
  range_start = -40
    .type = float
  range_stop = 40
    .type = float
  step = 2
    .type = float
}
"""

local_fix_params_str = """\
  number_of_macro_cycles = 1
    .type = int
    .short_caption = Run overall real-space refinement (RSR)
  real_space_refine_overall = False
    .type = bool
  validate_change = True
    .type = bool
  exclude_hydrogens = True
    .type = bool
    .short_caption = Exclude hydrogens from RSR
  filter_residual_map_value = 2.0
    .type = float
  filter_2fofc_map = None
    .type = float
  target_map = 2mFo-DFc
    .type = str
  residual_map = mFo-DFc
    .type = str
  model_map = Fc
    .type = str
"""

master_params_str = """\
fit_side_chains
  .short_caption = Sidechain rotamer fitting
  .style = menu_item auto_align
{
  %s
  use_dihedral_restraints = False
    .type = bool
  ignore_water_when_move_sidechains = True
    .type = bool
    .short_caption = Ignore water when moving sidechains
  residue_iteration
    .style = box auto_align
  {
    poor_cc_threshold = 0.9
      .type = float
      .short_caption = Poor CC threshold
    real_space_refine_rotamer = True
      .type = bool
    real_space_refine_max_iterations = 25
      .type = int
    real_space_refine_target_weight = 100.
      .type = float
    use_rotamer_iterator = True
      .type = bool
    torsion_grid_search = True
      .type = bool
    ignore_alt_conformers = True
      .type = bool
      .short_caption = Ignore alternate conformers
    %s
  }
}
"""% (local_fix_params_str, torsion_search_params_str)

def master_params():
  return iotbx.phil.parse(input_string = master_params_str)

def torsion_search_params():
  return iotbx.phil.parse(input_string = torsion_search_params_str)

def tardy_model_one_residue(residue, mon_lib_srv, log = None):
  if(log is None): log = sys.stdout
  comp_comp_id = mon_lib_srv.get_comp_comp_id_direct(comp_id=residue.resname)
  # create external clusters to address planes with missing atoms
  external_clusters = []
  for plane in comp_comp_id.get_planes():
    plane_atom_names = []
    for plane_atom in plane.plane_atoms:
      plane_atom_names.append(plane_atom.atom_id)
    plane_atom_i_seqs = []
    if 0: print plane_atom_names
    for i_atom, atom in enumerate(residue.atoms()):
      if(atom.name.strip() in plane_atom_names):
        plane_atom_i_seqs.append(i_atom)
    if(len(plane_atom_i_seqs)>0):
      external_clusters.append(plane_atom_i_seqs)
  #
  if 0: print "external_clusters1:",external_clusters
  residue_atoms = residue.atoms()
  atom_names = residue_atoms.extract_name()
  matched_atom_names = iotbx.pdb.atom_name_interpretation.interpreters[
    residue.resname].match_atom_names(atom_names=atom_names)
  mon_lib_atom_names = matched_atom_names.mon_lib_names()
  #
  rotamer_info = comp_comp_id.rotamer_info()
  bonds_to_omit = mmtbx.monomer_library.rotamer_utils.extract_bonds_to_omit(
    rotamer_info = rotamer_info)
  #
  external_edge_list = []
  secial_cases = {
    "ASP": [(["OD1"],["HD1","1HD","DD1","1DD"]),
            (["OD2"],["HD2","2HD","DD2","2DD"])],
    "GLU": [(["OE1"],["HE1","1HE","DE1","1DE"]),
            (["OE2"],["HE2","2HE","DE2","2DE"])]
  }
  if(residue.resname.strip().upper() in secial_cases.keys()):
    edge = []
    bonded_atom_names = secial_cases[residue.resname.strip().upper()]
    for ban in bonded_atom_names:
      for atom1 in residue_atoms:
        if(atom1.name.strip().upper() in ban[0]):
          for atom2 in residue_atoms:
            if(atom2.name.strip().upper() in ban[1]):
              assert atom1.i_seq != atom2.i_seq
              edge = [atom1.i_seq, atom2.i_seq]
      if(len(edge) > 0): external_edge_list.append(edge)
  #
  tardy_model = mmtbx.monomer_library.rotamer_utils.tardy_model(
    comp_comp_id       = comp_comp_id,
    input_atom_names   = atom_names,
    mon_lib_atom_names = mon_lib_atom_names,
    sites_cart         = residue_atoms.extract_xyz(),
    bonds_to_omit      = bonds_to_omit,
    constrain_dihedrals_with_sigma_less_than_or_equal_to = None,
    external_edge_list = external_edge_list,
    external_clusters  = external_clusters,
    return_none_if_unexpected_degrees_of_freedom = True)
  if(tardy_model is None):
    mes = "TARDY: cannot create tardy model for: %s (corrupted residue). Skipping it."
    print >> log, mes%residue.id_str(suppress_segid=1)[-12:]
    return None
  joint_dofs = tardy_model.degrees_of_freedom_each_joint()
  if(joint_dofs[0] != 0 or not joint_dofs[1:].all_eq(1)):
    mes = "TARDY error: unexpected degrees of freedom for %s. Skipping it."
    print >> log, mes%residue.id_str(suppress_segid=1)[-12:]
    return None
  return tardy_model

def axes_and_atoms_aa_specific(residue, mon_lib_srv,
                               remove_clusters_with_all_h=False, log=None):
  get_class = iotbx.pdb.common_residue_names_get_class
  if 0: print residue.id_str(suppress_segid=1)[-12:]
  if(not (get_class(residue.resname) == "common_amino_acid")): return None
  tardy_model = tardy_model_one_residue(residue = residue,
    mon_lib_srv = mon_lib_srv, log = log)
  if(tardy_model is None): return None
  clusters = tardy_model.tardy_tree.cluster_manager.clusters[1:]
  axes = tardy_model.tardy_tree.cluster_manager.hinge_edges[1:]
  assert len(clusters) == len(axes)
  if(len(axes)==0): return None
  if 0:
    print "clusters:", clusters
    print "    axes:", axes
    print
  #
  ic = 0
  axes_and_atoms_to_rotate = []
  while ic < len(axes):
    axis = axes[ic]
    cluster = clusters[ic]
    n_branches = 0
    for i_axis_, axis_ in enumerate(axes):
      if(i_axis_ > ic and axis[0] == axis_[0]): n_branches += 1
    if(n_branches == 0):
      atoms_to_rotate = []
      for ci in clusters[ic:]:
        atoms_to_rotate.extend(ci)
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
    else:
      assert n_branches <= 2, n_branches
      # branch 1
      icb = ic+1
      atoms_to_rotate = cluster
      while icb < len(axes):
        if(axis[1]==axes[icb][0]): atoms_to_rotate.extend(clusters[icb])
        icb += 1
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
      # branch 2
      axis = axes[ic]
      cluster = clusters[ic]
      icb = ic+1
      atoms_to_rotate = cluster
      while icb < len(axes):
        if(axis[0]==axes[icb][0]): atoms_to_rotate.extend(clusters[icb])
        icb += 1
      axes_and_atoms_to_rotate.append([axis, atoms_to_rotate])
      ic += 1
  #
  if 0:
    print axes_and_atoms_to_rotate
    print
  if(remove_clusters_with_all_h):
    tmp = []
    residue_atoms = residue.atoms()
    for axis, cluster in axes_and_atoms_to_rotate:
      count_h = 0
      for c_i in cluster:
        if(residue_atoms[c_i].element.strip().upper() in ["H","D"]): count_h +=1
      if(count_h != len(cluster)):
        tmp.append([axis, cluster])
    axes_and_atoms_to_rotate = tmp
  return axes_and_atoms_to_rotate

class residue_rsr_monitor(object):
  def __init__(self,
               residue_id_str = None,
               selection = None,
               sites_cart = None,
               twomfodfc = None,
               mfodfc = None,
               cc = None,
               residue_type = None,
               validate_iselection = None):
    adopt_init_args(self, locals())

class select_map(object):
  def __init__(self, unit_cell, target_map_data, model_map_data):
    adopt_init_args(self, locals())
    assert target_map_data.focus() == model_map_data.focus()
    assert target_map_data.all() == model_map_data.all()
    self.fft_n_real = target_map_data.focus()
    self.fft_m_real = target_map_data.all()

  def initialize_rotamers (self) :
    # XXX initialize classes needed for rotamer check
    self.sa, self.r, self.rot = [None]*3
    from mmtbx.rotamer.sidechain_angles import SidechainAngles
    from mmtbx.rotamer import rotamer_eval
    from mmtbx.validation.rotalyze import rotalyze
    self.sa = SidechainAngles(False)
    self.r = rotamer_eval.RotamerEval()
    self.rot = rotalyze()

  def select(self, sites_cart, atom_radius=2.0):
    return maptbx.grid_indices_around_sites(
      unit_cell  = self.unit_cell,
      fft_n_real = self.fft_n_real,
      fft_m_real = self.fft_m_real,
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), atom_radius))

  def get_cc(self, sites_cart, residue_iselection):
    sel_map = self.select(sites_cart = sites_cart)
    m1 = self.target_map_data.select(sel_map)
    m2 = self.model_map_data.select(sel_map)
    return flex.linear_correlation(x = m1, y = m2).coefficient()

  def get_map_sum(self, sites_cart, residue_iselection, map):
    sel_map = self.select(sites_cart = sites_cart, atom_radius=0.3)
    return flex.sum(map.select(sel_map))

  def is_refinement_needed(self, residue_group, residue, cc_limit, ignore_hd):
    result = False
    if([self.sa, self.r, self.rot].count(None)==0):
      all_dict = self.rot.construct_complete_sidechain(
        residue_group = residue_group)
      is_outlier, value = self.rot.evaluate_residue(residue_group, self.sa,
        self.r, all_dict)
      if(is_outlier): return True
    for atom in residue.atoms():
      if(not atom.element.strip().lower() in ["h","d"]):
        sel_map = self.select(sites_cart = flex.vec3_double([atom.xyz]))
        m1 = self.target_map_data.select(sel_map)
        m2 = self.model_map_data.select(sel_map)
        cc = flex.linear_correlation(x = m1, y = m2).coefficient()
        if(cc < cc_limit): return True
    return result


class refiner(object):
  def __init__(self,
               pdb_hierarchy,
               target_map,
               geometry_restraints_manager,
               real_space_target_weight,
               real_space_gradients_delta,
               max_iterations):
    self.target_map = target_map
    self.real_space_target_weight = real_space_target_weight
    self.real_space_gradients_delta = real_space_gradients_delta
    self.geometry_restraints_manager = geometry_restraints_manager
    self.pdb_hierarchy = pdb_hierarchy
    #self.pdb_atoms = pdb_hierarchy.atoms()
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)

  def refine_restrained(self, sites_cart_rsel, rsel, rs, use_lockit=False):
    assert rsel.size() == self.pdb_hierarchy.atoms_size()
    assert sites_cart_rsel.size() == rsel.count(True)
    geometry_restraints_manager = \
      self.geometry_restraints_manager.select(rsel)
    if(use_lockit): # XXX trying lockit...
      work_params = lockit.get_master_phil().extract()
      work_params.coordinate_refinement.lbfgs_max_iterations=50
      work_params.coordinate_refinement.compute_final_correlation=False
      work_params.coordinate_refinement.finishing_geometry_minimization.cycles_max=50
      result = lockit.run_coordinate_refinement(
        sites_cart                  = sites_cart_rsel,
        geometry_restraints_manager = geometry_restraints_manager,
        selection_variable          = rs,
        density_map                 = self.target_map,
        real_space_gradients_delta  = self.real_space_gradients_delta,
        work_params                 = work_params,
        home_restraints_list=[],
        work_scatterers=None,
        unit_cell=None,
        d_min=None,
        write_pdb_callback=None,
        log=None)
      return result.sites_cart.select(rs)
    result = maptbx.real_space_refinement_simple.lbfgs(
      selection_variable              = rs,
      sites_cart                      = sites_cart_rsel,
      density_map                     = self.target_map,
      geometry_restraints_manager     = geometry_restraints_manager,
      real_space_target_weight        = self.real_space_target_weight,
      real_space_gradients_delta      = self.real_space_gradients_delta,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)
    return result.sites_cart_variable

def target(sites_cart_residue, unit_cell, m):
  sites_frac_residue = unit_cell.fractionalize(sites_cart_residue)
  result = 0
  for rsf in sites_frac_residue:
    result += m.eight_point_interpolation(rsf)
  return result


class rotamer_evaluator(object):
  def __init__(self, sites_cart_start,
                     unit_cell,
                     two_mfo_dfc_map,
                     mfo_dfc_map):
    adopt_init_args(self, locals())
    t1 = target(self.sites_cart_start, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(self.sites_cart_start, self.unit_cell, self.mfo_dfc_map)
    self.t1 = t1
    self.t2 = t2
    self.t_start = t1+t2
    self.t_best = self.t_start
    self.t1_start = t1
    self.t1_best = self.t1_start
    self.t2_start = t2
    self.t2_best = self.t2_start

  def is_better(self, sites_cart):
    t1 = target(sites_cart, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(sites_cart, self.unit_cell, self.mfo_dfc_map)
    t = t1+t2#*3 # XXX very promising thing to do, but reaaly depends on resolution
    result = False
    size = sites_cart.size()
    if 1:
      if(t > self.t_best):
        if((t2 > 0 and self.t2_best > 0 and t2 > self.t2_best) or
           (t2 < 0 and self.t2_best < 0 and abs(t2)<abs(self.t2_best)) or
           (t2 > 0 and self.t2_best < 0)):
          result = True
          self.t2_best = t2
          self.t1_best = t1
          self.t_best = t
    return result

def include_residue_selection(selection, residue_iselection):
  size = selection.size()
  selection_ = selection.deep_copy()
  selection_ = selection_.iselection()
  selection_.extend(residue_iselection)
  new_sel = flex.bool(size, selection_)
  rs = flex.bool(size, residue_iselection).select(new_sel).iselection()
  return new_sel, rs

def generate_range(start, stop, step):
  assert abs(start) <= abs(stop)
  inc = start
  result = []
  while abs(inc) <= abs(stop):
    result.append(inc)
    inc += step
  return result

def rotate_point_around_axis(axis_point_1, axis_point_2, point, angle_deg):
  angle_rad = angle_deg*math.pi/180.
  xa,ya,za = axis_point_1
  xb,yb,zb = axis_point_2
  x,y,z = point
  xL,yL,zL = xb-xa,yb-ya,zb-za
  xLsq = xL**2
  yLsq = yL**2
  zLsq = zL**2
  dLsq = xLsq + yLsq + zLsq
  dL = math.sqrt(dLsq)
  ca = math.cos(angle_rad)
  dsa = math.sin(angle_rad)/dL
  oca = (1-ca)/dLsq
  xLyLo = xL*yL*oca
  xLzLo = xL*zL*oca
  yLzLo = yL*zL*oca
  xma,yma,zma = x-xa,y-ya,z-za
  m1 = xLsq*oca+ca
  m2 = xLyLo-zL*dsa
  m3 = xLzLo+yL*dsa
  m4 = xLyLo+zL*dsa
  m5 = yLsq*oca+ca
  m6 = yLzLo-xL*dsa
  m7 = xLzLo-yL*dsa
  m8 = yLzLo+xL*dsa
  m9 = zLsq*oca+ca
  x_new = xma*m1 + yma*m2 + zma*m3 + xa
  y_new = xma*m4 + yma*m5 + zma*m6 + ya
  z_new = xma*m7 + yma*m8 + zma*m9 + za
  return (x_new,y_new,z_new)

def torsion_search(residue_evaluator,
                   cluster_evaluators,
                   axes_and_atoms_to_rotate,
                   rotamer_sites_cart,
                   rotamer_id_best,
                   residue_sites_best,
                   params = None,
                   rotamer_id = None):
  if(params is None):
    params = torsion_search_params().extract().torsion_search
    params.range_start = 0
    params.range_stop = 360
    params.step = 1.0
  rotamer_sites_cart_ = rotamer_sites_cart.deep_copy()
  n_clusters = len(axes_and_atoms_to_rotate)
  c_counter = 0
  for cluster_evaluator, aa in zip(cluster_evaluators,axes_and_atoms_to_rotate):
    c_counter += 1
    axis = aa[0]
    atoms = aa[1]
    angle_deg_best = None
    angle_deg_good = None
    for angle_deg in generate_range(start = params.range_start, stop =
                                    params.range_stop, step = params.step):
      if(c_counter != n_clusters):
        new_xyz = flex.vec3_double([rotate_point_around_axis(
          axis_point_1 = rotamer_sites_cart[axis[0]],
          axis_point_2 = rotamer_sites_cart[axis[1]],
          point  = rotamer_sites_cart[atoms[0]],
          angle_deg = angle_deg)])
      else:
        new_xyz = flex.vec3_double()
        for atom in atoms:
          new_xyz.append(rotate_point_around_axis(
            axis_point_1 = rotamer_sites_cart[axis[0]],
            axis_point_2 = rotamer_sites_cart[axis[1]],
            point  = rotamer_sites_cart[atom],
            angle_deg = angle_deg))
      if(cluster_evaluator.is_better(sites_cart = new_xyz)):
        if(angle_deg_best is not None and
           abs(abs(angle_deg_best)-abs(angle_deg))>
           params.min_angle_between_solutions):
          angle_deg_good = angle_deg_best
        angle_deg_best = angle_deg
    if(angle_deg_best is not None):
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = rotamer_sites_cart[axis[0]],
          axis_point_2 = rotamer_sites_cart[axis[1]],
          point  = rotamer_sites_cart[atom],
          angle_deg = angle_deg_best)
        rotamer_sites_cart[atom] = new_xyz
    if(angle_deg_good is not None):
      for atom in atoms:
        new_xyz = rotate_point_around_axis(
          axis_point_1 = rotamer_sites_cart_[axis[0]],
          axis_point_2 = rotamer_sites_cart_[axis[1]],
          point  = rotamer_sites_cart_[atom],
          angle_deg = angle_deg_best)
        rotamer_sites_cart_[atom] = new_xyz
  for rsc in [rotamer_sites_cart, rotamer_sites_cart_]:
    if(residue_evaluator.is_better(sites_cart = rsc)):
      rotamer_id_best = rotamer_id
      residue_sites_best = rsc.deep_copy()
  return residue_sites_best, rotamer_id_best


def residue_iteration(pdb_hierarchy,
                      xray_structure,
                      selection,
                      target_map_data,
                      model_map_data,
                      residual_map_data,
                      mon_lib_srv,
                      rsr_manager,
                      optimize_hd,
                      params,
                      log):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  assert target_map_data.focus() == model_map_data.focus()
  assert target_map_data.all() == model_map_data.all()
  fmt1 = "                |--------START--------| |-----FINAL----|"
  fmt2 = "     residue   map_cc 2mFo-DFc mFo-DFc 2mFo-DFc mFo-DFc" \
    " rotamer n_rot max_moved"
  fmt3 = "  %12s%7.4f %8.2f %7.2f %8.2f %7.2f %7s %5d  %8.3f"
  print >> log, fmt1
  print >> log, fmt2
  unit_cell = xray_structure.unit_cell()
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    target_map_data = target_map_data,
    model_map_data = model_map_data)
  map_selector.initialize_rotamers()
  get_class = iotbx.pdb.common_residue_names_get_class
  n_other_residues = 0
  n_amino_acids_ignored = 0
  n_amino_acids_scored = 0
  sites_cart_start = xray_structure.sites_cart()
  result = []
  for model in pdb_hierarchy.models():
    for chain in model.chains():
      for residue_group in chain.residue_groups():
        conformers = residue_group.conformers()
        if(params.ignore_alt_conformers and len(conformers)>1): continue
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          if(get_class(residue.resname) == "common_amino_acid"):
            residue_iselection = residue.atoms().extract_i_seq()
            sites_cart_residue = xray_structure.sites_cart().select(residue_iselection)
            residue.atoms().set_xyz(new_xyz=sites_cart_residue)
            max_moved_dist = 0
            sites_cart_residue_start = sites_cart_residue.deep_copy()
            # XXX assume that "atoms" are the same in residue and residue_groups
            if(map_selector.is_refinement_needed(
               residue_group = residue_group,
               residue       = residue,
               cc_limit      = params.poor_cc_threshold,
               ignore_hd     = optimize_hd)):
              residue_id_str = residue.id_str(suppress_segid=1)[-12:]
              rsel, rs = include_residue_selection(
                selection          = selection,
                residue_iselection = residue_iselection)
              cc_start = map_selector.get_cc(
                sites_cart         = sites_cart_residue,
                residue_iselection = residue_iselection)
              rotamer_id_best = None
              rev = rotamer_evaluator(
                sites_cart_start = sites_cart_residue,
                unit_cell        = unit_cell,
                two_mfo_dfc_map  = target_map_data,
                mfo_dfc_map      = residual_map_data)
              residue_sites_best = sites_cart_residue.deep_copy()
              rm = residue_rsr_monitor(
                residue_id_str = residue_id_str,
                selection  = residue_iselection.deep_copy(),
                sites_cart = sites_cart_residue.deep_copy(),
                twomfodfc  = rev.t1_start,
                mfodfc     = rev.t2_start,
                cc         = cc_start)
              result.append(rm)
              axes_and_atoms_to_rotate = axes_and_atoms_aa_specific(
                residue     = residue,
                mon_lib_srv = mon_lib_srv,
                remove_clusters_with_all_h = optimize_hd,
                log         = log)
              if(axes_and_atoms_to_rotate is not None and
                 len(axes_and_atoms_to_rotate) > 0):
                # initialize criteria for first rotatable atom in each cluster
                rev_first_atoms = []
                for i_aa, aa in enumerate(axes_and_atoms_to_rotate):
                  if(i_aa == len(axes_and_atoms_to_rotate)-1):
                    sites_aa = flex.vec3_double()
                    for aa_ in aa[1]:
                      sites_aa.append(sites_cart_residue[aa_])
                  else:
                    sites_aa = flex.vec3_double([sites_cart_residue[aa[1][0]]])
                  rev_i = rotamer_evaluator(
                    sites_cart_start = sites_aa,
                    unit_cell        = unit_cell,
                    two_mfo_dfc_map  = target_map_data,
                    mfo_dfc_map      = residual_map_data)
                  rev_first_atoms.append(rev_i)
                # get rotamer iterator
                rotamer_iterator = lockit.get_rotamer_iterator(
                  mon_lib_srv         = mon_lib_srv,
                  residue             = residue,
                  atom_selection_bool = None)
                if(rotamer_iterator is None):
                  n_amino_acids_ignored += 1
                  n_rotamers = 0
                  print >> log, "No rotamers for: %s. Use torsion grid search."%\
                    residue_id_str
                  residue_sites_best, rotamer_id_best = torsion_search(
                    residue_evaluator        = rev,
                    cluster_evaluators       = rev_first_atoms,
                    axes_and_atoms_to_rotate = axes_and_atoms_to_rotate,
                    rotamer_sites_cart       = sites_cart_residue,
                    rotamer_id_best          = rotamer_id_best,
                    residue_sites_best       = residue_sites_best,
                    rotamer_id               = None,
                    params                   = None)
                else:
                  n_amino_acids_scored += 1
                  n_rotamers = 0
                  if(not params.use_rotamer_iterator):
                    if(params.torsion_grid_search):
                      residue_sites_best, rotamer_id_best = torsion_search(
                        residue_evaluator        = rev,
                        cluster_evaluators       = rev_first_atoms,
                        axes_and_atoms_to_rotate = axes_and_atoms_to_rotate,
                        rotamer_sites_cart       = sites_cart_residue,
                        rotamer_id_best          = rotamer_id_best,
                        residue_sites_best       = residue_sites_best,
                        rotamer_id               = None,
                        params                   = params.torsion_search)
                  else:
                    for rotamer, rotamer_sites_cart in rotamer_iterator:
                      n_rotamers += 1
                      if(params.torsion_grid_search):
                        residue_sites_best, rotamer_id_best = torsion_search(
                          residue_evaluator        = rev,
                          cluster_evaluators       = rev_first_atoms,
                          axes_and_atoms_to_rotate = axes_and_atoms_to_rotate,
                          rotamer_sites_cart       = rotamer_sites_cart,
                          rotamer_id_best          = rotamer_id_best,
                          residue_sites_best       = residue_sites_best,
                          rotamer_id               = rotamer.id,
                          params                   = params.torsion_search)
                      else:
                        if(rev.is_better(sites_cart = rotamer_sites_cart)):
                          rotamer_id_best = rotamer.id
                          residue_sites_best = rotamer_sites_cart.deep_copy()
                  residue.atoms().set_xyz(new_xyz=residue_sites_best)
                max_moved_dist = flex.max(flex.sqrt(
                  (sites_cart_residue_start-residue_sites_best).dot()))
                if(not params.real_space_refine_rotamer):
                  sites_cart_start = sites_cart_start.set_selected(
                    residue_iselection, residue_sites_best)
                else:
                  tmp = sites_cart_start.set_selected(
                    residue_iselection, residue_sites_best)
                  sites_cart_refined = rsr_manager.refine_restrained(
                    tmp.select(rsel), rsel, rs)
                  if(rev.is_better(sites_cart = sites_cart_refined)):
                    sites_cart_start = sites_cart_start.set_selected(
                      residue_iselection, sites_cart_refined)
                    residue.atoms().set_xyz(new_xyz=sites_cart_refined)
                    max_moved_dist = flex.max(flex.sqrt(
                      (sites_cart_residue_start - sites_cart_refined).dot()))
              if(abs(rev.t1_best-rev.t1_start) > 0.01 and
                 abs(rev.t2_best-rev.t2_start) > 0.01):
                print >> log, fmt3 % (
                  residue_id_str, cc_start, rev.t1_start, rev.t2_start,
                  rev.t1_best, rev.t2_best, rotamer_id_best, n_rotamers,
                  max_moved_dist)
  xray_structure.set_sites_cart(sites_cart_start)
  return result

def get_map_data(fmodel, map_type, resolution_factor=1./4, kick=False):
  if(kick):
    km = map_tools.kick_map(
      fmodel            = fmodel,
      map_type          = map_type,
      kick_sizes        = [0.0,0.3,0.5],
      number_of_kicks   = 50,
      real_map          = False,
      real_map_unpadded = True,
      update_bulk_solvent_and_scale = False,
      symmetry_flags    = maptbx.use_space_group_symmetry,
      average_maps      = False)
    map_data = km.map_data
    fft_map = km.fft_map
  else:
    map_obj = fmodel.electron_density_map()
    fft_map = map_obj.fft_map(resolution_factor = resolution_factor,
      map_type = map_type)
    fft_map.apply_sigma_scaling()
    map_data = fft_map.real_map_unpadded()
  return map_data,fft_map

def validate_changes(fmodel, residue_rsr_monitor, validate_method, log):
  assert (validate_method is None) or (hasattr(validate_method, "__call__"))
  xray_structure = fmodel.xray_structure
  target_map_data,fft_map_1 = get_map_data(
    fmodel = fmodel, map_type = "2mFo-DFc", kick=False)
  model_map_data,fft_map_2 = get_map_data(
    fmodel = fmodel, map_type = "Fc")
  residual_map_data,fft_map_3 = get_map_data(
    fmodel = fmodel, map_type = "mFo-DFc", kick=False)
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    target_map_data = target_map_data,
    model_map_data = model_map_data)
  map_selector.initialize_rotamers()
  sites_cart = xray_structure.sites_cart()
  sites_cart_result = sites_cart.deep_copy()
  unit_cell = xray_structure.unit_cell()
  print >> log, "Validate:"
  fmt1 = "                |--------START--------|  |--------FINAL--------|"
  fmt2 = "     residue    map_cc 2mFo-DFc mFo-DFc  map_cc 2mFo-DFc mFo-DFc"
  fmt3 = "  %12s %7.4f %8.2f %7.2f %7.4f %8.2f %7.2f"
  print >> log, fmt1
  print >> log, fmt2
  get_class = iotbx.pdb.common_residue_names_get_class
  for rm in residue_rsr_monitor:
    residue_name = rm.residue_id_str[1:4]
    if(get_class(name=residue_name) == "common_amino_acid"):
      sites_cart_residue = sites_cart.select(rm.selection)
      t1 = target(sites_cart_residue, unit_cell, target_map_data)
      t2 = target(sites_cart_residue, unit_cell, residual_map_data)
      cc = map_selector.get_cc(sites_cart = sites_cart_residue,
        residue_iselection = rm.selection)
      flag = ""
      ###
      if(abs(t2)<2.0): t2 = 0.
      if(abs(rm.mfodfc)<2.0): rm.mfodfc = 0.
      ###
      dmif1 = rm.mfodfc < 0 and t2 < 0 and t2 < rm.mfodfc
      dmif2 = rm.mfodfc > 0 and t2 < 0 and abs(t2) > abs(rm.mfodfc)
      dmif4 = rm.mfodfc<0 and t2<0 and (abs(rm.mfodfc)>5 and 2*abs(rm.mfodfc)<abs(t2))
      dmif41 = rm.mfodfc<0 and t2<0 and (abs(rm.mfodfc)>2 and 3*abs(rm.mfodfc)<abs(t2))
      dmif5 = abs(rm.mfodfc)<0.5 and t2<-5.
      dmif6 = rm.cc > cc and abs(rm.mfodfc)<0.5 and t2 < -5.
      dmif = dmif1 or dmif2
      if((cc < rm.cc or t1 < rm.twomfodfc) and
         (dmif or dmif4 or dmif5 or dmif6 or dmif41)):
        flag = " <<<"
      elif (validate_method is not None) :
        is_outlier = validate_method(
          res_type=rm.residue_type,
          i_seqs=rm.validate_iselection,
          sites_cart=sites_cart)
        if (is_outlier) :
          flag = " <<<"
      print >> log, fmt3 % (
        rm.residue_id_str, rm.cc, rm.twomfodfc, rm.mfodfc, cc,t1,t2), flag
      if(len(flag)>0):
        sites_cart_result = sites_cart_result.set_selected(
          rm.selection, rm.sites_cart)
  xray_structure.set_sites_cart(sites_cart_result)
  fmodel.update_xray_structure(xray_structure = xray_structure,
    update_f_calc=True, update_f_mask=True)
  print >> log, "r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free())

def rsr_overall (model,
                 fmodel,
                 params,
                 optimize_hd,
                 macro_cycle,
                 log) :
  assert model.xray_structure is fmodel.xray_structure
  rsr_params = mmtbx.refinement.real_space.master_params().extract()
  rsr_params.real_space_refinement.mode="diff_map"
  if(params.exclude_hydrogens and optimize_hd):
    hd_selection = fmodel.xray_structure.hd_selection()
    occupancies_cache= fmodel.xray_structure.scatterers().extract_occupancies()
    fmodel.xray_structure.set_occupancies(value=0, selection=hd_selection)
  mmtbx.refinement.real_space.run(
    fmodel = fmodel, # XXX neutron ?
    model  = model,
    params = rsr_params.real_space_refinement,
    log    = log)
  if(params.exclude_hydrogens and optimize_hd):
    fmodel.xray_structure.set_occupancies(value = occupancies_cache)
  assert model.xray_structure is fmodel.xray_structure
  fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)

def run(fmodel,
        geometry_restraints_manager,
        pdb_hierarchy,
        solvent_selection,
        optimize_hd = False,
        log = None,
        params = None):
  if(log is None): log = sys.stdout
  if(params is None): params = master_params().extract()
  print_statistics.make_sub_header(text="Fitting sidechains", out=log)
  mon_lib_srv = mmtbx.monomer_library.server.server()
  if(params.use_dihedral_restraints):
    sel = flex.bool(fmodel.xray_structure.scatterers().size(), True)
    geometry_restraints_manager = geometry_restraints_manager.select(sel)
    geometry_restraints_manager.remove_dihedrals_in_place(sel)
  restraints_manager = mmtbx.restraints.manager(
    geometry      = geometry_restraints_manager,
    normalization = True)
  model = mmtbx.model.manager(
    restraints_manager = restraints_manager,
    xray_structure = fmodel.xray_structure,
    pdb_hierarchy = pdb_hierarchy)
  backbone_selections = model.backbone_selections()
  if(params.ignore_water_when_move_sidechains):
    selection = ~model.solvent_selection()
  else:
    selection = flex.bool(model.xray_structure.scatterers().size(), True)
  selection &= backbone_selections
  fmt = "Macro-cycle %2d: r_work=%6.4f r_free=%6.4f"
  print >> log, fmt%(0, fmodel.r_work(), fmodel.r_free())
  for macro_cycle in range(1,params.number_of_macro_cycles+1):
    target_map_data, fft_map_1 = get_map_data(fmodel = fmodel,
      map_type = params.target_map, kick=False)
    model_map_data,fft_map_2 = get_map_data(fmodel = fmodel,
      map_type = params.model_map)
    residual_map_data,fft_map_3 = get_map_data(fmodel = fmodel,
      map_type = params.residual_map, kick=False)
    if(params.filter_residual_map_value is not None): #XXX use filtering....
      map_sel = flex.abs(residual_map_data) < params.filter_residual_map_value
      residual_map_data = residual_map_data.set_selected(map_sel, 0)
    if(params.filter_2fofc_map is not None):
      map_sel = flex.abs(target_map_data) < params.filter_2fofc_map
      target_map_data = target_map_data.set_selected(map_sel, 0)
    rsr_manager = refiner(
      pdb_hierarchy               = pdb_hierarchy,
      target_map                  = target_map_data,
      geometry_restraints_manager = geometry_restraints_manager,
      real_space_target_weight    = params.residue_iteration.real_space_refine_target_weight,
      real_space_gradients_delta  = fmodel.f_obs().d_min()/4,
      max_iterations              = params.residue_iteration.real_space_refine_max_iterations)
    residue_rsr_monitor = residue_iteration(
      pdb_hierarchy     = pdb_hierarchy,
      xray_structure    = fmodel.xray_structure,
      selection         = selection,
      target_map_data   = target_map_data,
      model_map_data    = model_map_data,
      residual_map_data = residual_map_data,
      mon_lib_srv       = mon_lib_srv,
      rsr_manager       = rsr_manager,
      optimize_hd       = optimize_hd,
      params            = params.residue_iteration,
      log               = log)
    fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
    print >> log, "1:", fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    del target_map_data, model_map_data, residual_map_data, fft_map_1, fft_map_2, fft_map_3
    if(params.real_space_refine_overall):
      rsr_overall(
        model=model,
        fmodel=fmodel,
        params=params,
        optimize_hd=optimize_hd,
        macro_cycle=macro_cycle,
        log=log)
      print >> log, "2:",fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    if(params.validate_change):
      validate_changes(fmodel = fmodel,
                       residue_rsr_monitor = residue_rsr_monitor,
                       validate_method=None,
                       log = log)

from cctbx.array_family import flex
from libtbx import adopt_init_args
from mmtbx.command_line import lockit
from mmtbx import real_space_correlation
from libtbx import adopt_init_args
from cctbx import maptbx
import scitbx.lbfgs
import iotbx.pdb
import mmtbx.monomer_library
import mmtbx.model
import mmtbx.refinement.real_space
from cctbx import miller
from mmtbx import map_tools
from libtbx.test_utils import approx_equal, not_approx_equal
from mmtbx import masks
from cctbx import maptbx
import time


class residue_rsr_monitor(object):
  def __init__(self,
               residue_id_str = None,
               selection = None,
               sites_cart = None,
               twomfodfc = None,
               mfodfc = None,
               cc = None):
    adopt_init_args(self, locals())

class select_map(object):
  def __init__(self, unit_cell, map_data_1, map_data_2):
    adopt_init_args(self, locals())
    assert map_data_1.focus() == map_data_2.focus()
    assert map_data_1.all() == map_data_2.all()
    self.fft_n_real = map_data_1.focus()
    self.fft_m_real = map_data_1.all()
    # XXX initialize classes needed for rotamer check
    # XXX Jeff/Nat: clean it so it is availabe in mmtbx and not from phenix
    self.sa, self.r, self.rot = [None]*3
    try:
      from mmtbx.rotamer.sidechain_angles import SidechainAngles
      from mmtbx.rotamer import rotamer_eval
      from phenix.validation.rotalyze import rotalyze
      self.sa = SidechainAngles(False)
      self.r = rotamer_eval.RotamerEval()
      self.rot = rotalyze()
    except ImportError: pass

  def select(self, sites_cart, atom_radius=2.0):
    return maptbx.grid_indices_around_sites(
      unit_cell  = self.unit_cell,
      fft_n_real = self.fft_n_real,
      fft_m_real = self.fft_m_real,
      sites_cart = sites_cart,
      site_radii = flex.double(sites_cart.size(), atom_radius))

  def get_cc(self, sites_cart, residue_iselection):
    sel_map = self.select(sites_cart = sites_cart)
    m1 = self.map_data_1.select(sel_map)
    m2 = self.map_data_2.select(sel_map)
    return flex.linear_correlation(x = m1, y = m2).coefficient()

  def get_map_sum(self, sites_cart, residue_iselection, map):
    sel_map = self.select(sites_cart = sites_cart, atom_radius=0.3)
    return flex.sum(map.select(sel_map))

  def is_refinement_needed(self, residue_group, residue, cc_limit = 0.9):
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
        m1 = self.map_data_1.select(sel_map)
        m2 = self.map_data_2.select(sel_map)
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
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)

  def refine_restrained(self, sites_cart_rsel, rsel, rs):
    assert rsel.size() == self.pdb_hierarchy.atoms_size()
    assert sites_cart_rsel.size() == rsel.count(True)
    geometry_restraints_manager = \
      self.geometry_restraints_manager.select(rsel)
    result = maptbx.real_space_refinement_simple.lbfgs(
      selection_variable              = rs,
      sites_cart                      = sites_cart_rsel,
      density_map                     = self.target_map,
      geometry_restraints_manager     = geometry_restraints_manager,
      real_space_target_weight        = self.real_space_target_weight,
      real_space_gradients_delta      = self.real_space_gradients_delta,
      lbfgs_termination_params        = self.lbfgs_termination_params,
      lbfgs_exception_handling_params = self.lbfgs_exception_handling_params)
    return result

def target(sites_cart_residue, unit_cell, m):
  sites_frac_residue = unit_cell.fractionalize(sites_cart_residue)
  result = 0
  for rsf in sites_frac_residue:
    result += m.eight_point_interpolation(rsf)
  return result

class rotomer_evaluator(object):
  def __init__(self, sites_cart_start,
                     unit_cell,
                     two_mfo_dfc_map,
                     mfo_dfc_map):
    adopt_init_args(self, locals())
    t1 = target(self.sites_cart_start, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(self.sites_cart_start, self.unit_cell, self.mfo_dfc_map)
    self.t_start = t1+t2
    self.t_best = self.t_start
    self.t1_start = t1
    self.t1_best = self.t1_start
    self.t2_start = t2
    self.t2_best = self.t2_start

  def is_better(self, sites_cart):
    t1 = target(sites_cart, self.unit_cell, self.two_mfo_dfc_map)
    t2 = target(sites_cart, self.unit_cell, self.mfo_dfc_map)
    t = t1+t2
    result = False
    if(t > self.t_best):
      #self.t_best = t # this gives lower Rfactor
      if(t2 > 0 and self.t2_best > 0 and t2 > self.t2_best):
        result = True
        self.t2_best = t2
        self.t1_best = t1
        self.t_best = t
      elif(t2 < 0 and self.t2_best < 0 and abs(t2)<abs(self.t2_best)):
        result = True
        self.t2_best = t2
        self.t1_best = t1
        self.t_best = t
      elif(t2 > 0 and self.t2_best < 0):#XXX and abs(t2)>abs(self.t2_best)):
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

def iterate_rotamers(pdb_hierarchy,
                     xray_structure,
                     selection,
                     map_data_1,
                     map_data_2,
                     map_data_3,
                     mon_lib_srv,
                     rsr_manager,
                     f_obs, fft_map,
                     poor_cc_threshold,
                     log,
                     real_space_refine = "restrained",
                     ignore_alt_conformers = True):
  assert map_data_1.focus() == map_data_2.focus()
  assert map_data_1.all() == map_data_2.all()
  assert real_space_refine in ["restrained","constrained"]
  fmt1 = "                |--------START--------| |-----FINAL----|"
  fmt2 = "     residue    map_cc 2mFo-DFc mFo-DFc 2mFo-DFc mFo-DFc" \
    " rotomer n_rotamers"
  fmt3 = "  %12s %7.4f %8.2f %7.2f %8.2f %7.2f %7s %10d"
  print >> log, fmt1
  print >> log, fmt2
  unit_cell = xray_structure.unit_cell()
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    map_data_1 = map_data_1,
    map_data_2 = map_data_2)
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
        if(ignore_alt_conformers and len(conformers)>1): continue
        for conformer in residue_group.conformers():
          residue = conformer.only_residue()
          # XXX assume that "atoms" are the same in residue and residue_groups
          if(map_selector.is_refinement_needed(residue_group = residue_group,
            residue = residue, cc_limit = poor_cc_threshold)):
            residue_id_str = residue.id_str(suppress_segid=1)[-12:]
            residue_iselection = residue.atoms().extract_i_seq()
            rsel, rs = include_residue_selection(
              selection          = selection,
              residue_iselection = residue_iselection)
            sites_cart = xray_structure.sites_cart().select(residue_iselection)
            cc_start = map_selector.get_cc(
              sites_cart         = sites_cart,
              residue_iselection = residue_iselection)
            rotamer_id_best = None
            rev = rotomer_evaluator(
              sites_cart_start = sites_cart,
              unit_cell        = unit_cell,
              two_mfo_dfc_map  = map_data_1,
              mfo_dfc_map      = map_data_3)
            residue_sites_best = sites_cart.deep_copy()
            rm = residue_rsr_monitor(
              residue_id_str = residue_id_str,
              selection  = residue_iselection.deep_copy(),
              sites_cart = sites_cart.deep_copy(),
              twomfodfc  = rev.t1_start,
              mfodfc     = rev.t2_start,
              cc         = cc_start)
            result.append(rm)
            if(get_class(residue.resname) != "common_amino_acid"):
              n_other_residues += 1
            else:
              rotamer_iterator = lockit.get_rotamer_iterator(
                mon_lib_srv         = mon_lib_srv,
                residue             = residue,
                atom_selection_bool = None)
              if(rotamer_iterator is None):
                n_amino_acids_ignored += 1
              else:
                n_amino_acids_scored += 1
                n_rotamers = 0
                for rotamer, rotamer_sites_cart in rotamer_iterator:
                  n_rotamers += 1
                  residue.atoms().set_xyz(new_xyz=rotamer_sites_cart)
                  tmp = sites_cart_start.set_selected(
                    residue_iselection, rotamer_sites_cart)
                  refined = rsr_manager.refine_restrained(
                    tmp.select(rsel),
                    rsel, rs)
                  if(rev.is_better(sites_cart = refined.sites_cart_variable)):
                    rotamer_id_best = rotamer.id
                    residue_sites_best = refined.sites_cart_variable.deep_copy()
                residue.atoms().set_xyz(new_xyz=residue_sites_best)
                if(rotamer_id_best is not None):
                  sites_cart_start = sites_cart_start.set_selected(
                      residue_iselection, residue_sites_best)
                else:
                  tmp = sites_cart_start.set_selected(
                    residue_iselection, residue_sites_best)
                  refined = rsr_manager.refine_restrained(
                    tmp.select(rsel),
                    rsel, rs)
                  if(rev.is_better(sites_cart = refined.sites_cart_variable)):
                    rotamer_sites_cart_refined = \
                      refined.sites_cart_variable.deep_copy()
                    sites_cart_start = sites_cart_start.set_selected(
                      residue_iselection, rotamer_sites_cart_refined)
                    residue.atoms().set_xyz(new_xyz=rotamer_sites_cart_refined)
            if(abs(rev.t1_best-rev.t1_start) > 0.01 and
               abs(rev.t2_best-rev.t2_start) > 0.01):
              print >> log, fmt3 % (
                residue_id_str, cc_start, rev.t1_start, rev.t2_start,
                rev.t1_best, rev.t2_best, rotamer_id_best, n_rotamers)
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

def validate(fmodel, residue_rsr_monitor, log):
  xray_structure = fmodel.xray_structure
  map_data_1,fft_map_1 = get_map_data(
    fmodel = fmodel, map_type = "2mFo-DFc", kick=False)
  map_data_2,fft_map_2 = get_map_data(
    fmodel = fmodel, map_type = "Fc")
  map_data_3,fft_map_3 = get_map_data(
    fmodel = fmodel, map_type = "mFo-DFc", kick=False)
  map_selector = select_map(
    unit_cell  = xray_structure.unit_cell(),
    map_data_1 = map_data_1,
    map_data_2 = map_data_2)
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
    residue_name = rm.residue_id_str.strip().split()[0][1:]
    if(get_class(name=residue_name) == "common_amino_acid"):
      sites_cart_residue = sites_cart.select(rm.selection)
      t1 = target(sites_cart_residue, unit_cell, map_data_1)
      t2 = target(sites_cart_residue, unit_cell, map_data_3)
      cc = map_selector.get_cc(sites_cart = sites_cart_residue,
        residue_iselection = rm.selection)
      flag = ""
      dmif1 = rm.mfodfc < 0 and t2 < 0 and t2 < rm.mfodfc
      dmif2 = rm.mfodfc > 0 and t2 < 0 and abs(t2) > abs(rm.mfodfc)
      dmif = dmif1 or dmif2
      if((cc < rm.cc or t1 < rm.twomfodfc) and dmif): flag = " <<<"
      print >> log, fmt3 % (
        rm.residue_id_str, rm.cc, rm.twomfodfc, rm.mfodfc, cc,t1,t2), flag
      if(len(flag)>0):
        sites_cart_result = sites_cart_result.set_selected(
          rm.selection, rm.sites_cart)
  xray_structure.set_sites_cart(sites_cart_result)
  fmodel.update_xray_structure(xray_structure = xray_structure,
    update_f_calc=True, update_f_mask=True)
  print "r_work=%6.4f r_free=%6.4f" % (fmodel.r_work(), fmodel.r_free())

def run(fmodel,
        geometry_restraints_manager,
        pdb_hierarchy,
        number_of_macro_cycles,
        do_not_use_dihedrals,
        solvent_selection,
        log,
        poor_cc_threshold=0.9,
        ignore_water = False,
        filter_residual_map_value = 2.0,
        filter_2fofc_map = None):
  mon_lib_srv = mmtbx.monomer_library.server.server()
  if(do_not_use_dihedrals):
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
  if(ignore_water):
    selection = ~model.solvent_selection()
  else:
    selection = flex.bool(model.xray_structure.scatterers().size(), True)
  selection &= backbone_selections
  fmt = "Macro-cycle %2d: r_work=%6.4f r_free=%6.4f"
  print >> log, fmt%(0, fmodel.r_work(), fmodel.r_free())
  for macro_cycle in range(1,number_of_macro_cycles+1):
    map_data_1,fft_map_1 = get_map_data(fmodel = fmodel, map_type = "2mFo-DFc", kick=False)
    map_data_2,fft_map_2 = get_map_data(fmodel = fmodel, map_type = "Fc")
    map_data_3,fft_map_3 = get_map_data(fmodel = fmodel, map_type = "mFo-DFc", kick=False)
    if(filter_residual_map_value is not None):
      map_sel = flex.abs(map_data_3) < filter_residual_map_value
      map_data_3 = map_data_3.set_selected(map_sel, 0)
    if(filter_2fofc_map is not None):
      map_sel = flex.abs(map_data_1) < filter_2fofc_map
      map_data_1 = map_data_1.set_selected(map_sel, 0)
    rsr_manager = refiner(
      pdb_hierarchy               = pdb_hierarchy,
      target_map                  = map_data_1,
      geometry_restraints_manager = geometry_restraints_manager,
      real_space_target_weight    = 25,
      real_space_gradients_delta  = fmodel.f_obs.d_min()/4,
      max_iterations              = 25)
    residue_rsr_monitor = iterate_rotamers(
      pdb_hierarchy     = pdb_hierarchy,
      xray_structure    = fmodel.xray_structure,
      selection         = selection,
      map_data_1        = map_data_1,
      map_data_2        = map_data_2,
      map_data_3        = map_data_3,
      mon_lib_srv       = mon_lib_srv,
      rsr_manager       = rsr_manager,
      f_obs=fmodel.f_obs,fft_map = fft_map_1,
      poor_cc_threshold = poor_cc_threshold,
      log               = log)
    fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
    print >> log, "1:", fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    del map_data_1, map_data_2, map_data_3, fft_map_1, fft_map_2, fft_map_3
    if 1: # XXX add option
      assert model.xray_structure is fmodel.xray_structure
      params = mmtbx.refinement.real_space.master_params().extract()
      params.real_space_refinement.mode="diff_map"
      mmtbx.refinement.real_space.run(
        fmodel = fmodel, # XXX neutron ?
        model  = model,
        params = params.real_space_refinement,
        log    = log)
      assert model.xray_structure is fmodel.xray_structure
      fmodel.update_xray_structure(update_f_calc=True, update_f_mask=True)
      print >> log, "2:",fmt%(macro_cycle, fmodel.r_work(), fmodel.r_free())
    #
    if 1: # XXX add option
      validate(fmodel=fmodel, residue_rsr_monitor=residue_rsr_monitor, log=log)

from cctbx.array_family import flex
from cctbx import xray
import math,sys
from libtbx import adopt_init_args
from cctbx import adptbx
import iotbx.xplor.map
import iotbx.phil
from libtbx.str_utils import format_value
from mmtbx import find_peaks
from mmtbx.refinement import minimization
import scitbx.lbfgs
import mmtbx.utils
from mmtbx import real_space_correlation

real_space_correlation_core_params_str = real_space_correlation.core_params_str
assert real_space_correlation_core_params_str.find("atom_radius = None") >= 0
real_space_correlation_core_params_str = \
  real_space_correlation_core_params_str.replace("atom_radius = None",
  "atom_radius = 1.5")

master_params_str = """\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
    .short_caption = Minimum resolution
  mode = *auto filter_only every_macro_cycle
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
    .short_caption = Mode
  n_cycles = 1
    .type = int
    .short_caption = Number of cycles
  output_residue_name = HOH
    .type=str
    .input_size = 50
  output_chain_id = S
    .type=str
    .input_size = 50
  output_atom_name = O
    .type=str
    .input_size = 50
  scattering_type = O
    .type=str
    .help = Defines scattering factors for newly added waters
    .expert_level=2
    .input_size = 50
  primary_map_type = mFobs-DFmodel
    .type=str
  primary_map_cutoff = 3.0
    .type=float
  secondary_map_and_map_cc_filter
  {
    cc_map_1_type = "Fc"
      .type = str
    cc_map_2_type = 2mFo-DFmodel
      .type = str
    %s
  }
  h_bond_min_mac = 1.8
    .type = float
    .short_caption = H-bond minimum for solvent-model
    .expert_level = 1
  h_bond_min_sol = 1.8
    .type = float
    .short_caption = H-bond minimum for solvent-solvent
    .expert_level = 1
  h_bond_max = 3.2
    .type = float
    .short_caption = Maximum H-bond length
    .expert_level = 1
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
    .short_caption = New solvent ADP type
    .expert_level = 1
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  b_iso_min = 1.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
    .short_caption = Minimum B-factor
    .expert_level = 1
  b_iso_max = 80.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
    .short_caption = Maximum B-factor
    .expert_level = 1
  anisotropy_min = 0.1
    .type = float
    .help = For solvent refined as anisotropic: remove is less than this value
    .short_caption = Minimum anisotropic B-factor
    .expert_level = 1
  b_iso = None
    .type=float
    .help = Initial B-factor value for newly added water
    .short_caption = Initial B-factor value
    .expert_level = 1
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
    .expert_level = 1
  occupancy_min = 0.1
    .type=float
    .help = Minimum occupancy value, waters with smaller value will be rejected
    .short_caption = Minimum occupancy
  occupancy_max = 1.0
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
    .short_caption = Maximum occupancy
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
    .short_caption = Initial occupancy value
  filter_at_start = True
    .type = bool
    .expert_level = 1
  ignore_final_filtering_step = False
    .type = bool
    .expert_level=2
  correct_drifted_waters = True
    .type = bool
    .expert_level=2
  use_kick_maps = False
    .type = bool
    .expert_level=2
    .help = Use Dusan's Turk kick maps for peak picking
  kick_map
    .help = parameters for kick maps
    .expert_level=2
  {
     kick_size = 0.5
       .type = float
     number_of_kicks = 100
       .type = int
  }
"""%real_space_correlation_core_params_str

def master_params():
  return iotbx.phil.parse(master_params_str)

class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     is_neutron_scat_table,
                     params = master_params().extract(),
                     find_peaks_params = None,
                     log = None):
    adopt_init_args(self, locals())
    assert self.fmodel.xray_structure is self.model.xray_structure
    if(self.params is None): self.params = master_params().extract()
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    if(self.params.mode == "filter_only"): self.filter_only = True
    else: self.filter_only = False
    if(self.log is None): self.log = sys.stdout
    assert self.model.xray_structure == self.fmodel.xray_structure
    self.sites = None
    self.heights = None
    assert self.fmodel.xray_structure is self.model.xray_structure
    if(self.find_peaks_params.max_number_of_peaks is None):
      if(self.model.solvent_selection().count(False) > 0):
        self.find_peaks_params.max_number_of_peaks = \
          self.model.solvent_selection().count(False)
      else:
        self.find_peaks_params.max_number_of_peaks = \
          self.model.xray_structure.scatterers().size()
    assert self.fmodel.xray_structure is self.model.xray_structure
    self.move_solvent_to_the_end_of_atom_list()
    if(not self.is_water_last()):
      raise RuntimeError("Water picking failed: solvent must be last.")
    self.show(message = "Start model:")
    assert self.fmodel.xray_structure is self.model.xray_structure
    if(self.params.filter_at_start):
      self.filter_solvent()
      self.show(message = "Filtered:")
    assert self.fmodel.xray_structure is self.model.xray_structure
    if(not self.filter_only):
      assert self.params.primary_map_type is not None
      peaks = self.find_peaks(
        map_type     = self.params.primary_map_type,
        map_cutoff   = self.params.primary_map_cutoff,
        use_kick_map = self.params.use_kick_maps).peaks_mapped()
      if(peaks is not None):
        self.sites, self.heights = peaks.sites, peaks.heights
        self.add_new_solvent()
        assert self.fmodel.xray_structure is self.model.xray_structure
        self.show(message = "Just added new:")
        if(self.params.filter_at_start):
          self.filter_solvent()
          self.show(message = "Filtered:")
    assert self.fmodel.xray_structure is self.model.xray_structure
    #
    if(self.filter_only):
      self.correct_drifted_waters(map_cutoff =
          self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold)
    #
    if([self.sites, self.heights].count(None)==0):
      if(not self.filter_only and self.params.correct_drifted_waters):
        self.correct_drifted_waters(map_cutoff =
          self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold)
      for i in xrange(self.params.n_cycles):
        self.refine_adp()
        self.refine_occupancies()
    self.show(message = "Before filtering:")
    assert self.fmodel.xray_structure is self.model.xray_structure
    self.filter_solvent()
    assert self.fmodel.xray_structure is self.model.xray_structure
    ###
    if(self.params.secondary_map_and_map_cc_filter.cc_map_2_type is not None):
      self.find_peaks_2fofc()
      self.show(message = "%s map selection:"%
        self.params.secondary_map_and_map_cc_filter.cc_map_2_type)
    ###
    assert self.fmodel.xray_structure is self.model.xray_structure
    self.show(message = "Final:")
    self.move_solvent_to_the_end_of_atom_list()
    self.convert_water_adp()
    assert self.fmodel.xray_structure is self.model.xray_structure

  def convert_water_adp(self):
    sol_sel = self.model.solvent_selection().iselection()
    if(self.params.new_solvent == "isotropic"):
      self.model.xray_structure.convert_to_isotropic(selection = sol_sel)
    elif(self.params.new_solvent == "anisotropic"):
      selection_aniso = self.model.solvent_selection().deep_copy()
      selection_iso = self.model.solvent_selection().deep_copy()
      selection_aniso.set_selected(self.model.xray_structure.hd_selection(),
        False)
      selection_iso.set_selected(selection_aniso, False)
      selection_iso.set_selected(self.model.xray_structure.hd_selection(), True)
      self.model.xray_structure.convert_to_anisotropic(selection = selection_aniso)
      self.model.xray_structure.convert_to_isotropic(selection = selection_iso.iselection())
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)

  def move_solvent_to_the_end_of_atom_list(self):
    solsel = flex.bool(self.model.solvent_selection().count(False), False)
    solsel.extend(flex.bool(self.model.solvent_selection().count(True), True))
    xrs_sol =  self.model.xray_structure.select(self.model.solvent_selection())
    if(xrs_sol.hd_selection().count(True) == 0):
      self.reset_solvent(
        solvent_selection      = solsel,
        solvent_xray_structure = xrs_sol)
    self.model.renumber_water()
    self.fmodel.xray_structure = self.model.xray_structure

  def is_water_last(self):
    result = True
    sol_sel = self.model.solvent_selection()
    i_sol_sel = sol_sel.iselection()
    i_mac_sel = (~sol_sel).iselection()
    if(i_sol_sel.size() > 0 and i_mac_sel.size() > 0):
      if(flex.min_default(i_sol_sel,0)-flex.max_default(i_mac_sel,0) != 1):
        result = False
    return result

  def filter_solvent(self):
    sol_sel = self.model.solvent_selection()
    hd_sel = self.model.xray_structure.hd_selection()
    selection = self.model.xray_structure.all_selection()
    scatterers = self.model.xray_structure.scatterers()
    occ = scatterers.extract_occupancies()
    b_isos = scatterers.extract_u_iso_or_u_equiv(
      self.model.xray_structure.unit_cell()) * math.pi**2*8
    anisotropy = scatterers.anisotropy(unit_cell =
      self.model.xray_structure.unit_cell())
    #
    distance_iselection = mmtbx.utils.select_water_by_distance(
      sites_frac_all      = self.model.xray_structure.sites_frac(),
      element_symbols_all = self.model.xray_structure.scattering_types(),
      water_selection_o   = sol_sel.set_selected(hd_sel, False).iselection(),
      dist_max            = self.params.h_bond_max,
      dist_min_mac        = self.params.h_bond_min_mac,
      dist_min_sol        = self.params.h_bond_min_sol,
      unit_cell           = self.model.xray_structure.unit_cell())
    distance_selection = flex.bool(scatterers.size(), distance_iselection)
    #
    selection &= distance_selection
    selection &= b_isos >= self.params.b_iso_min
    selection &= b_isos <= self.params.b_iso_max
    selection &= occ >= self.params.occupancy_min
    selection &= occ <= self.params.occupancy_max
    # XXX selection &= anisotropy > self.params.anisotropy_min
    selection.set_selected(hd_sel, True)
    selection.set_selected(~sol_sel, True)
    # ISOR
    anisosel = anisotropy < self.params.anisotropy_min
    anisosel.set_selected(hd_sel, False)
    anisosel.set_selected(~sol_sel, False)
    self.model.xray_structure.convert_to_isotropic(selection =
      anisosel.iselection())
    self.model.xray_structure.convert_to_anisotropic(selection =
      anisosel.iselection())
    #
    xht = self.model.xh_connectivity_table()
    if(xht is not None):
      for ti in xht:
        if(not selection[ti[0]]): selection[ti[1]]=False
        if(selection[ti[0]]): selection[ti[1]]=True
    if(selection.size() != selection.count(True)):
      self.model = self.model.select(selection)
      self.fmodel.update_xray_structure(
        xray_structure = self.model.xray_structure,
        update_f_calc  = True)

  def reset_solvent(self, solvent_selection, solvent_xray_structure):
    assert solvent_selection.count(True) == \
      solvent_xray_structure.scatterers().size()
    self.model = self.model.remove_solvent()
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)

  def show(self, message):
    print >> self.log, message
    sol_sel = self.model.solvent_selection()
    xrs_mac_h = self.model.xray_structure.select(~sol_sel)
    xrs_sol_h = self.model.xray_structure.select(sol_sel)
    hd_sol = self.model.xray_structure.hd_selection().select(sol_sel)
    hd_mac = self.model.xray_structure.hd_selection().select(~sol_sel)
    xrs_sol = xrs_sol_h.select(~hd_sol)
    xrs_mac = xrs_mac_h.select(~hd_mac)
    scat = xrs_sol.scatterers()
    occ = scat.extract_occupancies()
    b_isos = scat.extract_u_iso_or_u_equiv(
      self.model.xray_structure.unit_cell()) * math.pi**2*8
    smallest_distances = xrs_mac.closest_distances(
      sites_frac      = xrs_sol.sites_frac(),
      distance_cutoff = self.find_peaks_params.map_next_to_model.\
        max_model_peak_dist).smallest_distances
    number = format_value("%-7d",scat.size())
    b_min  = format_value("%-7.2f", flex.min_default( b_isos, None))
    b_max  = format_value("%-7.2f", flex.max_default( b_isos, None))
    b_ave  = format_value("%-7.2f", flex.mean_default(b_isos, None))
    bl_min = format_value("%-7.2f", self.params.b_iso_min).strip()
    bl_max = format_value("%-7.2f", self.params.b_iso_max).strip()
    o_min  = format_value("%-7.2f", flex.min_default(occ, None))
    o_max  = format_value("%-7.2f", flex.max_default(occ, None))
    ol_min = format_value("%-7.2f", self.params.occupancy_min).strip()
    ol_max = format_value("%-7.2f", self.params.occupancy_max).strip()
    d_min  = format_value("%-7.2f", flex.min_default(smallest_distances, None))
    d_max  = format_value("%-7.2f", flex.max_default(smallest_distances, None))
    dl_min = format_value("%-7.2f",
      self.find_peaks_params.map_next_to_model.min_model_peak_dist).strip()
    dl_max = format_value("%-7.2f",
      self.find_peaks_params.map_next_to_model.max_model_peak_dist).strip()
    ani_min = format_value("%-7.2f", flex.min_default(scat.anisotropy(
      unit_cell = xrs_sol.unit_cell()), None))
    ani_min_l = format_value("%-7.2f",self.params.anisotropy_min).strip()
    print >>self.log,"  number           = %s"%number
    print >>self.log,"  b_iso_min        = %s (limit = %s)"%(b_min, bl_min)
    print >>self.log,"  b_iso_max        = %s (limit = %s)"%(b_max, bl_max)
    print >>self.log,"  b_iso_mean       = %s             "%(b_ave)
    print >>self.log,"  anisotropy_min   = %s (limit = %s)"%(ani_min,ani_min_l)
    print >>self.log,"  occupancy_min    = %s (limit = %s)"%(o_min, ol_min)
    print >>self.log,"  occupancy_max    = %s (limit = %s)"%(o_max, ol_max)
    print >>self.log,"  dist_sol_mol_min = %s (limit = %s)"%(d_min, dl_min)
    print >>self.log,"  dist_sol_mol_max = %s (limit = %s)"%(d_max, dl_max)

  def find_peaks(self, map_type, map_cutoff, use_kick_map=False):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    mnm = mmtbx.map_names(map_name_string = map_type)
    save_k_part, save_b_part = None, None
    if(abs(mnm.k) == abs(mnm.n) and self.fmodel.k_part()!=0):
      save_k_part = self.fmodel.k_part()
      save_b_part = self.fmodel.b_part()
      self.fmodel.update_core(k_part=0, b_part=0)
    result = find_peaks.manager(
      fmodel          = self.fmodel,
      map_type        = map_type,
      map_cutoff      = map_cutoff,
      params          = self.find_peaks_params,
      use_kick_map    = use_kick_map,
      kick_map_params = self.params.kick_map,
      log             = self.log)
    if(abs(mnm.k) == abs(mnm.n) and save_k_part is not None):
      self.fmodel.update_core(k_part=save_k_part, b_part=save_b_part)
    return result

  def correct_drifted_waters(self, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    find_peaks_params_drifted.resolution_factor = \
      self.find_peaks_params.resolution_factor
    peaks = find_peaks.manager(fmodel     = self.fmodel,
                               map_type   = "2mFobs-DFmodel",
                               map_cutoff = map_cutoff,
                               params     = find_peaks_params_drifted,
                               log        = self.log).peaks_mapped()
    if(peaks is not None and self.fmodel.r_work() > 0.01):
      sites_frac, heights = peaks.sites, peaks.heights
      model_sites_frac = self.model.xray_structure.sites_frac()
      solvent_selection = self.model.solvent_selection()
      mmtbx.utils.correct_drifted_waters(
        sites_frac_all   = model_sites_frac,
        sites_frac_peaks = sites_frac,
        water_selection  = solvent_selection,
        unit_cell        = self.model.xray_structure.unit_cell())
      self.model.xray_structure.set_sites_frac(sites_frac = model_sites_frac)
      self.fmodel.update_xray_structure(
        xray_structure = self.model.xray_structure,
        update_f_calc  = True)

  def find_peaks_2fofc(self):
    if(self.fmodel.twin): # XXX Make it possible when someone consolidates fmodels.
      print >> self.log, "Map CC and map value based filtering is disabled for twin refinement."
      return
    print >> self.log, "Before RSCC filtering: ", \
      self.model.solvent_selection().count(True)
    assert self.fmodel.xray_structure is self.model.xray_structure
    assert len(list(self.model.pdb_hierarchy.atoms_with_labels())) == \
      self.model.xray_structure.scatterers().size()
    from mmtbx import real_space_correlation
    par = self.params.secondary_map_and_map_cc_filter
    selection = self.model.solvent_selection()
    rscc_and_map_result = real_space_correlation.simple(
      fmodel                = self.fmodel,
      pdb_hierarchy         = self.model.pdb_hierarchy,
      map_1_name            = par.cc_map_1_type,
      map_2_name            = par.cc_map_2_type,
      diff_map_name         = None,
      number_of_grid_points = par.number_of_grid_points,
      atom_radius           = par.atom_radius,
      details_level         = "atom",
      selection             = selection,
      show                  = False,
      set_cc_to_zero_if_n_grid_points_less_than = par.set_cc_to_zero_if_n_grid_points_less_than,
      poor_cc_threshold                         = par.poor_cc_threshold,
      poor_map_1_value_threshold                = par.poor_map_value_threshold,
      poor_map_2_value_threshold                = par.poor_map_value_threshold)
    scatterers = self.model.xray_structure.scatterers()
    for rcc_res in rscc_and_map_result:
      try:
        i_seqs = [rcc_res.i_seq]
      except:
        i_seqs = rcc_res.residue.selection
      for i_seq in i_seqs:
        assert selection[i_seq]
        if(rcc_res.poor_flag and not
           rcc_res.scatterer.element_symbol().strip().upper() in ["H","D"]):
          selection[i_seq] = False
    sol_sel = self.model.solvent_selection()
    hd_sel = self.model.xray_structure.hd_selection()
    selection.set_selected(hd_sel, True)
    selection.set_selected(~sol_sel, True)
    xht = self.model.xh_connectivity_table()
    if(xht is not None):
      for ti in xht:
        if(not selection[ti[0]]): selection[ti[1]]=False
        if(selection[ti[0]]): selection[ti[1]]=True
    if(selection.size() != selection.count(True)):
      self.model = self.model.select(selection)
      self.fmodel.update_xray_structure(
        xray_structure = self.model.xray_structure,
        update_f_calc  = True)
    print >> self.log, "After RSCC filtering: ", \
      self.model.solvent_selection().count(True)

  def add_new_solvent(self):
    if(self.params.b_iso is None):
      sol_sel = self.model.solvent_selection()
      xrs_mac_h = self.model.xray_structure.select(~sol_sel)
      hd_mac = self.model.xray_structure.hd_selection().select(~sol_sel)
      xrs_mac = xrs_mac_h.select(~hd_mac)
      b = xrs_mac.extract_u_iso_or_u_equiv() * math.pi**2*8
      b_solv = flex.mean_default(b, None)
      if(b_solv is not None and b_solv < self.params.b_iso_min or
         b_solv > self.params.b_iso_max):
        b_solv = (self.params.b_iso_min + self.params.b_iso_max) / 2.
    else:
      b_solv = self.params.b_iso
    if(self.params.new_solvent == "isotropic"):
      new_scatterers = flex.xray_scatterer(
        self.sites.size(),
        xray.scatterer(occupancy       = self.params.occupancy,
                       b               = b_solv,
                       scattering_type = self.params.scattering_type))
    elif(self.params.new_solvent == "anisotropic"):
      u_star = adptbx.u_iso_as_u_star(self.model.xray_structure.unit_cell(),
        adptbx.b_as_u(b_solv))
      new_scatterers = flex.xray_scatterer(
        self.sites.size(),
        xray.scatterer(
          occupancy       = self.params.occupancy,
          u               = u_star,
          scattering_type = self.params.scattering_type))
    else: raise RuntimeError
    new_scatterers.set_sites(self.sites)
    solvent_xray_structure = xray.structure(
      special_position_settings = self.model.xray_structure,
      scatterers                = new_scatterers)
    xrs_sol = self.model.xray_structure.select(self.model.solvent_selection())
    xrs_mac = self.model.xray_structure.select(~self.model.solvent_selection())
    xrs_sol = xrs_sol.concatenate(other = solvent_xray_structure)
    sol_sel = flex.bool(xrs_mac.scatterers().size(), False)
    sol_sel.extend( flex.bool(xrs_sol.scatterers().size(), True) )
    self.model.add_solvent(
      solvent_xray_structure = solvent_xray_structure,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id,
      refine_occupancies     = self.params.refine_occupancies,
      refine_adp             = self.params.new_solvent)
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)

  def refine_adp(self):
    if(not self.filter_only and self.params.refine_adp and
       self.model.refinement_flags.individual_adp and
       self.model.solvent_selection().count(True) > 0):
      self.fmodels.update_xray_structure(
         xray_structure = self.model.xray_structure,
         update_f_calc  = True,
         update_f_mask  = True)
      print >> self.log, \
        "ADP refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free())
      if(self.params.new_solvent == "anisotropic"):
        selection_aniso = self.model.solvent_selection().deep_copy()
        selection_iso = self.model.solvent_selection().deep_copy()
        selection_aniso.set_selected(self.model.xray_structure.hd_selection(),
          False)
        selection_iso.set_selected(selection_aniso, False)
        selection_iso.set_selected(self.model.xray_structure.hd_selection(),
          False)
        #
        self.model.set_refine_individual_adp(selection_aniso = selection_aniso,
                                             selection_iso   = selection_iso)
      else:
        selection_iso = self.model.solvent_selection().deep_copy()
        selection_iso.set_selected(self.model.xray_structure.hd_selection(),
          False)
        self.model.set_refine_individual_adp(selection_iso = selection_iso)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = 25)
      minimized = minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = self.fmodels,
        model                    = self.model,
        is_neutron_scat_table    = self.is_neutron_scat_table,
        refine_adp               = True,
        lbfgs_termination_params = lbfgs_termination_params)
      print >> self.log,\
        "ADP refinement (water only), final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free())

  def refine_occupancies(self):
    if(not self.filter_only and self.params.refine_occupancies and
       self.model.refinement_flags.occupancies and
       self.model.solvent_selection().count(True) > 0):
      self.fmodels.update_xray_structure(
         xray_structure = self.model.xray_structure,
         update_f_calc  = True,
         update_f_mask  = True)
      print >> self.log,\
        "occupancy refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free())
      self.fmodels.fmodel_xray().xray_structure.scatterers().flags_set_grads(
        state = False)
      i_selection = self.model.solvent_selection().iselection()
      self.fmodels.fmodel_xray().xray_structure.scatterers(
        ).flags_set_grad_occupancy(iselection = i_selection)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations = 25)
      minimized = mmtbx.refinement.minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = self.fmodels,
        model                    = self.model,
        is_neutron_scat_table    = self.is_neutron_scat_table,
        lbfgs_termination_params = lbfgs_termination_params)
      self.fmodels.fmodel_xray().xray_structure.adjust_occupancy(
        occ_max   = self.params.occupancy_max,
        occ_min   = self.params.occupancy_min,
        selection = i_selection)
      #
      print >> self.log,\
        "occupancy refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free())

def show_histogram(data,
                   n_slots,
                   out=None,
                   prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix
    histogram = flex.histogram(data    = data,
                               n_slots = n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print >> out, "%7.3f - %7.3f: %d" % (low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff
    out.flush()
    return histogram

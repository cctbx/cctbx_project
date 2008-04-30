from iotbx import pdb
from cctbx.array_family import flex
from cctbx import xray
import math,sys
from mmtbx import max_lik
from mmtbx.max_lik import maxlik
from mmtbx import masks
from mmtbx import bulk_solvent
import mmtbx.restraints
from cctbx import crystal
from mmtbx.max_lik import max_like_non_uniform
from libtbx import adopt_init_args
import mmtbx.f_model
from libtbx import introspection
from cctbx import maptbx
from cctbx import adptbx
from cctbx import crystal
import libtbx.load_env
import iotbx.xplor.map
from scitbx import matrix
from libtbx.test_utils import approx_equal
import iotbx.phil
from libtbx.str_utils import format_value
from mmtbx import find_peaks
from mmtbx.refinement import minimization
import scitbx.lbfgs
import mmtbx.utils

master_params = iotbx.phil.parse("""\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
  mode = *auto filter_only every_macro_cycle
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
  output_residue_name = HOH
    .type=str
  output_chain_id = S
    .type=str
  output_atom_name = O
    .type=str
  b_iso_min = 1.0
    .type=float
    .help = Minimum B-factor value, waters with smaller value will be rejected
  b_iso_max = 80.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
  anisotropy_min = 0.1
    .type = float
    .help = For solvent refined as anisotropic: remove is less than this value
  b_iso = None
    .type=float
    .help = Initial B-factor value for newly added water
  scattering_type = O
    .type=str
    .help = Defines scattering factors for newly added waters
    .expert_level=2
  occupancy_min = 0.1
    .type=float
    .help = Minimum occupancy value, waters with smaller value will be rejected
  occupancy_max = 1.0
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
  primary_map_type = mFobs-DFmodel
    .type=str
  primary_map_cutoff = 3.0
    .type=float
  secondary_map_type = 2mFobs-DFmodel
    .type=str
  secondary_map_cutoff = 1.0
    .type=float
  h_bond_min_mac = 1.8
    .type = float
  h_bond_min_sol = 1.8
    .type = float
  h_bond_max = 3.2
    .type = float
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
  filter_at_start = True
    .type = bool
  n_cycles = 1
    .type = int
  ignore_final_filtering_step = False
    .type = bool
    .expert_level=2
  correct_drifted_waters = True
    .type = bool
    .expert_level=2
""")

class water_ids(object):
  residue_names = ["HOH","SOL","SOLV","WAT","DOD","TIP3"]
  atom_names = ["O","OH2","H","H1","H2","D","D1","D2"]
  element_types = ["O","H","D","", " "]


class manager(object):
  def __init__(self, fmodel,
                     fmodels,
                     model,
                     params = master_params.extract(),
                     find_peaks_params = None,
                     log    = None):
    adopt_init_args(self, locals())
    if(self.params is None): self.params = master_params.extract()
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    if(self.params.mode == "filter_only"): self.filter_only = True
    else: self.filter_only = False
    if(self.log is None): self.log = sys.stdout
    assert self.model.xray_structure == self.fmodel.xray_structure
    self.sites = None
    self.heights = None
    if(self.find_peaks_params.max_number_of_peaks is None):
      if(self.model.solvent_selection().count(False) > 0):
        self.find_peaks_params.max_number_of_peaks = \
          self.model.solvent_selection().count(False)
      else:
        self.find_peaks_params.max_number_of_peaks = \
          self.model.xray_structure.scatterers().size()
    self.move_solvent_to_the_end_of_atom_list()
    if(not self.is_water_last()):
      raise RuntimeError("Water picking failed: solvent must be last.")
    self.show(message = "Start model:")
    if(self.params.filter_at_start):
      self.filter_solvent()
      self.show(message = "Filtered:")
    if(not self.filter_only):
      assert self.params.primary_map_type is not None
      peaks = self.find_peaks(
        map_type   = self.params.primary_map_type,
        map_cutoff = self.params.primary_map_cutoff).peaks_mapped()
      self.sites, self.heights = peaks.sites, peaks.heights
      self.add_new_solvent()
      self.show(message = "Just added new:")
      if(self.params.filter_at_start):
        self.filter_solvent()
        self.show(message = "Filtered:")
    #
    if(not self.filter_only and self.params.correct_drifted_waters):
      self.correct_drifted_waters(map_cutoff = self.params.secondary_map_cutoff)
    #
    for i in xrange(self.params.n_cycles):
      self.refine_adp()
      self.refine_occupancies()
    #
    if(not self.filter_only):
      if(self.params.secondary_map_type is not None):
        self.find_peaks_2fofc()
        self.show(message = "2Fo-Fc map selection:")
    self.show(message = "Before filtering:")
    self.filter_solvent()
    self.show(message = "Final:")
    self.move_solvent_to_the_end_of_atom_list()
    self.convert_water_adp()

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
      #
      self.model.xray_structure.convert_to_anisotropic(selection = selection_aniso)
      self.model.xray_structure.convert_to_isotropic(selection = selection_iso.iselection())

  def move_solvent_to_the_end_of_atom_list(self):
    solsel = flex.bool(self.model.solvent_selection().count(False), False)
    solsel.extend(flex.bool(self.model.solvent_selection().count(True), True))
    xrs_sol =  self.model.xray_structure.select(self.model.solvent_selection())
    if(xrs_sol.hd_selection().count(True) == 0):
      self.reset_solvent(
        solvent_selection      = solsel,
        solvent_xray_structure = xrs_sol)
    self.model.renumber_water()

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
    self.model = self.model.select(selection)

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

  def find_peaks(self, map_type, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    return find_peaks.manager(fmodel     = self.fmodel,
                              map_type   = map_type,
                              map_cutoff = map_cutoff,
                              params     = self.find_peaks_params,
                              log        = self.log)

  def correct_drifted_waters(self, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.xray_structure,
      update_f_calc  = True)
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    peaks = find_peaks.manager(fmodel     = self.fmodel,
                               map_type   = "2mFobs-DFmodel",
                               map_cutoff = map_cutoff,
                               params     = find_peaks_params_drifted,
                               log        = self.log).peaks_mapped()
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
    fft_map = self.fmodel.electron_density_map(
      map_type          = self.params.secondary_map_type,
      resolution_factor = self.find_peaks_params.resolution_factor,
      symmetry_flags    = maptbx.use_space_group_symmetry)
    fft_map.apply_sigma_scaling()
    fft_map_data = fft_map.real_map_unpadded()
    selection = self.model.solvent_selection()
    scatterers = self.model.xray_structure.scatterers()
    for i_seq, sel_i in enumerate(selection):
      if(sel_i):
        ed_val = maptbx.eight_point_interpolation(fft_map_data,
          scatterers[i_seq].site)
        if(ed_val < self.params.secondary_map_cutoff):
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
    self.model = self.model.select(selection)

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

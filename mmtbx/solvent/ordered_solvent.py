from __future__ import absolute_import, division, print_function
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
from cctbx import miller
from cctbx import maptbx
from libtbx.test_utils import approx_equal
from six.moves import range

output_params_str = """
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
"""

h_bond_params_str = """
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
"""

adp_occ_params_str = """
  new_solvent = *isotropic anisotropic
    .type = choice
    .help = Based on the choice, added solvent will have isotropic or \
            anisotropic b-factors
    .short_caption = New solvent ADP type
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
"""

master_params_str = """\
  low_resolution = 2.8
    .type = float
    .help = Low resolution limit for water picking (at lower resolution water \
            will not be picked even if requessted)
    .short_caption = Minimum resolution
  mode = *second_half filter_only every_macro_cycle every_macro_cycle_after_first
    .type=choice
    .help = Choices for water picking strategy: auto - start water picking \
            after ferst few macro-cycles, filter_only - remove water only, \
            every_macro_cycle - do water update every macro-cycle
    .short_caption = Mode
  n_cycles = 1
    .type = int
    .short_caption = Number of cycles
  %s
  primary_map_type = mFobs-DFmodel
    .type=str
    .help = Map used to identify candidate water sites - by default this is \
      the standard difference map.
  primary_map_cutoff = 3.0
    .type=float
    .short_caption = Primary map cutoff (sigma)
  secondary_map_and_map_cc_filter
    .short_caption = Secondary map filter
    .style = auto_align box
  {
    cc_map_1_type = "Fc"
      .type = str
      .short_caption = Model map type for CC calculation
    cc_map_2_type = 2mFo-DFmodel
      .type = str
      .short_caption = Experimental map type for CC calculation
    poor_cc_threshold = 0.7
      .type = float
      .short_caption = Minimum map correlation
    poor_map_value_threshold = 1.0
      .type = float
      .short_caption = Minimum map value (sigma)
  }
  %s
  refine_adp = True
    .type = bool
    .help = Refine ADP for newly placed solvent.
    .short_caption = Refine new solvent ADPs
    .expert_level = 1
  refine_occupancies = False
    .type = bool
    .help = Refine solvent occupancies.
    .expert_level = 1
  %s
  filter_at_start = True
    .type = bool
    .expert_level = 1
  ignore_final_filtering_step = False
    .type = bool
    .expert_level=2
  correct_drifted_waters = True
    .type = bool
    .expert_level=2
""" % (output_params_str, h_bond_params_str, adp_occ_params_str)

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
    assert (0 < len(self.params.output_atom_name) <= 4)
    assert (len(self.params.output_chain_id) <= 2)
    assert (0 < len(self.params.output_residue_name) <= 3)
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    if(self.params is None): self.params = master_params().extract()
    if(self.find_peaks_params is None):
      self.find_peaks_params = find_peaks.master_params.extract()
    if(self.params.mode == "filter_only"): self.filter_only = True
    else: self.filter_only = False
    if(self.log is None): self.log = sys.stdout
    assert self.model.get_xray_structure() == self.fmodel.xray_structure
    self.sites = None
    self.heights = None
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    if(self.find_peaks_params.max_number_of_peaks is None):
      if(self.model.solvent_selection().count(False) > 0):
        self.find_peaks_params.max_number_of_peaks = \
          self.model.solvent_selection().count(False)
      else:
        self.find_peaks_params.max_number_of_peaks = \
          self.model.get_number_of_atoms()
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    self.move_solvent_to_the_end_of_atom_list()
    if(not self.is_water_last()):
      raise RuntimeError("Water picking failed: solvent must be last.")
    self.show(message = "Start model:")
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    if(self.params.filter_at_start):
      self.filter_solvent()
      self.show(message = "Filtered:")
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    if(not self.filter_only):
      assert self.params.primary_map_type is not None
      peaks = self.find_peaks(
        map_type   = self.params.primary_map_type,
        map_cutoff = self.params.primary_map_cutoff).peaks_mapped()
      if(peaks is not None):
        self.sites, self.heights = peaks.sites, peaks.heights
        self.add_new_solvent()
        assert self.fmodel.xray_structure is self.model.get_xray_structure()
        self.show(message = "Just added new:")
        if(self.params.filter_at_start):
          self.filter_solvent()
          self.show(message = "Filtered:")
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    #
    if(self.filter_only):
      self.correct_drifted_waters(map_cutoff =
          self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold)
    #
    if([self.sites, self.heights].count(None)==0):
      if(not self.filter_only and self.params.correct_drifted_waters):
        self.correct_drifted_waters(map_cutoff =
          self.params.secondary_map_and_map_cc_filter.poor_map_value_threshold)
      for i in range(self.params.n_cycles):
        self.refine_adp()
        self.refine_occupancies()
    self.show(message = "Before filtering:")
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    self.filter_solvent()
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    ###
    if(self.params.secondary_map_and_map_cc_filter.cc_map_2_type is not None):
      self.find_peaks_2fofc()
      self.show(message = "%s map selection:"%
        self.params.secondary_map_and_map_cc_filter.cc_map_2_type)
    ###
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    self.show(message = "Final:")
    self.move_solvent_to_the_end_of_atom_list()
    self.convert_water_adp()
    assert self.fmodel.xray_structure is self.model.get_xray_structure()

  def convert_water_adp(self):
    hd_sel     = self.model.get_hd_selection()
    not_hd_sel = ~hd_sel
    sol_sel    = self.model.solvent_selection()
    not_sol_sel= ~sol_sel
    selection_aniso = self.model.get_xray_structure().use_u_aniso().deep_copy()
    if(self.params.new_solvent == "anisotropic"):
      selection_aniso.set_selected(sol_sel, True)
    selection_iso   = self.model.get_xray_structure().use_u_iso().deep_copy()
    selection_aniso.set_selected(not_sol_sel, False)
    selection_iso  .set_selected(not_sol_sel, False)
    if(not self.is_neutron_scat_table):
      selection_aniso.set_selected(hd_sel, False)
      selection_iso.set_selected(hd_sel, False)
    selection_aniso.set_selected(selection_iso, False)
    selection_iso.set_selected(selection_aniso, False)
    self.model.get_xray_structure().convert_to_anisotropic(selection = selection_aniso)
    self.model.get_xray_structure().convert_to_isotropic(selection = selection_iso.iselection())
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True)

  def move_solvent_to_the_end_of_atom_list(self):
    solsel = flex.bool(self.model.solvent_selection().count(False), False)
    solsel.extend(flex.bool(self.model.solvent_selection().count(True), True))
    xrs_sol =  self.model.get_xray_structure().select(self.model.solvent_selection())
    if(xrs_sol.hd_selection().count(True) == 0):
      self.reset_solvent(
        solvent_selection      = solsel,
        solvent_xray_structure = xrs_sol)
    self.model.renumber_water()
    self.fmodel.xray_structure = self.model.get_xray_structure()

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
    hd_sel = self.model.get_hd_selection()
    selection = self.model.get_xray_structure().all_selection()
    scatterers = self.model.get_xray_structure().scatterers()
    occ = scatterers.extract_occupancies()
    b_isos = scatterers.extract_u_iso_or_u_equiv(
      self.model.get_xray_structure().unit_cell()) * math.pi**2*8
    anisotropy = scatterers.anisotropy(unit_cell =
      self.model.get_xray_structure().unit_cell())
    #
    distance_iselection = mmtbx.utils.select_water_by_distance(
      sites_frac_all      = self.model.get_xray_structure().sites_frac(),
      element_symbols_all = self.model.get_xray_structure().scattering_types(),
      water_selection_o   = sol_sel.set_selected(hd_sel, False).iselection(),
      dist_max            = self.params.h_bond_max,
      dist_min_mac        = self.params.h_bond_min_mac,
      dist_min_sol        = self.params.h_bond_min_sol,
      unit_cell           = self.model.get_xray_structure().unit_cell())
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
    self.model.get_xray_structure().convert_to_isotropic(selection =
      anisosel.iselection())
    self.model.get_xray_structure().convert_to_anisotropic(selection =
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
        xray_structure = self.model.get_xray_structure(),
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
    print(message, file=self.log)
    sol_sel = self.model.solvent_selection()
    xrs_mac_h = self.model.get_xray_structure().select(~sol_sel)
    xrs_sol_h = self.model.get_xray_structure().select(sol_sel)
    hd_sol = self.model.get_hd_selection().select(sol_sel)
    hd_mac = self.model.get_hd_selection().select(~sol_sel)
    xrs_sol = xrs_sol_h.select(~hd_sol)
    xrs_mac = xrs_mac_h.select(~hd_mac)
    scat = xrs_sol.scatterers()
    occ = scat.extract_occupancies()
    b_isos = scat.extract_u_iso_or_u_equiv(
      self.model.get_xray_structure().unit_cell()) * math.pi**2*8
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
    print("  number           = %s"%number, file=self.log)
    print("  b_iso_min        = %s (limit = %s)"%(b_min, bl_min), file=self.log)
    print("  b_iso_max        = %s (limit = %s)"%(b_max, bl_max), file=self.log)
    print("  b_iso_mean       = %s             "%(b_ave), file=self.log)
    print("  anisotropy_min   = %s (limit = %s)"%(ani_min,ani_min_l), file=self.log)
    print("  occupancy_min    = %s (limit = %s)"%(o_min, ol_min), file=self.log)
    print("  occupancy_max    = %s (limit = %s)"%(o_max, ol_max), file=self.log)
    print("  dist_sol_mol_min = %s (limit = %s)"%(d_min, dl_min), file=self.log)
    print("  dist_sol_mol_max = %s (limit = %s)"%(d_max, dl_max), file=self.log)

  def find_peaks(self, map_type, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True)
    return find_peaks.manager(
      fmodel     = self.fmodel,
      map_type   = map_type,
      map_cutoff = map_cutoff,
      params     = self.find_peaks_params,
      log        = self.log)

  def correct_drifted_waters(self, map_cutoff):
    self.fmodel.update_xray_structure(
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True)
    find_peaks_params_drifted = find_peaks.master_params.extract()
    find_peaks_params_drifted.map_next_to_model.min_model_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.min_peak_peak_dist=0.01
    find_peaks_params_drifted.map_next_to_model.max_model_peak_dist=0.5
    find_peaks_params_drifted.peak_search.min_cross_distance=0.5
    find_peaks_params_drifted.resolution_factor = \
      self.find_peaks_params.resolution_factor
    peaks = find_peaks.manager(
      fmodel         = self.fmodel,
      map_type       = "2mFobs-DFmodel",
      map_cutoff     = map_cutoff,
      params         = find_peaks_params_drifted,
      log            = self.log).peaks_mapped()
    if(peaks is not None and self.fmodel.r_work() > 0.01):
      sites_frac, heights = peaks.sites, peaks.heights
      model_sites_frac = self.model.get_xray_structure().sites_frac()
      solvent_selection = self.model.solvent_selection()
      mmtbx.utils.correct_drifted_waters(
        sites_frac_all   = model_sites_frac,
        sites_frac_peaks = sites_frac,
        water_selection  = solvent_selection,
        unit_cell        = self.model.get_xray_structure().unit_cell())
      self.model.get_xray_structure().set_sites_frac(sites_frac = model_sites_frac)
      self.fmodel.update_xray_structure(
        xray_structure = self.model.get_xray_structure(),
        update_f_calc  = True)

  def find_peaks_2fofc(self):
    if(self.fmodel.twin): # XXX Make it possible when someone consolidates fmodels.
      print("Map CC and map value based filtering is disabled for twin refinement.", file=self.log)
      return
    print("Before RSCC filtering: ", \
      self.model.solvent_selection().count(True), file=self.log)
    assert self.fmodel.xray_structure is self.model.get_xray_structure()
    assert len(list(self.model.get_hierarchy().atoms_with_labels())) == \
      self.model.get_number_of_atoms()
    par = self.params.secondary_map_and_map_cc_filter
    selection = self.model.solvent_selection()
    # filter by map cc and value
    e_map_obj = self.fmodel.electron_density_map()
    coeffs_1 = e_map_obj.map_coefficients(
      map_type     = par.cc_map_1_type,
      fill_missing = False,
      isotropize   = True)
    coeffs_2 = e_map_obj.map_coefficients(
      map_type     = par.cc_map_2_type,
      fill_missing = False,
      isotropize   = True)
    fft_map_1 = coeffs_1.fft_map(resolution_factor = 1./4)
    fft_map_1.apply_sigma_scaling()
    map_1 = fft_map_1.real_map_unpadded()
    fft_map_2 = miller.fft_map(
      crystal_gridding     = fft_map_1,
      fourier_coefficients = coeffs_2)
    fft_map_2.apply_sigma_scaling()
    map_2 = fft_map_2.real_map_unpadded()
    sites_cart = self.fmodel.xray_structure.sites_cart()
    sites_frac = self.fmodel.xray_structure.sites_frac()
    scatterers = self.model.get_xray_structure().scatterers()
    assert approx_equal(self.model.get_xray_structure().sites_frac(), sites_frac)
    unit_cell = self.fmodel.xray_structure.unit_cell()
    for i, sel_i in enumerate(selection):
      if(sel_i):
        sel = maptbx.grid_indices_around_sites(
          unit_cell  = unit_cell,
          fft_n_real = map_1.focus(),
          fft_m_real = map_1.all(),
          sites_cart = flex.vec3_double([sites_cart[i]]),
          site_radii = flex.double([1.5]))
        cc = flex.linear_correlation(x=map_1.select(sel),
          y=map_2.select(sel)).coefficient()
        map_value_1 = map_1.eight_point_interpolation(sites_frac[i])
        map_value_2 = map_2.eight_point_interpolation(sites_frac[i])
        if((cc < par.poor_cc_threshold or
           map_value_1 < par.poor_map_value_threshold or
           map_value_2 < par.poor_map_value_threshold) and not
           scatterers[i].element_symbol().strip().upper() in ["H","D"]):
          selection[i]=False
    #
    sol_sel = self.model.solvent_selection()
    hd_sel = self.model.get_hd_selection()
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
        xray_structure = self.model.get_xray_structure(),
        update_f_calc  = True)
    print("After RSCC filtering: ", \
      self.model.solvent_selection().count(True), file=self.log)

  def add_new_solvent(self):
    if(self.params.b_iso is None):
      sol_sel = self.model.solvent_selection()
      xrs_mac_h = self.model.get_xray_structure().select(~sol_sel)
      hd_mac = self.model.get_hd_selection().select(~sol_sel)
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
      u_star = adptbx.u_iso_as_u_star(self.model.get_xray_structure().unit_cell(),
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
      special_position_settings = self.model.get_xray_structure(),
      scatterers                = new_scatterers)
    xrs_sol = self.model.get_xray_structure().select(self.model.solvent_selection())
    xrs_mac = self.model.get_xray_structure().select(~self.model.solvent_selection())
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
      xray_structure = self.model.get_xray_structure(),
      update_f_calc  = True)

  def refine_adp(self):
    if(not self.filter_only and self.params.refine_adp and
       self.model.refinement_flags.individual_adp and
       self.model.solvent_selection().count(True) > 0):
      self.fmodels.update_xray_structure(
         xray_structure = self.model.get_xray_structure(),
         update_f_calc  = True,
         update_f_mask  = True)
      print("ADP refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)
      # set refinement flags (not exercised!)
      hd_sel     = self.model.get_hd_selection()
      not_hd_sel = ~hd_sel
      sol_sel    = self.model.solvent_selection()
      not_sol_sel= ~sol_sel
      selection_aniso = self.model.get_xray_structure().use_u_aniso().deep_copy()
      if(self.params.new_solvent == "anisotropic"):
        selection_aniso.set_selected(sol_sel, True)
      selection_iso   = self.model.get_xray_structure().use_u_iso().deep_copy()
      selection_aniso.set_selected(not_sol_sel, False)
      selection_iso  .set_selected(not_sol_sel, False)
      if(not self.is_neutron_scat_table):
        selection_aniso.set_selected(hd_sel, False)
        selection_iso.set_selected(hd_sel, False)
      selection_aniso.set_selected(selection_iso, False)
      selection_iso.set_selected(selection_aniso, False)
      self.model.set_refine_individual_adp(
        selection_aniso = selection_aniso, selection_iso = selection_iso)
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
          max_iterations = 25)
      minimized = minimization.lbfgs(
        restraints_manager       = None,
        fmodels                  = self.fmodels,
        model                    = self.model,
        is_neutron_scat_table    = self.is_neutron_scat_table,
        refine_adp               = True,
        lbfgs_termination_params = lbfgs_termination_params)
      print("ADP refinement (water only), final r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

  def refine_occupancies(self):
    if(not self.filter_only and self.params.refine_occupancies and
       self.model.refinement_flags.occupancies and
       self.model.solvent_selection().count(True) > 0):
      self.fmodels.update_xray_structure(
         xray_structure = self.model.get_xray_structure(),
         update_f_calc  = True,
         update_f_mask  = True)
      print("occupancy refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)
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
      print("occupancy refinement (water only), start r_work=%6.4f r_free=%6.4f"%(
        self.fmodel.r_work(), self.fmodel.r_free()), file=self.log)

def show_histogram(data,
                   n_slots,
                   out=None,
                   prefix=""):
    if (out is None): out = sys.stdout
    print(prefix, file=out)
    histogram = flex.histogram(data    = data,
                               n_slots = n_slots)
    low_cutoff = histogram.data_min()
    for i,n in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print("%7.3f - %7.3f: %d" % (low_cutoff, high_cutoff, n), file=out)
      low_cutoff = high_cutoff
    out.flush()
    return histogram

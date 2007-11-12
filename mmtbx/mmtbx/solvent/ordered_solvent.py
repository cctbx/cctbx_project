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
import libtbx.load_env
import iotbx.xplor.map
from scitbx import matrix
from libtbx.test_utils import approx_equal
import iotbx.phil
from libtbx.str_utils import format_value

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
  b_iso_max = 50.0
    .type=float
    .help = Maximum B-factor value, waters with bigger value will be rejected
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
  occupancy_max = 1.2
    .type=float
    .help = Maximum occupancy value, waters with bigger value will be rejected
  occupancy = 1.0
    .type=float
    .help = Initial occupancy value for newly added water
  bulk_solvent_mask_exclusions = True
    .type = bool
    .help = Do water selection based on bulk-solvent mask
  use_sigma_scaled_maps = True
    .type=bool
    .help = Use sigma scales maps for water pick picking
  primary_map_type = m*Fobs-D*Fmodel
    .type=str
  primary_map_k = None
    .type=float
  primary_map_n = None
    .type=float
  primary_map_cutoff = 3.0
    .type=float
  secondary_map_type = 2m*Fobs-D*Fmodel
    .type=str
  secondary_map_k = None
    .type=float
  secondary_map_n = None
    .type=float
  secondary_map_cutoff = 1.0
    .type=float
  peak_map_matching_tolerance = 2.0
    .type=float
  resolution_factor = 1./4.
    .type=float
  min_solv_macromol_dist = 1.8
    .type=float
  max_solv_macromol_dist = 6.0
    .type=float
  min_solv_solv_dist = 1.8
    .type=float
  max_number_of_peaks = None
    .type=int
  verbose = 1
    .type=int
  peak_search
    .expert_level=1
  {
    peak_search_level = 1
      .type=int
    max_peaks = 0
      .type=int
    interpolate = True
      .type=bool
    min_distance_sym_equiv = 1.e-6
      .type=float
    general_positions_only = False
      .type=bool
    min_cross_distance = 2.0
      .type=float
  }
""")


class manager(object):
  def __init__(self, fmodel,
                     model,
                     solvent_selection,
                     params = master_params.extract(),
                     log    = None):
    adopt_init_args(self, locals())
    #
    self.b_iso_min                   = self.params.b_iso_min
    self.b_iso_max                   = self.params.b_iso_max
    self.b_iso                       = self.params.b_iso
    self.scattering_type             = self.params.scattering_type
    self.occupancy_min               = self.params.occupancy_min
    self.occupancy_max               = self.params.occupancy_max
    self.occupancy                   = self.params.occupancy
    self.use_sigma_scaled_maps       = self.params.use_sigma_scaled_maps
    self.bulk_solvent_mask_exclusions= self.params.bulk_solvent_mask_exclusions
    self.primary_map_type            = self.params.primary_map_type
    self.primary_map_k               = self.params.primary_map_k
    self.primary_map_n               = self.params.primary_map_n
    self.primary_map_cutoff          = self.params.primary_map_cutoff
    self.secondary_map_type          = self.params.secondary_map_type
    self.secondary_map_k             = self.params.secondary_map_k
    self.secondary_map_n             = self.params.secondary_map_n
    self.secondary_map_cutoff        = self.params.secondary_map_cutoff
    self.peak_map_matching_tolerance = self.params.peak_map_matching_tolerance
    self.resolution_factor           = self.params.resolution_factor
    self.min_solv_macromol_dist      = self.params.min_solv_macromol_dist
    self.max_solv_macromol_dist      = self.params.max_solv_macromol_dist
    self.min_solv_solv_dist          = self.params.min_solv_solv_dist
    self.max_number_of_peaks         = self.params.max_number_of_peaks
    if(self.params.mode == "filter_only"): self.filter_only = True
    else: self.filter_only = False
    self.verbose                     = self.params.verbose
    self.peak_search_level     = self.params.peak_search.peak_search_level
    self.max_peaks             = self.params.peak_search.max_peaks
    self.interpolate           = self.params.peak_search.interpolate
    self.min_distance_sym_equiv= self.params.peak_search.min_distance_sym_equiv
    self.general_positions_only= self.params.peak_search.general_positions_only
    self.min_cross_distance    = self.params.peak_search.min_cross_distance
    #
    if (self.log is None): self.log = sys.stdout
    self.xray_structure = self.fmodel.xray_structure.deep_copy_scatterers()
    self.sites = None
    self.heights = None
    if(self.max_number_of_peaks is None and
                                      self.solvent_selection.count(False) > 0):
       self.max_number_of_peaks= int(self.solvent_selection.count(False)/10*10)
    else:
       self.max_number_of_peaks=1000
    if(self.verbose > 0): self.show(message = "Start model:")
    self.filter_solvent()
    if(self.filter_only == False):
       self.find_peaks()
       hd_sel = self.xray_structure.hd_selection()
       result_filter = filter_sites_and_map_next_to_model(
         xray_structure = self.xray_structure.select(~hd_sel),
         sites          = self.sites,
         max_dist       = self.max_solv_macromol_dist,
         min_dist       = self.min_solv_macromol_dist,
         log            = self.log)
       self.sites = result_filter.sites
       self.heights = self.heights.select(result_filter.selection)
       self.filter_close_peak_peak_contacts()
       self.create_solvent_xray_structure()
       self.update_xray_structure()
       self.filter_solvent()
    if(self.verbose > 0): self.show(message = "Final model:")

  def filter_solvent(self):
    hd_sel = self.xray_structure.hd_selection()
    hd_sel = hd_sel.select(~self.solvent_selection)
    xrs_sol = self.xray_structure.select(self.solvent_selection)
    xrs_mac_h = self.xray_structure.select(~self.solvent_selection)
    xrs_mac = xrs_mac_h.select(~hd_sel)
    selection = xrs_sol.all_selection()
    scat_sol = xrs_sol.scatterers()
    occ_sol = scat_sol.extract_occupancies()
    b_isos_sol = scat_sol.extract_u_iso_or_u_equiv(
      self.xray_structure.unit_cell()) * math.pi**2*8
    result = xrs_mac.closest_distances(sites_frac = xrs_sol.sites_frac(),
      distance_cutoff = self.max_solv_macromol_dist)
    selection &= b_isos_sol >= self.b_iso_min
    selection &= b_isos_sol <= self.b_iso_max
    selection &= occ_sol >= self.occupancy_min
    selection &= occ_sol <= self.occupancy_max
    selection &= result.smallest_distances <= self.max_solv_macromol_dist
    selection &= result.smallest_distances >= self.min_solv_macromol_dist
    xrs_sol = xrs_sol.select(selection)
    sol_sel = flex.bool(xrs_mac_h.scatterers().size(), False)
    sol_sel.extend( flex.bool(xrs_sol.scatterers().size(), True) )
    self.model.remove_solvent()
    self.model.add_solvent(
      solvent_selection      = sol_sel,
      solvent_xray_structure = xrs_sol,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id)
    self.xray_structure = self.model.xray_structure
    self.solvent_selection = self.model.solvent_selection

  def show(self, message):
    print >> self.log, message
    xrs_mac = self.xray_structure.select(~self.solvent_selection)
    xrs_sol = self.xray_structure.select(self.solvent_selection)
    scat = xrs_sol.scatterers()
    occ = scat.extract_occupancies()
    b_isos = scat.extract_u_iso_or_u_equiv(
      self.xray_structure.unit_cell()) * math.pi**2*8
    smallest_distances = xrs_mac.closest_distances(
      sites_frac      = xrs_sol.sites_frac(),
      distance_cutoff = self.max_solv_macromol_dist).smallest_distances
    number = format_value("%-7d",scat.size())
    b_min  = format_value("%-7.2f", flex.min_default( b_isos, None))
    b_max  = format_value("%-7.2f", flex.max_default( b_isos, None))
    b_ave  = format_value("%-7.2f", flex.mean_default(b_isos, None))
    bl_min = format_value("%-7.2f", self.b_iso_min).strip()
    bl_max = format_value("%-7.2f", self.b_iso_max).strip()
    o_min  = format_value("%-7.2f", flex.min_default(occ, None))
    o_max  = format_value("%-7.2f", flex.max_default(occ, None))
    ol_min = format_value("%-7.2f", self.occupancy_min).strip()
    ol_max = format_value("%-7.2f", self.occupancy_max).strip()
    d_min  = format_value("%-7.2f", flex.min_default(smallest_distances, None))
    d_max  = format_value("%-7.2f", flex.max_default(smallest_distances, None))
    dl_min = format_value("%-7.2f", self.min_solv_macromol_dist).strip()
    dl_max = format_value("%-7.2f", self.max_solv_macromol_dist).strip()
    print >> self.log,"  number           = %s"%number
    print >> self.log,"  b_iso_min        = %s (limit = %s)"%(b_min, bl_min)
    print >> self.log,"  b_iso_max        = %s (limit = %s)"%(b_max, bl_max)
    print >> self.log,"  b_iso_mean       = %s             "%(b_ave)
    print >> self.log,"  occupancy_min    = %s (limit = %s)"%(o_min, ol_min)
    print >> self.log,"  occupancy_max    = %s (limit = %s)"%(o_max, ol_max)
    print >> self.log,"  dist_sol_mol_min = %s (limit = %s)"%(d_min, dl_min)
    print >> self.log,"  dist_sol_mol_max = %s (limit = %s)"%(d_max, dl_max)

  def find_peaks(self):
    out = self.log
    self.sites, self.heights = self.find_peaks_helper(
      map_type = self.primary_map_type,
      k        = self.primary_map_k,
      n        = self.primary_map_n,
      peak_search_parameters = maptbx.peak_search_parameters(
        peak_search_level      = self.peak_search_level,
        max_peaks              = self.max_peaks,
        peak_cutoff            = self.primary_map_cutoff,
        interpolate            = self.interpolate,
        min_distance_sym_equiv = self.min_distance_sym_equiv,
        general_positions_only = self.general_positions_only,
        min_cross_distance     = self.min_cross_distance))
    if(self.verbose > 0):
       print >> out, "Peak search:"
       print >> out, "   maximum allowed number of peaks:          ", \
                                                       self.max_number_of_peaks
       print >> out, "   number of peaks from "+self.primary_map_type+\
                                                     " map: ",self.sites.size()
    if(self.secondary_map_type is not None):
       sites_2nd, heights_2nd = self.find_peaks_helper(
         map_type = self.secondary_map_type,
         k        = self.secondary_map_k,
         n        = self.secondary_map_n,
         peak_search_parameters = maptbx.peak_search_parameters(
           peak_search_level      = self.peak_search_level,
           max_peaks              = self.max_peaks,
           peak_cutoff            = self.secondary_map_cutoff,
           interpolate            = self.interpolate,
           min_distance_sym_equiv = self.min_distance_sym_equiv,
           general_positions_only = self.general_positions_only,
           min_cross_distance     = self.min_cross_distance))
       if(self.verbose > 0):
          print >> out, "   number of peaks from "+self.secondary_map_type+\
                                                      " map:", sites_2nd.size()
       clustering_manager = max_lik.peak_clustering(
                                     r1f    = self.sites,
                                     r2f    = sites_2nd,
                                     h1     = self.heights,
                                     h2     = heights_2nd,
                                     uc     = self.xray_structure.unit_cell(),
                                     cutoff = self.peak_map_matching_tolerance)
       self.sites = clustering_manager.sites()
       self.heights = clustering_manager.heights()
       if(self.verbose > 0):
          print >> out, "   number of peaks after clustering:         ", \
                                                              self.sites.size()

  def find_peaks_helper(self, map_type, k, n, peak_search_parameters):
    fft_map = self.fmodel.electron_density_map(
                           map_type          = map_type,
                           k                 = k,
                           n                 = n,
                           resolution_factor = self.resolution_factor,
                           symmetry_flags    = maptbx.use_space_group_symmetry)
    self.gridding_n_real = fft_map.n_real()
    if(self.use_sigma_scaled_maps): fft_map.apply_sigma_scaling()
    fft_map_data = fft_map.real_map_unpadded()
    crystal_gridding_tags = fft_map.tags()
    cluster_analysis = crystal_gridding_tags.peak_search(
        parameters = peak_search_parameters,
        map        = fft_map_data).all(max_clusters = self.max_number_of_peaks)
    sites = cluster_analysis.sites()
    heights = cluster_analysis.heights()
    return sites, heights

  def filter_close_peak_peak_contacts(self):
    out = self.log
    n_peaks_old = self.sites.size()
    sites = self.sites.deep_copy()
    heights = self.heights.deep_copy()
    k = 1
    for i,hi in zip(self.sites,self.heights):
        for j,hj in zip(self.sites,self.heights):
            if(abs(hi-hj) > 1.e-3):
               d = self.xray_structure.unit_cell().distance(i,j)
               if(d < self.min_solv_solv_dist and d > 1.e-3):
                  k = k + 1
                  if(hi > hj):
                     sel = heights == hj
                     heights = heights.select(~sel)
                     sites = sites.select(~sel)
                  else:
                     sel = heights == hi
                     heights = heights.select(~sel)
                     sites = sites.select(~sel)
    self.sites = sites
    self.heights = heights
    if(self.verbose > 0):
       n_peaks_new = n_peaks_old - self.sites.size()
       print >> out, "Peak filtering (peak - peak close contact elimination):"
       print >> out, "   peaks rejected:                 ", n_peaks_new
       print >> out, "   total number of peaks selected: ", self.sites.size()

  def create_solvent_xray_structure(self):
    if(self.b_iso is None):
       b = self.xray_structure.extract_u_iso_or_u_equiv() * math.pi**2*8
       b_solv = flex.mean_default(b, None)
       if(b_solv is not None and b_solv < self.b_iso_min or b_solv > self.b_iso_max):
          b_solv = (self.b_iso_min + self.b_iso_max) / 2.
    else:
       b_solv = self.b_iso
    new_scatterers = flex.xray_scatterer(self.sites.size(),
                                     xray.scatterer(
                                       occupancy = self.occupancy,
                                       b               = b_solv,
                                       scattering_type = self.scattering_type))
    new_scatterers.set_sites(self.sites)
    self.solvent_xray_structure = xray.structure(
                               special_position_settings = self.xray_structure,
                               scatterers                = new_scatterers)

  def update_xray_structure(self):
    xrs_sol = self.xray_structure.select(self.solvent_selection)
    xrs_mac = self.xray_structure.select(~self.solvent_selection)
    xrs_sol = xrs_sol.concatenate(other = self.solvent_xray_structure)
    sol_sel = flex.bool(xrs_mac.scatterers().size(), False)
    sol_sel.extend( flex.bool(xrs_sol.scatterers().size(), True) )
    self.model.remove_solvent()
    self.model.add_solvent(
      solvent_selection      = sol_sel,
      solvent_xray_structure = xrs_sol,
      residue_name           = self.params.output_residue_name,
      atom_name              = self.params.output_atom_name,
      chain_id               = self.params.output_chain_id)
    self.xray_structure = self.model.xray_structure
    self.solvent_selection = self.model.solvent_selection

class filter_sites_and_map_next_to_model(object):
  def __init__(self, xray_structure, sites, max_dist, min_dist, log):
    initial_number_of_sites = sites.size()
    print >> log, "Filter by distance & map next to the model:"
    result = xray_structure.closest_distances(sites_frac = sites,
      distance_cutoff = max_dist)
    smallest_distances_sq = result.smallest_distances_sq
    smallest_distances = result.smallest_distances
    new_sites = result.sites_frac
    not_too_far = smallest_distances_sq <= max_dist**2
    not_too_close = smallest_distances_sq >= min_dist**2
    self.selection = (not_too_far & not_too_close)
    self.sites = new_sites.select(self.selection)
    smallest_distances = flex.sqrt(smallest_distances_sq)
    d_min = flex.min_default(smallest_distances, 0)
    d_max = flex.max_default(smallest_distances, 0)
    print >> log, "   mapped sites are within: %5.3f - %5.3f " % (d_min, d_max)
    print >> log, "   number of sites selected in [dist_min=%5.2f, " \
      "dist_max=%5.2f]: %d from: %d" % (min_dist, max_dist, self.sites.size(),
      initial_number_of_sites)
    smallest_distances =flex.sqrt(smallest_distances_sq.select(self.selection))
    d_min = flex.min_default(smallest_distances, 0)
    d_max = flex.max_default(smallest_distances, 0)
    print >> log, "   mapped sites are within: %5.3f - %5.3f " % (d_min, d_max)

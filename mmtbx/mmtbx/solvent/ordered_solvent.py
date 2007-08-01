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

def show_histogram_(data,
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
    print
    print >> out, flex.max(data)
    print >> out, flex.min(data)
    print >> out, flex.mean(data)
    return histogram

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
    if(self.verbose > 0): self.show_current_state(header = "Start model:")
    self.check_existing_solvent()
    if(self.filter_only == False):
       self.find_peaks()
       if(solvent_selection.count(False) > 0):
          if(self.bulk_solvent_mask_exclusions):
             self.filter_peaks_with_bulk_solvent_mask()
       self.sites, self.heights, dummy = self.filter_by_distance(
                                 self.xray_structure, self.sites, self.heights)
       self.filter_close_peak_peak_contacts()
       self.create_solvent_xray_structure()
       self.update_xray_structure()
    if(self.filter_only == False):
       self.solvent_selection.extend(flex.bool(self.sites.size(), True))
    self.final_distance_check()
    if(self.verbose > 0): self.show_current_state(header = "Final model:")

  def final_distance_check(self):
    # Introduced as a quick fix for emerged problem.
    # Please carefull in making changes here: the effect of using of this
    # function is proven to be essential for selected cases but NOT exercised
    # in routine regression tests. This will be done (along with some code
    # cleaning) in the next revision of ordered_solvent.py code.
    print >> self.log, "*** WARNING from ordered solvent picking ***" # XXX
    print >> self.log, "*** WARNING Entering into suboptimal calculation..." # XXX
    print >> self.log, "*** WARNING Calculations may take some time..." # XXX
    print >> self.log
    sol_sel = self.solvent_selection
    mac_sel = ~self.solvent_selection
    xrs_sol = self.xray_structure.select(sol_sel)
    xrs_mac = self.xray_structure.select(mac_sel)
    sol_mac_distances = xrs_sol.closest_distances(other = xrs_mac,
      max_distance_cutoff = self.max_solv_macromol_dist)
    not_too_far = sol_mac_distances != -1
    not_too_close = sol_mac_distances >= self.min_solv_macromol_dist
    selection_good = (not_too_far & not_too_close)
    print >> self.log, \
      "Additional number of solvent removed by distance cutoffs = ", \
      selection_good.count(False)
    #
    xrs_sol = self.xray_structure.select(self.solvent_selection)
    xrs_mac = self.xray_structure.select(~self.solvent_selection)
    xrs_sol = xrs_sol.select(selection_good)
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

  def check_existing_solvent(self):
    scatterers      = self.xray_structure.scatterers()
    non_solvent_sel = self.solvent_selection.select(~self.solvent_selection)
    solvent_sel     = self.solvent_selection.select(self.solvent_selection)
    non_solvent_scatterers = scatterers.select(~self.solvent_selection)
    solvent_scatterers = scatterers.select(self.solvent_selection)
    occupancies = solvent_scatterers.extract_occupancies()
    b_isos = solvent_scatterers.extract_u_iso_or_u_equiv(
                                self.xray_structure.unit_cell()) * math.pi**2*8
    new_solvent_sel  = b_isos >= self.b_iso_min
    new_solvent_sel &= b_isos <= self.b_iso_max
    new_solvent_sel &= occupancies >= self.occupancy_min
    new_solvent_sel &= occupancies <= self.occupancy_max
    ############################
    xrs = self.xray_structure.deep_copy_scatterers()
    xrs.erase_scatterers()
    xrs.add_scatterers(non_solvent_scatterers)
    # XXX work-around to allow filling of an empty unit cell
    #if(solvent_sel.count(False) > 0):
    #   dummy, dummy, sel_by_dist = self.filter_by_distance(
    #                                   xrs, solvent_scatterers.extract_sites())
    #else:
    #   sel_by_dist = flex.bool(solvent_sel.size(), True)
    dummy, dummy, sel_by_dist = self.filter_by_distance(
                                       xrs, solvent_scatterers.extract_sites())
    ############################
    new_solvent_sel &= sel_by_dist
    ###
    reduce_original_sel = (~self.solvent_selection) | flex.bool(
                            size = scatterers.size(),
                            iselection = flex.size_t(xrange(scatterers.size()))
                                         .select(self.solvent_selection)
                                         .select(new_solvent_sel))
    ###
    new_solvent_scatterers = solvent_scatterers.select(new_solvent_sel)
    non_solvent_sel.extend(new_solvent_sel.select(new_solvent_sel))
    self.solvent_selection = non_solvent_sel
    self.xray_structure.erase_scatterers()
    self.xray_structure.add_scatterers(non_solvent_scatterers)
    self.xray_structure.add_scatterers(new_solvent_scatterers)
    if(self.verbose > 0):
       st = "b_iso = [%.2f,%.2f] & q = [%.2f,%.2f]:"% \
                                               (self.b_iso_min,self.b_iso_max,\
                                         self.occupancy_min,self.occupancy_max)
       self.show_current_state(header= "Initial solvent selection with "+st)
    ########

    self.model.update(selection = reduce_original_sel)
    #############
    self.xray_structure = self.model.xray_structure.deep_copy_scatterers()
    self.solvent_selection = self.model.solvent_selection.deep_copy()
    #################
    assert approx_equal(self.solvent_selection, self.model.solvent_selection)

  def show_current_state(self, header):
    out = self.log
    scatterers = self.xray_structure.scatterers()
    non_solvent_scatterers = scatterers.select(~self.solvent_selection)
    solvent_scatterers     = scatterers.select(self.solvent_selection)
    q_solv = solvent_scatterers.extract_occupancies()
    q_prot = non_solvent_scatterers.extract_occupancies()
    b_solv = solvent_scatterers.extract_u_iso_or_u_equiv(
                                self.xray_structure.unit_cell()) * math.pi**2*8
    b_prot = non_solvent_scatterers.extract_u_iso_or_u_equiv(
                                self.xray_structure.unit_cell()) * math.pi**2*8
    n_solv = self.solvent_selection.count(True)
    n_prot = self.solvent_selection.count(False)
    n_tot  = self.xray_structure.scatterers().size()
    st1 = (n_solv,n_prot,n_tot)
    if(n_solv > 0):
       st2 = (flex.max(b_solv),flex.max_default(b_prot, None))
       st3 = (flex.min(b_solv),flex.min_default(b_prot, None))
       st4 = (flex.mean(b_solv),flex.mean_default(b_prot, None))
       st5 = (flex.max(q_solv),flex.min_default(q_prot, None))
       st6 = (flex.min(q_solv),flex.max_default(q_prot, None))
    else:
       st2 = (0,flex.max_default(b_prot, None))
       st3 = (0,flex.min_default(b_prot, None))
       st4 = (0,flex.mean_default(b_prot, None))
       st5 = (0,flex.min_default(q_prot, None))
       st6 = (0,flex.max_default(q_prot, None))
    print >> out, header
    print >> out, "                   solvent non-solvent       total"
    print >> out, "   number    =%12d%12d%12d" % st1
    try: print >> out, "   b_iso_max =%12.2f%12.2f" % st2
    except: print >> out, "   b_iso_max =%12.2f%12s" % st2
    try: print >> out, "   b_iso_min =%12.2f%12.2f" % st3
    except: print >> out, "   b_iso_min =%12.2f%12s" % st3
    try: print >> out, "   b_iso_ave =%12.2f%12.2f" % st4
    except: print >> out, "   b_iso_ave =%12.2f%12s" % st4
    try: print >> out, "   q_min     =%12.2f%12.2f" % st5
    except: print >> out, "   q_min     =%12.2f%12s" % st5
    try: print >> out, "   q_max     =%12.2f%12.2f" % st6
    except: print >> out, "   q_max     =%12.2f%12s" % st6

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
    self.tags = crystal_gridding_tags.tags()
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

  def filter_peaks_with_bulk_solvent_mask(self):
    out = self.log
    print >> out, "Peak filtering (bulk solvent mask exclusions):"
    bulk_solvent_mask = masks.bulk_solvent(
                               xray_structure           = self.xray_structure,
                               gridding_n_real          = self.gridding_n_real,
                               solvent_radius           = 0.0,
                               shrink_truncation_radius = 0.0)
    mask_data = bulk_solvent_mask.data.deep_copy()
    print >> out, "   initial mask:"
    print >> out, "      number of 0:     %d" % mask_data.count(1)
    print >> out, "      number of 1:     %d" % mask_data.count(0)
    print >> out, "      solvent fraction: %6.3f" % \
                                     bulk_solvent_mask.contact_surface_fraction
    n_overlaps = self.tags.apply_symmetry_to_mask(mask_data)
    print >> out, "   symmetrized mask:"
    print >> out, "      number of 0:     %d" % mask_data.count(1)
    print >> out, "      number of 1:     %d" % mask_data.count(0)
    print >> out, "      solvent fraction: %6.3f"% (float(mask_data.count(1))/\
                                    matrix.col(self.gridding_n_real).product())
    print >> out, "   number of overlaps in au = ", n_overlaps
    mask_data_double =  mask_data.as_double()
    bs = 0
    sites = flex.vec3_double()
    heights = flex.double()
    for site, height in zip(self.sites, self.heights):
        v = maptbx.value_at_closest_grid_point(mask_data_double, site)
        if(v == 0.0): bs += 1
        if(v > 0.0):
           sites.append(site)
           heights.append(height)
    self.sites = sites
    self.heights = heights
    if(self.verbose > 0):
       print >> out, "   peaks rejected (inside macromolecule region):",bs
       print >> out, "   total number of peaks selected:              ",\
                                                              self.sites.size()

  def filter_by_distance(self, xray_structure, sites, heights = None):
    initial_number_of_sites = sites.size()
    out = self.log
    print >> out, "Peak filtering by distance & mapping next to the model:"
    distance_cutoff = self.max_solv_macromol_dist
    asu_mappings = xray_structure.asu_mappings(buffer_thickness = \
                                                               distance_cutoff)
    asu_mappings.process_sites_frac(sites, min_distance_sym_equiv = \
                                  xray_structure.min_distance_sym_equiv())
    pair_generator = crystal.neighbors_fast_pair_generator(
                                             asu_mappings    = asu_mappings,
                                             distance_cutoff = distance_cutoff)
    n_xray = xray_structure.scatterers().size()
    water_next_to_protein_sites = sites.deep_copy()
    if(heights is not None):
       height_next_to_protein_sites = heights.deep_copy()
    smallest_distances_sq = flex.double(sites.size(), distance_cutoff**2+1.)
    for pair in pair_generator:
        if (pair.i_seq < n_xray):
          if (pair.j_seq < n_xray): continue
          # i_seq = protein
          # j_seq = water
          rt_mx_i = asu_mappings.get_rt_mx_i(pair)
          rt_mx_j = asu_mappings.get_rt_mx_j(pair)
          rt_mx_ji = rt_mx_i.inverse().multiply(rt_mx_j)
          i_seq_water = pair.j_seq - n_xray
          water_site = rt_mx_ji * sites[i_seq_water]
          if(heights is not None):
             water_height = heights[i_seq_water]
        else:
          if (pair.j_seq >= n_xray): continue
          # i_seq = water
          # j_seq = protein
          rt_mx_i = asu_mappings.get_rt_mx_i(pair)
          rt_mx_j = asu_mappings.get_rt_mx_j(pair)
          rt_mx_ij = rt_mx_j.inverse().multiply(rt_mx_i)
          i_seq_water = pair.i_seq - n_xray
          water_site = rt_mx_ij * sites[i_seq_water]
          if(heights is not None):
             water_height = heights[i_seq_water]
        if (smallest_distances_sq[i_seq_water] > pair.dist_sq):
          smallest_distances_sq[i_seq_water] = pair.dist_sq
          water_next_to_protein_sites[i_seq_water] = water_site
          if(heights is not None):
             height_next_to_protein_sites[i_seq_water] = water_height
    not_too_far = smallest_distances_sq <= distance_cutoff**2
    not_too_close = smallest_distances_sq >= self.min_solv_macromol_dist**2
    selection = (not_too_far & not_too_close)
    sites = water_next_to_protein_sites.select(selection)
    if(heights is not None):
       heights = height_next_to_protein_sites.select(selection)
    smallest_distances = flex.sqrt(smallest_distances_sq)
    d_min = flex.min_default(smallest_distances, 0)
    d_max = flex.max_default(smallest_distances, 0)
    print >> out, "   mapped sites are within: %5.3f - %5.3f " % (d_min, d_max)
    print >> out, "   number of peaks selected in [dist_min=%5.2f, " \
          "dist_max=%5.2f]: %d from: %d" % (self.min_solv_macromol_dist,\
          distance_cutoff, sites.size(), initial_number_of_sites)
    smallest_distances = flex.sqrt(smallest_distances_sq.select(selection))
    d_min = flex.min_default(smallest_distances, 0)
    d_max = flex.max_default(smallest_distances, 0)
    print >> out, "   mapped sites are within: %5.3f - %5.3f " % (d_min, d_max)
    return sites, heights, selection

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
    self.model.add_solvent(
                      solvent_selection      = self.solvent_selection,
                      solvent_xray_structure = self.solvent_xray_structure,
                      residue_name           = self.params.output_residue_name,
                      atom_name              = self.params.output_atom_name,
                      chain_id               = self.params.output_chain_id)
    self.xray_structure = self.model.xray_structure
    self.solvent_selection = self.model.solvent_selection

def show_histogram(data = None,
                   n_slots = None):
    histogram = flex.histogram(data    = data,
                               n_slots = n_slots)
    low_cutoff = histogram.data_min()
    for (i,n) in enumerate(histogram.slots()):
      high_cutoff = histogram.data_min() + histogram.slot_width() * (i+1)
      print "%20.6f - %20.6f: %5d" % \
             (low_cutoff, high_cutoff, n)
      low_cutoff = high_cutoff

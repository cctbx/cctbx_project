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


def get_manager(fmodel, model, solvent_selection, log = None, params = None):
  if(params is not None):
     return manager(
       fmodel                             = fmodel,
       model                              = model,
       solvent_selection                  = solvent_selection,
       b_iso_min                          = params.b_iso_min,
       b_iso_max                          = params.b_iso_max,
       b_iso                              = params.b_iso,
       scattering_type                    = params.scattering_type,
       occupancy_min                      = params.occupancy_min,
       occupancy_max                      = params.occupancy_max,
       occupancy                          = params.occupancy,
       use_sigma_scaled_maps              = params.use_sigma_scaled_maps,
       primary_map_type                   = params.primary_map_type,
       primary_map_k                      = params.primary_map_k,
       primary_map_n                      = params.primary_map_n,
       primary_map_cutoff                 = params.primary_map_cutoff,
       secondary_map_type                 = params.secondary_map_type,
       secondary_map_k                    = params.secondary_map_k,
       secondary_map_n                    = params.secondary_map_n,
       secondary_map_cutoff               = params.secondary_map_cutoff,
       peak_map_matching_tolerance        = params.peak_map_matching_tolerance,
       resolution_factor                  = params.resolution_factor,
       min_solv_macromol_dist             = params.min_solv_macromol_dist,
       max_solv_macromol_dist             = params.max_solv_macromol_dist,
       max_number_of_peaks                = params.max_number_of_peaks,
       solvent_pdb_file_name              = None,
       filter_only                        = params.filter_only,
       verbose                            = params.verbose,
       log                                = log,
       peak_search_level           = params.peak_search.peak_search_level,
       max_peaks                   = params.peak_search.max_peaks,
       interpolate                 = params.peak_search.interpolate,
       min_distance_sym_equiv      = params.peak_search.min_distance_sym_equiv,
       general_positions_only      = params.peak_search.general_positions_only,
       min_cross_distance          = params.peak_search.min_cross_distance,
       params                      = params)
  else:
     return manager(fmodel, solvent_selection, params = params)

class manager(object):
  def __init__(self, fmodel,
                     model,
                     solvent_selection,
                     b_iso_min                   = 1.0,
                     b_iso_max                   = 50.0,
                     b_iso                       = None,
                     scattering_type             = "O",
                     occupancy_min               = 0.1,
                     occupancy_max               = 2.0,
                     occupancy                   = 1.0,
                     use_sigma_scaled_maps       = True,
                     primary_map_type            = "m*Fobs-D*Fmodel",
                     primary_map_k               = None,
                     primary_map_n               = None,
                     primary_map_cutoff          = 3.0,
                     secondary_map_type          = "2m*Fobs-D*Fmodel",
                     secondary_map_k             = None,
                     secondary_map_n             = None,
                     secondary_map_cutoff        = 1.0,
                     peak_map_matching_tolerance = 2.0,
                     resolution_factor           = 1./4.,
                     min_solv_macromol_dist      = 1.8,
                     max_solv_macromol_dist      = 6.0,
                     min_solv_solv_dist          = 1.8,
                     max_number_of_peaks         = None,
                     solvent_pdb_file_name       = None,
                     filter_only                 = False,
                     verbose                     = 1,
                     log                         = None,
                     peak_search_level           = 1,
                     max_peaks                   = 0,
                     peak_cutoff                 = 1.0,
                     interpolate                 = True,
                     min_distance_sym_equiv      = 1.e-6,
                     general_positions_only      = False,
                     min_cross_distance          = 2.0,
                     params                      = None):
    adopt_init_args(self, locals())
    if (self.log is None): self.log = sys.stdout
    self.xray_structure = self.fmodel.xray_structure.deep_copy_scatterers()
    self.sites = None
    self.heights = None
    if(self.max_number_of_peaks is None):
       self.max_number_of_peaks= int(self.solvent_selection.count(False)/10*10)
    if(self.verbose > 0): self.show_current_state(header = "Start model:")
    self.check_existing_solvent()
    if(self.filter_only == False):
       self.find_peaks()
       self.filter_peaks_with_bulk_solvent_mask()
       self.sites, self.heights, dummy = self.filter_by_distance(
                                 self.xray_structure, self.sites, self.heights)
       self.filter_close_peak_peak_contacts()
       self.create_solvent_xray_structure()
       self.update_xray_structure()
    if(self.solvent_pdb_file_name is not None):
       open(self.solvet_pdb_file_name, "w").write(self.xray_structure.as_pdb_file())
    if(self.filter_only == False):
       self.solvent_selection.extend(flex.bool(self.sites.size(), True))
    if(self.verbose > 0): self.show_current_state(header = "Final model:")

  def check_existing_solvent(self):
    #for sc, s in zip(self.xray_structure.scatterers(), self.model.solvent_selection):
    #  print sc.element_symbol(), s
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

    ###
    #for sc1, s1, sc2, s2, aal in zip(self.xray_structure.scatterers(), self.solvent_selection,
    #         self.model.xray_structure.scatterers(), self.model.solvent_selection, self.model.atom_attributes_list):
    #    print sc1.element_symbol(), s1, sc2.element_symbol(), s2, aal

    assert approx_equal(self.solvent_selection, self.model.solvent_selection)

  def show_current_state(self, header):
    out = self.log
    scatterers = self.xray_structure.scatterers()
    print scatterers.size(), self.solvent_selection.size()
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
       st2 = (flex.max(b_solv),flex.max(b_prot))
       st3 = (flex.min(b_solv),flex.min(b_prot))
       st4 = (flex.mean(b_solv),flex.mean(b_prot))
       st5 = (flex.max(q_solv),flex.min(q_prot))
       st6 = (flex.min(q_solv),flex.max(q_prot))
    else:
       st2 = (0,flex.max(b_prot))
       st3 = (0,flex.min(b_prot))
       st4 = (0,flex.mean(b_prot))
       st5 = (0,flex.min(q_prot))
       st6 = (0,flex.max(q_prot))
    print >> out, header
    print >> out, "                   solvent non-solvent       total"
    print >> out, "   number    =%12d%12d%12d" % st1
    print >> out, "   b_iso_max =%12.2f%12.2f" % st2
    print >> out, "   b_iso_min =%12.2f%12.2f" % st3
    print >> out, "   b_iso_ave =%12.2f%12.2f" % st4
    print >> out, "   q_min     =%12.2f%12.2f" % st5
    print >> out, "   q_max     =%12.2f%12.2f" % st6

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
       b_solv = flex.mean(b)
       if(b_solv < self.b_iso_min or b_solv > self.b_iso_max):
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

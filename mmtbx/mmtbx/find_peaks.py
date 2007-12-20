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

master_params = iotbx.phil.parse("""\
  use_sigma_scaled_maps = True
    .type=bool
    .help = Use sigma scales maps for water pick picking
  map_type = mFobs-DFmodel
    .type=str
  map_cutoff = 3.0
    .type=float
  resolution_factor = 1./4.
    .type=float
  map_next_to_model
  {
    min_model_peak_dist = 1.8
      .type=float
    max_model_peak_dist = 6.0
      .type=float
    min_peak_peak_dist = 1.8
      .type=float
    use_hydrogens = False
      .type = bool
  }
  max_number_of_peaks = None
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

class peaks_holder(object):
  def __init__(self, heights, sites):
    self.heights = heights
    self.sites = sites

class manager(object):
  def __init__(self, fmodel, params = None, log = None):
    adopt_init_args(self, locals())
    if(self.log is None): self.log = sys.stdout
    if(self.params is None): self.params = master_params.extract()
    peak_search_parameters = maptbx.peak_search_parameters(
      peak_search_level      = self.params.peak_search.peak_search_level,
      max_peaks              = self.params.peak_search.max_peaks,
      peak_cutoff            = self.params.map_cutoff,
      interpolate            = self.params.peak_search.interpolate,
      min_distance_sym_equiv = self.params.peak_search.min_distance_sym_equiv,
      general_positions_only = self.params.peak_search.general_positions_only,
      min_cross_distance     = self.params.peak_search.min_cross_distance)
    fft_map = self.fmodel.electron_density_map(
      map_type          = self.params.map_type,
      resolution_factor = self.params.resolution_factor,
      symmetry_flags    = maptbx.use_space_group_symmetry)
    gridding_n_real = fft_map.n_real()
    map_units = "e/A**3"
    if(self.params.use_sigma_scaled_maps):
      fft_map.apply_sigma_scaling()
      map_units = "sigma"
    fft_map_data = fft_map.real_map_unpadded()
    crystal_gridding_tags = fft_map.tags()
    cluster_analysis = crystal_gridding_tags.peak_search(
      parameters = peak_search_parameters,
      map = fft_map_data).all(max_clusters = self.params.max_number_of_peaks)
    self.peaks_ = peaks_holder(heights = cluster_analysis.heights(),
                               sites   = cluster_analysis.sites())
    print >>self.log,"Number of peaks found at %s map (map cutoff=%s %s)= %s"%(
      self.params.map_type, format_value("%-5.2f", self.params.map_cutoff),
      map_units, format_value("%-12d", self.peaks_.sites.size()))

  def peaks(self):
    return self.peaks_

  def peaks_mapped(self):
   max_dist = self.params.map_next_to_model.max_model_peak_dist
   min_dist = self.params.map_next_to_model.min_model_peak_dist
   xray_structure = self.fmodel.xray_structure.deep_copy_scatterers()
   if(not self.params.map_next_to_model.use_hydrogens):
     hd_sel = xray_structure.hd_selection()
     xray_structure = xray_structure.select(~hd_sel)
   initial_number_of_sites = self.peaks_.sites.size()
   print >> self.log, "Filter by distance & map next to the model:"
   result = xray_structure.closest_distances(sites_frac = self.peaks_.sites,
     distance_cutoff = max_dist)
   smallest_distances_sq = result.smallest_distances_sq
   smallest_distances = result.smallest_distances
   new_sites = result.sites_frac
   new_heights = self.peaks_.heights.select(~result.remove_selection)
   in_box = smallest_distances_sq > 0
   not_too_far = smallest_distances_sq <= max_dist**2
   not_too_close = smallest_distances_sq >= min_dist**2
   selection = (not_too_far & not_too_close & in_box)
   peaks = peaks_holder(heights = new_heights.select(selection),
                        sites   = new_sites.select(selection))
   sd = flex.sqrt(smallest_distances_sq.select(in_box))
   d_min = flex.min_default(sd, 0)
   d_max = flex.max_default(sd, 0)
   print >> self.log, "   mapped sites are within: %5.3f - %5.3f"%(d_min,d_max)
   print >> self.log, "   number of sites selected in [dist_min=%5.2f, " \
     "dist_max=%5.2f]: %d from: %d" % (min_dist, max_dist, peaks.sites.size(),
     initial_number_of_sites)
   smallest_distances = flex.sqrt(smallest_distances_sq.select(selection))
   d_min = flex.min_default(smallest_distances, 0)
   d_max = flex.max_default(smallest_distances, 0)
   print >> self.log, "   mapped sites are within: %5.3f - %5.3f"%(d_min,d_max)
   return peaks

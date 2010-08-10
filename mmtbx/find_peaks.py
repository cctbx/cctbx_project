from cctbx.array_family import flex
import math, sys
from libtbx import adopt_init_args
from cctbx import maptbx
from cctbx import crystal
import libtbx.load_env
from libtbx.test_utils import approx_equal
import iotbx.phil
from libtbx.str_utils import format_value
from mmtbx.refinement import print_statistics

master_params = iotbx.phil.parse("""\
  use_sigma_scaled_maps = True
    .type=bool
    .help = Default is sigma scaled map, map in absolute scale is used \
            otherwise.
  resolution_factor = 1./4.
    .type=float
  map_next_to_model
    .expert_level=2
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
    .expert_level=1
  peak_search
    .expert_level=1
  {
    peak_search_level = 1
      .type=int
    max_peaks = 0
      .type=int
      .short_caption=Maximum peaks
    interpolate = True
      .type=bool
    min_distance_sym_equiv = None
      .type=float
      .short_caption=Minimum distance between symmetry-equivalent atoms
    general_positions_only = False
      .type=bool
    min_cross_distance = 1.8
      .type=float
      .short_caption=Minimum cross distance
    min_cubicle_edge = 5
      .type=float
      .short_caption=Minimum edge length of cubicles used for \
        fast neighbor search
      .expert_level=2
  }
""")

class peaks_holder(object):
  def __init__(self, heights, sites, iseqs_of_closest_atoms = None):
    self.heights = heights
    self.sites = sites
    self.iseqs_of_closest_atoms = iseqs_of_closest_atoms

class manager(object):
  def __init__(self, fmodel, map_type, map_cutoff, params = None, log = None,
                     use_kick_map = False, kick_map_params = None,
                     silent = False):
    adopt_init_args(self, locals())
    self.mapped = False
    self.peaks_ = None
    if(self.log is None): self.log = sys.stdout
    if(self.params is None): self.params = master_params.extract()
    if(use_kick_map):
      from mmtbx import map_tools
      km = map_tools.kick_map(
        fmodel            = self.fmodel,
        map_type          = map_type)
      fft_map = km.map_coeffs.fft_map(
        resolution_factor=self.params.resolution_factor,
        symmetry_flags=maptbx.use_space_group_symmetry)
      fft_map.apply_sigma_scaling()
      fft_map_data = fft_map.real_map_unpadded()
      map_units = "sigma"
    else:
      fft_map = self.fmodel.electron_density_map().\
        fft_map(resolution_factor = self.params.resolution_factor,
                symmetry_flags    = maptbx.use_space_group_symmetry,
                map_type          = self.map_type)
      if(self.params.use_sigma_scaled_maps):
        fft_map.apply_sigma_scaling()
        map_units = "sigma"
      else:
        fft_map.apply_volume_scaling()
        map_units = "e/A**3"
      fft_map_data = fft_map.real_map_unpadded()
    gridding_n_real = fft_map.n_real()
    crystal_gridding_tags = fft_map.tags()
    max_number_of_peaks = self.params.max_number_of_peaks
    if(self.params.max_number_of_peaks is None):
      max_number_of_peaks = self.fmodel.xray_structure.scatterers().size() * 5
    negative = False
    if(self.map_cutoff < 0):
      self.map_cutoff *= -1
      negative = True
      fft_map_data = fft_map_data * (-1.)
    min_distance_sym_equiv = self.params.peak_search.min_distance_sym_equiv
    if(min_distance_sym_equiv is None):
      min_distance_sym_equiv = \
        self.fmodel.xray_structure.min_distance_sym_equiv()
    peak_search_parameters = maptbx.peak_search_parameters(
      peak_search_level      = self.params.peak_search.peak_search_level,
      max_peaks              = self.params.peak_search.max_peaks,
      peak_cutoff            = self.map_cutoff,
      interpolate            = self.params.peak_search.interpolate,
      min_distance_sym_equiv = min_distance_sym_equiv,
      general_positions_only = self.params.peak_search.general_positions_only,
      min_cross_distance     = self.params.peak_search.min_cross_distance,
      min_cubicle_edge       = self.params.peak_search.min_cubicle_edge)
    if(self.fmodel.r_work() > 0.00001 and self.fmodel.r_free() > 0.00001):
      cluster_analysis = crystal_gridding_tags.peak_search(
        parameters = peak_search_parameters,
        map = fft_map_data).all(max_clusters = max_number_of_peaks)
      heights = cluster_analysis.heights()
      if(negative):
        heights *= -1.
      self.peaks_ = peaks_holder(heights = heights,
                                 sites   = cluster_analysis.sites())
      if(not self.silent):
        print >>self.log,"Number of peaks found at %s map (map cutoff=%s %s)= %s"%(
          self.map_type, format_value("%-5.2f", self.map_cutoff).strip(),
          map_units, format_value("%-12d", self.peaks_.sites.size()))

  def peaks(self):
    return self.peaks_

  def peaks_mapped(self):
    if(self.peaks_ is None): return None
    assert self.mapped == False
    max_dist = self.params.map_next_to_model.max_model_peak_dist
    min_dist = self.params.map_next_to_model.min_model_peak_dist
    xray_structure = self.fmodel.xray_structure.deep_copy_scatterers()
    use_selection = None
    if(not self.params.map_next_to_model.use_hydrogens):
      use_selection = ~xray_structure.hd_selection()
    initial_number_of_sites = self.peaks_.sites.size()
    if(not self.silent):
      print >> self.log, "Filter by distance & map next to the model:"
    result = xray_structure.closest_distances(sites_frac = self.peaks_.sites,
      distance_cutoff = max_dist, use_selection = use_selection)
    smallest_distances_sq = result.smallest_distances_sq
    smallest_distances = result.smallest_distances
    in_box = smallest_distances_sq > 0
    not_too_far = smallest_distances_sq <= max_dist**2
    not_too_close = smallest_distances_sq >= min_dist**2
    selection = (not_too_far & not_too_close & in_box)
    iseqs_of_closest_atoms = result.i_seqs.select(selection)
    peaks = peaks_holder(
      heights                = self.peaks_.heights.select(selection),
      sites                  = result.sites_frac.select(selection),
      iseqs_of_closest_atoms = iseqs_of_closest_atoms)
    sd = flex.sqrt(smallest_distances_sq.select(in_box))
    d_min = flex.min_default(sd, 0)
    d_max = flex.max_default(sd, 0)
    if(not self.silent):
      print >> self.log,"   mapped sites are within: %5.3f - %5.3f"%(d_min,d_max)
      print >> self.log, "   number of sites selected in [dist_min=%5.2f, " \
        "dist_max=%5.2f]: %d from: %d" % (min_dist, max_dist, peaks.sites.size(),
        initial_number_of_sites)
    smallest_distances = flex.sqrt(smallest_distances_sq.select(selection))
    d_min = flex.min_default(smallest_distances, 0)
    d_max = flex.max_default(smallest_distances, 0)
    if(not self.silent):
      print >> self.log,"   mapped sites are within: %5.3f - %5.3f"%(d_min,d_max)
    self.mapped = True
    self.peaks_ = peaks
    return peaks

  def show_mapped(self, pdb_atoms):
    if(self.peaks_ is None): return None
    peaks = self.peaks()
    if(peaks.iseqs_of_closest_atoms is None):
      raise RuntimeError("iseqs_of_closest_atoms is None")
    scatterers = self.fmodel.xray_structure.scatterers()
    assert scatterers.size() == pdb_atoms.size()
    assert peaks.sites.size() == peaks.heights.size()
    assert peaks.heights.size() == peaks.iseqs_of_closest_atoms.size()
    print >> self.log
    dist = self.fmodel.xray_structure.unit_cell().distance
    for i in flex.sort_permutation(data=peaks.iseqs_of_closest_atoms):
      s = peaks.sites[i]
      h = peaks.heights[i]
      i_seq = peaks.iseqs_of_closest_atoms[i]
      sc = scatterers[i_seq]
      d = dist(s, sc.site)
      element = sc.element_symbol()
      print >> self.log, "peak= %8.3f closest distance to %s = %8.3f" % (
        h, pdb_atoms[i_seq].id_str(), d)
      assert d <= self.params.map_next_to_model.max_model_peak_dist
      assert d >= self.params.map_next_to_model.min_model_peak_dist

def show_highest_peaks_and_deepest_holes(fmodel,
                                         pdb_atoms,
                                         map_type,
                                         map_cutoff_plus,
                                         map_cutoff_minus,
                                         log = None):
  if(log is None): log = sys.stdout
  print_statistics.make_header(
    "residual map %s: highest peaks and deepst holes"%map_type, out = log)
  fp_params = master_params.extract()
  fp_params.map_next_to_model.min_model_peak_dist = 0.7
  fp_params.map_next_to_model.min_peak_peak_dist = 0.7
  fp_params.peak_search.min_cross_distance = 0.7
  fp_params.map_next_to_model.use_hydrogens = True
  for par in [(map_cutoff_plus,"peaks"), (map_cutoff_minus,"holes")]:
    print_statistics.make_sub_header(par[1], out = log)
    result = manager(fmodel     = fmodel,
                     map_type   = "mFobs-DFmodel",
                     map_cutoff = par[0],
                     params     = fp_params,
                     log        = log)
    result.peaks_mapped()
    result.show_mapped(pdb_atoms = pdb_atoms)

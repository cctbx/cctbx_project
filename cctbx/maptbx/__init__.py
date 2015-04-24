from __future__ import division
import cctbx.sgtbx

import boost.python
ext = boost.python.import_ext("cctbx_maptbx_ext")
from cctbx_maptbx_ext import *

from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx import matrix
from scitbx.python_utils import dicts
from libtbx import adopt_init_args
import libtbx.load_env
import math
import sys, os
import scitbx.math

debug_peak_cluster_analysis = os.environ.get(
  "CCTBX_MAPTBX_DEBUG_PEAK_CLUSTER_ANALYSIS", "")

def value_at_closest_grid_point(map, x_frac):
  return map[closest_grid_point(map.accessor(), x_frac)]

flex.int.value_at_closest_grid_point = value_at_closest_grid_point
flex.double.value_at_closest_grid_point = value_at_closest_grid_point
flex.double.eight_point_interpolation = eight_point_interpolation
flex.double.tricubic_interpolation = tricubic_interpolation

def cc_peak(cutoff, map_1=None,map_2=None, map_coeffs_1=None,map_coeffs_2=None):
  """
  Compute CCpeak as described in
    Acta Cryst. (2014). D70, 2593-2606
    Metrics for comparison of crystallographic maps
    A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams
  """
  from cctbx import miller
  assert [map_1,map_2].count(None) in [0,2]
  assert [map_coeffs_1,map_coeffs_2].count(None) in [0,2]
  if([map_1,map_2].count(None)==0):
    # Maps are assumed to be quantile rank scaled (HE).
    return ext.cc_peak(map_1=map_1, map_2=map_2, cutoff=cutoff)
  elif([map_coeffs_1,map_coeffs_2].count(None)==0):
    d_min = min(map_coeffs_1.d_min(), map_coeffs_2.d_min())
    crystal_gridding = map_coeffs_1.crystal_gridding(
      d_min             = d_min,
      symmetry_flags    = use_space_group_symmetry,
      resolution_factor = 0.25)
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = map_coeffs_1)
    map_1 = fft_map.real_map_unpadded()
    fft_map = miller.fft_map(
      crystal_gridding     = crystal_gridding,
      fourier_coefficients = map_coeffs_2)
    map_2 = fft_map.real_map_unpadded()
    m1_he = volume_scale(map = map_1,  n_bins = 10000).map_data()
    m2_he = volume_scale(map = map_2,  n_bins = 10000).map_data()
    return ext.cc_peak(map_1=m1_he, map_2=m2_he, cutoff=cutoff)
  else:
    raise RuntimeError("Combination of inputs not supported.")

def map_accumulator(n_real, use_max_map, smearing_b=5, max_peak_scale=2,
                    smearing_span=10, use_exp_table=True):
  """
  Good defaults for 2mFo-DFc type maps:
    smearing_b=1, max_peak_scale=100, smearing_span=5
  """
  return ext.map_accumulator(n_real=n_real, smearing_b=smearing_b,
    max_peak_scale=max_peak_scale, smearing_span=smearing_span,
      use_exp_table=use_exp_table, use_max_map=use_max_map)

def peak_volume_estimate(map_data, sites_cart, crystal_symmetry, cutoff,
      atom_radius=1.5):
  v = flex.double()
  sites_frac = crystal_symmetry.unit_cell().fractionalize(sites_cart)
  for sc, sf in zip(sites_cart, sites_frac):
    if(map_data.value_at_closest_grid_point(sf)>=cutoff):
      sel = grid_indices_around_sites(
        unit_cell  = crystal_symmetry.unit_cell(),
        fft_n_real = map_data.focus(),
        fft_m_real = map_data.all(),
        sites_cart = flex.vec3_double([sc]),
        site_radii = flex.double([atom_radius]*1))
      v.append((map_data.select(sel)>=cutoff).count(True))
  r = flex.min_default(v, None)
  if(r==0): return None
  return r

def truncate(map_data, by_sigma_less_than, scale_by, set_value=0):
  """
  Trunate map inplace by standard deviation (sigma) while scale it with
  specified scale, such as volume (scale_by=1/volume) or sigma
  (scale_by=1/standard_deviation). Input map_data is expected to be unscaled (
  right out of FT).
  """
  sigma = statistics(map_data).sigma()
  if(sigma == 0):
    map_data = map_data*scale_by
    return
  ext.truncate(
    map_data           = map_data,
    standard_deviation = sigma,
    by_sigma_less_than = by_sigma_less_than,
    scale_by           = scale_by,
    set_value          = set_value)

def mask(xray_structure, n_real, solvent_radius):
  xrs_p1 = xray_structure.expand_to_p1(sites_mod_positive=True)
  from cctbx.masks import vdw_radii_from_xray_structure
  radii = vdw_radii_from_xray_structure(xray_structure = xrs_p1)
  radii = radii + solvent_radius
  return ext.mask(
    sites_frac = xrs_p1.sites_frac(),
    unit_cell  = xrs_p1.unit_cell(),
    n_real     = n_real,
    radii      = radii)

class statistics(ext.statistics):

  def __init__(self, map):
    ext.statistics.__init__(self, map)

class _(boost.python.injector, ext.statistics):

  def show_summary(self, f=None, prefix=""):
    if (f is None): f = sys.stdout
    print >> f, prefix + "max %.6g" % (self.max())
    print >> f, prefix + "min %.6g" % (self.min())
    print >> f, prefix + "mean %.6g" % (self.mean())
    print >> f, prefix + "sigma %.6g" % (self.sigma())

use_space_group_symmetry = sgtbx.search_symmetry_flags(
  use_space_group_symmetry=True)

class _(boost.python.injector, ext.histogram) :
  """
  Injector for extending cctbx.maptbx.histogram
  """
  # XXX make a method of scitbx
  def get_percentile_cutoffs (self, map, vol_cutoff_plus_percent,
      vol_cutoff_minus_percent) :
    """
    For the double-step filtration in cctbx.miller (used as part of the
    procedure for replacing missing F-obs in maps), we need to calculate upper
    and lower cutoffs for the data based on percentile values.  This can be
    done in just a few lines of code by using flex.sort_permutation over the
    entire map, but this has a huge memory overhead (and possibly computational
    overhead as well).  Since we are only interested in subsets of values at
    the extreme ends of the distribution, we can perform the sort for these
    subsets instead, which should cut down on memory use.

    Returns the upper and lower map value cutoffs (as Python floats).
    """
    map_values = map.as_1d()
    size = map_values.size()
    # upper limit
    i_bin_plus = -1
    for i_bin, value in enumerate(self.v_values()) :
      if ((value*100) <= vol_cutoff_plus_percent) :
        i_bin_plus = i_bin - 1
        break
    assert (i_bin_plus >= 0)
    cutoffp_lower_limit = self.arguments()[i_bin_plus]
    top_values = map_values.select(map_values >= cutoffp_lower_limit)
    i_upper = min(int(size * (vol_cutoff_plus_percent / 100.)),
                  top_values.size())
    s = flex.sort_permutation(top_values)
    top_values_sorted = top_values.select(s)
    del s
    assert (top_values_sorted.size() >= i_upper)
    cutoffp = top_values_sorted[-i_upper]
    del top_values
    del top_values_sorted
    # lower limit
    i_bin_minus = -1
    for i_bin, value in enumerate(self.c_values()) :
      if ((value*100) > vol_cutoff_minus_percent) :
        i_bin_minus = i_bin
        break
    assert (i_bin_minus >= 0)
    cutoffm_upper_limit = self.arguments()[i_bin_minus]
    bottom_values = map_values.select(map_values <= cutoffm_upper_limit)
    i_lower = min(int(size * (vol_cutoff_minus_percent / 100.)),
                  bottom_values.size() - 1)
    s = flex.sort_permutation(bottom_values)
    bottom_values_sorted = bottom_values.select(s)
    del s
    assert (bottom_values_sorted.size() > i_lower)
    cutoffm = bottom_values_sorted[i_lower]
    del bottom_values
    del bottom_values_sorted
    return cutoffp, cutoffm

class peak_list(ext.peak_list):

  def __init__(self, data,
                     tags,
                     peak_search_level=1,
                     max_peaks=0,
                     peak_cutoff=None,
                     interpolate=True):
    if (peak_cutoff is None):
      ext.peak_list.__init__(self,
        data, tags, peak_search_level, max_peaks, interpolate)
    else:
      ext.peak_list.__init__(self,
        data, tags, peak_search_level, peak_cutoff, max_peaks, interpolate)

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling=True):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy(
  to_map=structure_factors_to_map,
  from_map=structure_factors_from_map)

class crystal_gridding(object):

  def __init__(self, unit_cell,
                     d_min=None,
                     resolution_factor=None,
                     step=None,
                     symmetry_flags=None,
                     space_group_info=None,
                     mandatory_factors=None,
                     max_prime=5,
                     assert_shannon_sampling=True,
                     pre_determined_n_real=None):
    if (pre_determined_n_real is None):
      assert [d_min, step].count(None) == 1
      if (step is not None):
        d_min = step*2
        resolution_factor = 0.5
      elif (resolution_factor is None):
        resolution_factor = 1/3
      if (symmetry_flags is not None): assert space_group_info is not None
      if (mandatory_factors is None): mandatory_factors = (1,1,1)
      assert len(mandatory_factors) == 3
    else:
      assert d_min is None
      assert step is None
      assert mandatory_factors is None
    adopt_init_args(self, locals(), hide=True)
    if (pre_determined_n_real is not None):
      self._n_real = pre_determined_n_real
    elif (symmetry_flags is not None):
      self._n_real = determine_gridding(
        unit_cell, d_min, resolution_factor,
        symmetry_flags, space_group_info.type(),
        mandatory_factors, max_prime, assert_shannon_sampling)
    else:
      self._n_real = determine_gridding(
        unit_cell, d_min, resolution_factor,
        mandatory_factors, max_prime, assert_shannon_sampling)

  def _copy_constructor(self, other):
    self._unit_cell = other._unit_cell
    self._d_min = other._d_min
    self._resolution_factor = other._resolution_factor
    self._symmetry_flags = other._symmetry_flags
    self._space_group_info = other._space_group_info
    self._mandatory_factors = other._mandatory_factors
    self._max_prime = other._max_prime
    self._n_real = other._n_real

  def unit_cell(self):
    return self._unit_cell

  def d_min(self):
    return self._d_min

  def resolution_factor(self):
    return self._resolution_factor

  def symmetry_flags(self):
    return self._symmetry_flags

  def space_group_info(self):
    return self._space_group_info

  def change_space_group(self, space_group_info):
    assert (space_group_info.group().refine_gridding(self.n_real())
            == self.n_real())
    self._space_group_info = space_group_info

  def mandatory_factors(self):
    return self._mandatory_factors

  def max_prime(self):
    return self._max_prime

  def n_real(self):
    return self._n_real

  def space_group(self):
    assert self.space_group_info() is not None
    return self.space_group_info().group()

  def crystal_symmetry(self):
    assert self.space_group_info() is not None
    return crystal.symmetry(
      unit_cell=self.unit_cell(),
      space_group_info=self.space_group_info())

  def n_grid_points(self):
    result = 1
    for n in self.n_real():
      result *= n
    return result

  def tags(self):
    return crystal_gridding_tags(self)

class crystal_gridding_tags(crystal_gridding):

  def __init__(self, gridding):
    crystal_gridding._copy_constructor(self, gridding)
    assert gridding.symmetry_flags() is not None
    self._tags = grid_tags(dim=self.n_real())
    self._tags.build(
      space_group_type=self.space_group_info().type(),
      symmetry_flags=self.symmetry_flags())
    assert self._tags.n_grid_misses() == 0

  def tags(self):
    return self._tags

  def peak_search(self, parameters, map, verify_symmetry=True):
    if (parameters is None):
      parameters = peak_search_parameters()
    if (verify_symmetry and libtbx.env.full_testing):
      assert self._tags.verify(map)
    if (map.accessor().is_padded()):
      map = copy(map, flex.grid(map.focus()))
    grid_peaks = peak_list(
      data=map,
      tags=self._tags.tag_array(),
      peak_search_level=parameters.peak_search_level(),
      max_peaks=parameters.max_peaks(),
      peak_cutoff=parameters.peak_cutoff(),
      interpolate=parameters.interpolate())
    if (parameters.min_distance_sym_equiv() is None):
      return grid_peaks
    return peak_cluster_analysis(
      peak_list=grid_peaks,
      special_position_settings=crystal.special_position_settings(
        crystal_symmetry=self.crystal_symmetry(),
        min_distance_sym_equiv=parameters.min_distance_sym_equiv()),
      general_positions_only=parameters.general_positions_only(),
      effective_resolution=parameters.effective_resolution(),
      significant_height_fraction=parameters.significant_height_fraction(),
      cluster_height_fraction=parameters.cluster_height_fraction(),
      min_cross_distance=parameters.min_cross_distance(),
      max_clusters=parameters.max_clusters(),
      min_cubicle_edge=parameters.min_cubicle_edge())

class boxes(object):
  """
  Split box defined by n_real into boxes where each box is a fraction of the
  whole box.
  """
  def __init__(self,
               n_real,
               fraction=None,
               log=None,
               max_boxes=2000,
               prefix=""):
    self.n_real = n_real
    i=0
    n_boxes = 1.e+9
    n_boxes_ = []
    while n_boxes>max_boxes:
      ba,bb,bc = \
        min(10+i,max(3,int(n_real[0]*fraction))), \
        min(10+i,max(3,int(n_real[1]*fraction))), \
        min(10+i,max(3,int(n_real[2]*fraction)))
      n_boxes = self._generate_boxes(ba,bb,bc)
      if(n_boxes_.count(n_boxes)>3): break
      n_boxes_.append(n_boxes)
      i += 1
    assert n_boxes == len(self.starts)
    if(log):
      print >> log, prefix, "n1,n2,n3 (n_real)  :", n_real
      print >> log, prefix, "points per box edge:", ba,bb,bc
      print >> log, prefix, "number of boxes    :", len(self.starts)

  def _generate_boxes(self, ba,bb,bc):
    def regroup(be):
      maxe = be[len(be)-1][1]
      step = int(maxe/len(be))
      result = []
      for i in xrange(len(be)):
        if(i==0):
          l = 0
          r = step
        elif(i==len(be)-1):
          l = i*step
          r = maxe
        else:
          l = i*step
          r = (i+1)*step
        result.append([l,r])
      return result
    be = []
    for i, b in enumerate([ba,bb,bc]):
      be_ = self._box_edges(n_real_1d = self.n_real[i], step=b)
      be_ = regroup(be_)
      be.append(be_)
    self.starts = []
    self.ends = []
    for i in be[0]:
      for j in be[1]:
        for k in be[2]:
          self.starts.append([i[0],j[0],k[0]])
          self.ends.append([i[1],j[1],k[1]])
    return len(self.starts)

  def _box_edges(self, n_real_1d, step):
    limits = []
    for i in range(0, n_real_1d, step): limits.append(i)
    limits.append(n_real_1d)
    box_1d = []
    for i in xrange(len(limits)):
      if(i==0):               box_1d.append([limits[0],  limits[1]])
      elif(i!=len(limits)-1): box_1d.append([limits[i],limits[i+1]])
    return box_1d

class peak_search_parameters(object):

  def __init__(self, peak_search_level=1,
                     max_peaks=0,
                     peak_cutoff=None,
                     interpolate=True,
                     min_distance_sym_equiv=None,
                     general_positions_only=False,
                     effective_resolution=None,
                     significant_height_fraction=None,
                     cluster_height_fraction=None,
                     min_cross_distance=None,
                     max_clusters=None,
                     min_cubicle_edge=5):
    adopt_init_args(self, locals(), hide=True)

  def _copy_constructor(self, other):
    self._peak_search_level = other._peak_search_level
    self._max_peaks = other._max_peaks
    self._peak_cutoff = other._peak_cutoff
    self._interpolate = other._interpolate
    self._min_distance_sym_equiv = other._min_distance_sym_equiv
    self._general_positions_only = other._general_positions_only
    self._effective_resolution = other._effective_resolution
    self._significant_height_fraction = other._significant_height_fraction
    self._cluster_height_fraction = other._cluster_height_fraction
    self._min_cross_distance = other._min_cross_distance
    self._max_clusters = other._max_clusters
    self._min_cubicle_edge = other._min_cubicle_edge

  def peak_search_level(self):
    return self._peak_search_level

  def max_peaks(self):
    return self._max_peaks

  def peak_cutoff(self):
    return self._peak_cutoff

  def interpolate(self):
    return self._interpolate

  def min_distance_sym_equiv(self):
    return self._min_distance_sym_equiv

  def general_positions_only(self):
    return self._general_positions_only

  def effective_resolution(self):
    return self._effective_resolution

  def significant_height_fraction(self):
    return self._significant_height_fraction

  def cluster_height_fraction(self):
    return self._cluster_height_fraction

  def min_cross_distance(self):
    return self._min_cross_distance

  def max_clusters(self):
    return self._max_clusters

  def min_cubicle_edge(self):
    return self._min_cubicle_edge

class cluster_site_info(object):

  def __init__(self, peak_list_index, grid_index, grid_height, site, height):
    self.peak_list_index = peak_list_index
    self.grid_index = grid_index
    self.grid_height = grid_height
    self.site = site
    self.height = height

class peak_cluster_analysis(object):

  def __init__(self, peak_list,
                     special_position_settings,
                     general_positions_only=False,
                     effective_resolution=None,
                     significant_height_fraction=None,
                     cluster_height_fraction=None,
                     min_cross_distance=None,
                     max_clusters=None,
                     min_cubicle_edge=5):
    if (effective_resolution is not None):
      if (significant_height_fraction is None):
          significant_height_fraction = 1/5
      if (cluster_height_fraction is None):
          cluster_height_fraction = 1/3
    if (min_cross_distance is None):
        min_cross_distance = special_position_settings.min_distance_sym_equiv()
    adopt_init_args(self, locals(), hide=True)
    assert self._min_cross_distance is not None
    self._gridding = peak_list.gridding()
    if (effective_resolution is not None):
      self._is_processed = flex.bool(peak_list.size(), False)
    else:
      self._is_processed = None
    if (   effective_resolution is not None
        or debug_peak_cluster_analysis == "use_old"):
      self._site_cluster_analysis = None
    else:
      self._site_cluster_analysis = \
        self._special_position_settings.site_cluster_analysis(
          min_cross_distance=self._min_cross_distance,
          min_self_distance
            =self._special_position_settings.min_distance_sym_equiv(),
          general_positions_only=self._general_positions_only,
          min_cubicle_edge=self._min_cubicle_edge)
    self._peak_list_indices = flex.size_t()
    self._peak_list_index = 0
    self._sites = flex.vec3_double()
    self._heights = flex.double()
    self._fixed_site_indices = flex.size_t()

  def next(self):
    if (self._effective_resolution is not None):
      return self.next_with_effective_resolution()
    else:
      return self.next_site_cluster_analysis()

  def all(self,max_clusters=None):
    if (self._effective_resolution is not None):
      return self.all_with_effective_resolution(max_clusters=max_clusters)
    else:
      return self.all_site_cluster_analysis(max_clusters=max_clusters)

  def __iter__(self):
    while 1:
      site_info = self.next()
      if site_info is None: break
      yield site_info

  def peak_list(self):
    return self._peak_list

  def special_position_settings(self):
    return self._special_position_settings

  def general_positions_only(self):
    return self._general_positions_only

  def effective_resolution(self):
    return self._effective_resolution

  def significant_height_fraction(self):
    return self._significant_height_fraction

  def cluster_height_fraction(self):
    return self._cluster_height_fraction

  def min_cross_distance(self):
    return self._min_cross_distance

  def max_clusters(self):
    return self._max_clusters

  def site_cluster_analysis(self):
    return self._site_cluster_analysis

  def peak_list_indices(self):
    return self._peak_list_indices

  def fixed_site_indices(self):
    return self._fixed_site_indices

  def sites(self):
    return self._sites

  def heights(self):
    return self._heights

  def max_grid_height(self):
    if (self._peak_list.size() == 0):
      return None
    return self._peak_list.heights()[0]

  def append_fixed_site(self, site, height=0):
    if (self._site_cluster_analysis is not None):
      self._site_cluster_analysis.insert_fixed_site_frac(original_site=site)
    self._fixed_site_indices.append(self._sites.size())
    self._sites.append(site)
    self._heights.append(height)
    self._peak_list_indices.append(self._peak_list.size())

  def discard_last(self):
    assert self._peak_list_indices.size() > 0
    if (self._site_cluster_analysis is not None):
      self._site_cluster_analysis.discard_last()
    self._peak_list_indices.pop_back()
    self._sites.pop_back()
    self._heights.pop_back()

  def next_site_cluster_analysis(self):
    while 1:
      peak_list_index = self._peak_list_index
      if (peak_list_index >= self._peak_list.size()): return None
      self._peak_list_index += 1
      site_symmetry = self._special_position_settings.site_symmetry(
        site=self._peak_list.sites()[peak_list_index])
      site = site_symmetry.exact_site()
      if (not self._site_cluster_analysis.process_site_frac(
                original_site=site,
                site_symmetry_ops=site_symmetry)): continue
      height = self._peak_list.heights()[peak_list_index]
      self._peak_list_indices.append(peak_list_index)
      self._sites.append(site)
      self._heights.append(height)
      return cluster_site_info(
        peak_list_index=peak_list_index,
        grid_index=self._peak_list.grid_indices(peak_list_index),
        grid_height=self._peak_list.grid_heights()[peak_list_index],
        site=site,
        height=height)

  def all_site_cluster_analysis(self, max_clusters=None):
    if (max_clusters is None):
      max_clusters = self._max_clusters
    assert max_clusters is not None
    while 1:
      if (self._sites.size() >= max_clusters): break
      if (self.next_site_cluster_analysis() is None): break
    return self

  def next_with_effective_resolution(self):
    while 1:
      peak_list_index = self._peak_list_index
      if (peak_list_index >= self._peak_list.size()): return None
      self._peak_list_index += 1
      if (self._is_processed is not None):
        if (self._is_processed[peak_list_index]): continue
        self._is_processed[peak_list_index] = True
      grid_index = self._peak_list.grid_indices(peak_list_index)
      grid_height = self._peak_list.grid_heights()[peak_list_index]
      site = self._peak_list.sites()[peak_list_index]
      height = self._peak_list.heights()[peak_list_index]
      site_symmetry = self._special_position_settings.site_symmetry(site)
      if (    self._general_positions_only
          and not site_symmetry.is_point_group_1()):
        continue
      site = site_symmetry.exact_site()
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      keep = True
      if (self._sites.size() > 250):
        import warnings
        warnings.warn(
          message="This function should not be used for"
                  " processing a large number of peaks.",
          category=RuntimeWarning)
      for s in self._sites:
        dist = sgtbx.min_sym_equiv_distance_info(equiv_sites, s).dist()
        if (dist < self._min_cross_distance):
          keep = False
          break
      if (keep == True):
        if (    self._effective_resolution is not None
            and (   self._heights.size() == 0
                 or height <   self._heights[0]
                             * self._significant_height_fraction)):
            site, height = self._accumulate_significant(
              site, height, site_symmetry, equiv_sites)
        self._peak_list_indices.append(peak_list_index)
        self._sites.append(site)
        self._heights.append(height)
        return cluster_site_info(
          peak_list_index=peak_list_index,
          grid_index=grid_index,
          grid_height=grid_height,
          site=site,
          height=height)

  def _accumulate_significant(self, site, height, site_symmetry, equiv_sites):
    unit_cell = self.special_position_settings().unit_cell()
    orth = unit_cell.orthogonalize
    frac = unit_cell.fractionalize
    sum_w_sites = matrix.col(orth(site)) * height
    sum_w = height
    height_cutoff = height * self._cluster_height_fraction
    for i in xrange(self._peak_list_index, self._peak_list.size()):
      if (self._is_processed[i]): continue
      other_height = self._peak_list.heights()[i]
      if (other_height < height_cutoff): break
      other_site = self._peak_list.sites()[i]
      other_site_symmetry = self._special_position_settings.site_symmetry(
        other_site)
      if (    self._general_positions_only
          and not other_site_symmetry.is_point_group_1()):
        self._is_processed[i] = True
        continue
      other_site = other_site_symmetry.exact_site()
      dist_info = sgtbx.min_sym_equiv_distance_info(equiv_sites, other_site)
      dist = dist_info.dist()
      if (dist < self._min_cross_distance):
        self._is_processed[i] = True
        close_site = dist_info.apply(flex.vec3_double([other_site]))[0]
        close_site = site_symmetry.special_op() * close_site
        sum_w_sites += matrix.col(orth(close_site)) * other_height
        sum_w += other_height
    return frac(sum_w_sites / sum_w), height

  def all_with_effective_resolution(self, max_clusters=None):
    if (max_clusters is None):
      max_clusters = self._max_clusters
    assert max_clusters is not None
    while 1:
      if (self._sites.size() >= max_clusters): break
      if (self.next_with_effective_resolution() is None): break
    return self

def region_density_correlation(
      large_unit_cell,
      large_d_min,
      large_density_map,
      sites_cart,
      site_radii,
      work_scatterers):
  sites_frac_large = large_unit_cell.fractionalize(sites_cart)
  large_frac_min = sites_frac_large.min()
  large_frac_max = sites_frac_large.max()
  large_n_real = large_density_map.focus()
  from scitbx import fftpack
  from libtbx.math_utils import ifloor, iceil
  large_ucp = large_unit_cell.parameters()
  small_n_real = [0,0,0]
  small_origin_in_large_grid = [0,0,0]
  small_abc = [0,0,0]
  sites_frac_shift = [0,0,0]
  for i in xrange(3):
    grid_step = large_ucp[i] / large_n_real[i]
    buffer = large_d_min / grid_step
    grid_min = ifloor(large_frac_min[i] * large_n_real[i] - buffer)
    grid_max = iceil(large_frac_max[i] * large_n_real[i] + buffer)
    min_grid = grid_max - grid_min + 1
    small_n_real[i] = fftpack.adjust_gridding(min_grid=min_grid, max_prime=5)
    if (small_n_real[i] < large_n_real[i]):
      shift_min = (small_n_real[i] - min_grid) // 2
      small_origin_in_large_grid[i] = grid_min - shift_min
      small_abc[i] = small_n_real[i] * grid_step
      sites_frac_shift[i] = small_origin_in_large_grid[i] / large_n_real[i]
    else:
      small_n_real[i] = large_n_real[i]
      small_origin_in_large_grid[i] = 0
      small_abc[i] = large_ucp[i]
      sites_frac_shift[i] = 0
  sites_cart_shift = large_unit_cell.orthogonalize(sites_frac_shift)
  sites_cart_small = sites_cart - sites_cart_shift
  from cctbx import xray
  small_xray_structure = xray.structure(
    crystal_symmetry=crystal.symmetry(
      unit_cell=tuple(small_abc)+large_ucp[3:],
      space_group_symbol="P1"),
    scatterers=work_scatterers)
  small_xray_structure.set_sites_cart(sites_cart=sites_cart_small)
  small_f_calc = small_xray_structure.structure_factors(
    d_min=large_d_min).f_calc()
  small_gridding = crystal_gridding(
    unit_cell=small_f_calc.unit_cell(),
    space_group_info=small_f_calc.space_group_info(),
    pre_determined_n_real=small_n_real)
  from cctbx import miller
  small_fft_map = miller.fft_map(
    crystal_gridding=small_gridding,
    fourier_coefficients=small_f_calc)
  small_fft_map.apply_sigma_scaling()
  small_map = small_fft_map.real_map_unpadded()
  grid_indices = grid_indices_around_sites(
    unit_cell=small_xray_structure.unit_cell(),
    fft_n_real=small_n_real,
    fft_m_real=small_n_real,
    sites_cart=sites_cart_small,
    site_radii=site_radii)
  small_copy_from_large_map = copy(
    map_unit_cell=large_density_map,
    first=small_origin_in_large_grid,
    last=matrix.col(small_origin_in_large_grid)
       + matrix.col(small_n_real)
       - matrix.col((1,1,1)))
  assert small_copy_from_large_map.all() == small_map.all()
  corr = flex.linear_correlation(
    x=small_map.select(grid_indices),
    y=small_copy_from_large_map.select(grid_indices))
  if (not corr.is_well_defined()):
    return None
  return corr.coefficient()

def ccv(map_1, map_2, modified, centered, cutoff=None, n_bins=10000):
  if(modified):
    map_1 = volume_scale(map=map_1, n_bins=n_bins).map_data()
    map_2 = volume_scale(map=map_2, n_bins=n_bins).map_data()
  if(cutoff is not None):
    map_1 = map_1 - cutoff
    map_2 = map_2 - cutoff
    s1 = map_1 < 0
    s2 = map_2 < 0
    map_1 = map_1.set_selected(s1, 0)
    map_2 = map_2.set_selected(s2, 0)
    def corr(x, y, centered):
      s1 = x > 0
      s2 = y > 0
      s = s1 | s2
      s = s.iselection()
      x_ = x.select(s)
      y_ = y.select(s)
      return flex.linear_correlation(x = x_, y = y_,
        subtract_mean = centered).coefficient()
    return corr(x = map_1, y = map_2, centered = centered)
  else:
    return flex.linear_correlation(x = map_1.as_1d(), y = map_2.as_1d(),
      subtract_mean = centered).coefficient()

class spherical_variance_around_point (object) :
  def __init__ (self,
      real_map,
      unit_cell,
      site_cart,
      radius,
      n_points=40,
      spline_interpolation=True,
      write_sphere_points_to_pdb_file=None) :
    self.site_cart = site_cart
    self.radius = radius
    assert n_points>0
    sphere_points = []
    x, y, z = site_cart
    # reference: "Distributing many points on a sphere" by E.B. Saff and
    #     A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
    # derived from http://packinon.sourceforge.net/py_progs/pg_saff.html
    for k in range(1,n_points+1):
      h = -1 + 2 * (k - 1) / float(n_points - 1)
      theta = math.acos(h)
      if (k == 1) or (k == n_points) :
        phi = 0
      else:
        phi = (old_phi + 3.6/math.sqrt(n_points*(1-h*h))) % (2*math.pi)
      sphere_points.append((
        x + math.sin(phi)*math.sin(theta),
        y + math.cos(theta),
        z + math.cos(phi)*math.sin(theta) ))
      old_phi = phi
    map_values = flex.double()
    for point in sphere_points :
      site_frac = unit_cell.fractionalize(site_cart=point)
      value_at_point = real_map.tricubic_interpolation(site_frac)
      map_values.append(value_at_point)
    self.min = flex.min(map_values)
    self.max = flex.max(map_values)
    self.mean = flex.mean(map_values)
    self.standard_deviation = map_values.standard_deviation_of_the_sample()
    if (write_sphere_points_to_pdb_file is not None) :
      f = open(write_sphere_points_to_pdb_file, "w")
      for i, point in enumerate(sphere_points) :
        f.write(
          "HETATM    1  O   HOH A   1     %7.3f %7.3f %7.3f  1.00 20.00\n"%
          point)
      f.close()

  def show (self, out=None, prefix="") :
    if (out is None) : out = sys.stdout
    print >> out, "%sMap values around point [%g, %g, %g], radius=%g:" % \
      (prefix, self.site_cart[0], self.site_cart[1], self.site_cart[2],
       self.radius)
    print >> out, "%s  min=%.2f  max=%.2f  mean=%.2f  stddev=%.2f" % \
      (prefix, self.min, self.max, self.mean, self.standard_deviation)

def principal_axes_of_inertia (
    real_map,
    site_cart,
    unit_cell,
    radius) :
  assert (radius > 0)
  import scitbx.math
  from libtbx.utils import n_dim_index_from_one_dim
  sel = grid_indices_around_sites(
    unit_cell=unit_cell,
    fft_n_real=real_map.focus(),
    fft_m_real=real_map.all(),
    sites_cart=flex.vec3_double([site_cart]),
    site_radii=flex.double([radius]))
  grid_to_cart = grid2cart(real_map.focus(),
    unit_cell.orthogonalization_matrix())
  sites = flex.vec3_double()
  values = flex.double()
  site_frac = unit_cell.fractionalize(site_cart=site_cart)
  weight_scale = real_map.tricubic_interpolation(site_frac)
  if (weight_scale <= 0) :
    weight_scale = flex.max(real_map.as_1d())
  for i_seq in sel :
    uvw = n_dim_index_from_one_dim(i_seq, real_map.all())
    xyz = grid_to_cart(uvw)
    sites.append(xyz)
    site_frac = unit_cell.fractionalize(site_cart=xyz)
    # XXX weights of zero risk crashing the underlying routine
    values.append(max(0.000001, real_map[i_seq]) / weight_scale)
  assert (flex.max(values) >= 0)
  return scitbx.math.principal_axes_of_inertia(
    points=sites,
    weights=values)

class local_scale(object):
  def __init__(
        self,
        crystal_gridding,
        crystal_symmetry,
        f_map=None,
        map_data=None,
        miller_array=None,
        d_min = None): #XXX =1: more features and noise
    # process inputs
    assert [f_map, map_data].count(None) == 1
    if(f_map is not None):
      import cctbx.miller
      fft_map = cctbx.miller.fft_map(
        crystal_gridding     = crystal_gridding,
        fourier_coefficients = f_map)
      fft_map.apply_sigma_scaling()
      map_data = fft_map.real_map_unpadded()
    #
    self.map_result = None
    self.map_coefficients = None
    b = boxes(
      n_real   = crystal_gridding.n_real(),
      fraction = 0.03)
    # Loop over boxes, fill map_result with one box at a time
    self.map_result = flex.double(flex.grid(b.n_real))
    for s,e in zip(b.starts, b.ends):
      box = copy(map_data, s, e)
      box.reshape(flex.grid(box.all()))
      mi,ma,me = box.as_1d().min_max_mean().as_tuple()
      if(mi < ma):
        box = volume_scale(map = box, n_bins = 1000).map_data()
      set_box(
        map_data_from = box,
        map_data_to   = self.map_result,
        start         = s,
        end           = e)
    sd = self.map_result.sample_standard_deviation()
    self.map_result = self.map_result/sd
    if(miller_array is not None):
      complete_set = miller_array
      if(d_min is not None):
        d_min = miller_array.d_min()
        complete_set = miller_array.complete_set(d_min=d_min)
      self.map_coefficients = complete_set.structure_factors_from_map(
        map            = self.map_result,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)

def map_peak_3d_as_2d(
      map_data,
      unit_cell,
      center_cart,
      radius,
      step = 0.01,
      s_angle_sampling_step = 10,
      t_angle_sampling_step = 10):
  rho_1d = flex.double()
  dist = flex.double()
  radius = int(radius*100)+1
  step = int(step*100)
  for r in range(0,radius,step):
    r = r/100.
    dist.append(r)
    rho = flex.double()
    for s in xrange(0,360,s_angle_sampling_step):
      for t in xrange(0,360,t_angle_sampling_step):
        xc,yc,zc = scitbx.math.point_on_sphere(r=r, s_deg=s, t_deg=t,
          center=center_cart)
        xf,yf,zf = unit_cell.fractionalize([xc,yc,zc])
        rho.append(map_data.eight_point_interpolation([xf,yf,zf]))
    rho_1d.append(flex.mean(rho))
  return dist, rho_1d

class positivity_constrained_density_modification(object):
  def __init__(self, f, f_000, n_cycles=100, resolution_factor=0.25, d_min=None,
               crystal_gridding=None, complete_set=None):
    self.f = f
    self.d_min = d_min
    self.map = None
    self.crystal_gridding = crystal_gridding
    from cctbx import miller
    if(self.d_min is None): self.d_min = self.f.d_min()
    if(complete_set is None):
      complete_set = self.f.complete_set(d_min = self.d_min)
    if(self.crystal_gridding is None):
      self.crystal_gridding = self.f.crystal_gridding(
        d_min                   = d_min,
        resolution_factor       = resolution_factor,
        grid_step               = None,
        symmetry_flags          = None,
        mandatory_factors       = None,
        max_prime               = 5,
        assert_shannon_sampling = True)
    self.f_mod = self.f.deep_copy()
    for i in xrange(n_cycles):
      fft_map = miller.fft_map(
        crystal_gridding     = self.crystal_gridding,
        fourier_coefficients = self.f_mod,
        f_000                = f_000)
      if(f_000 is not None): fft_map.apply_volume_scaling()
      self.map = fft_map.real_map_unpadded()
      convert_to_non_negative(self.map, 0)
      self.f_mod = complete_set.structure_factors_from_map(
        map            = self.map,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)
      self.f_mod = self.f.complete_with(other = self.f_mod, scale=True,
        replace_phases=True)
      #self.assert_equal()

  def assert_equal(self):
    from libtbx.test_utils import approx_equal
    x, y = self.f, self.f_mod
    x,y = x.common_sets(y)
    x = abs(x).data()
    y = abs(y).data()
    assert approx_equal(x, y)

def d_min_from_map(map_data, unit_cell, resolution_factor=1./2.):
  a,b,c = unit_cell.parameters()[:3]
  nx,ny,nz = map_data.all()
  d1,d2,d3 = \
    a/nx/resolution_factor,\
    b/ny/resolution_factor,\
    c/nz/resolution_factor
  return max(d1,d2,d3)

def resolution_from_map_and_model(map_data, xray_structure):
  """
  Given map and model estimate resolution by maximizing map CC(map, model-map).
  """
  from cctbx import miller
  xrs = xray_structure
  sel = grid_indices_around_sites(
    unit_cell  = xrs.unit_cell(),
    fft_n_real = map_data.focus(),
    fft_m_real = map_data.all(),
    sites_cart = xrs.sites_cart(),
    site_radii = flex.double(xrs.scatterers().size(),5))
  map_data_selected_as_1d = map_data.select(sel).as_1d()
  cg = crystal_gridding(
    unit_cell             = xrs.unit_cell(),
    space_group_info      = xrs.space_group_info(),
    pre_determined_n_real = map_data.accessor().all())
  def optimize(d_min_min, d_min_max, step):
    d_min = d_min_min
    d_min_result = None
    cc_best=-1.e6
    while d_min < d_min_max+1.e-6:
      f_calc = xrs.structure_factors(d_min=d_min).f_calc()
      fft_map = miller.fft_map(
        crystal_gridding     = cg,
        fourier_coefficients = f_calc)
      map_data_ = fft_map.real_map_unpadded()
      cc = flex.linear_correlation(
        x=map_data_selected_as_1d,
        y=map_data_.select(sel).as_1d()).coefficient()
      if(cc > cc_best):
        cc_best = cc
        d_min_result = d_min
      d_min += step
    return d_min_result
  # coarse estimate
  d_min_max = d_min_from_map(
    map_data=map_data, unit_cell=xrs.unit_cell(), resolution_factor=0.2)
  d_min_min = d_min_from_map(
    map_data=map_data, unit_cell=xrs.unit_cell(), resolution_factor=0.5)
  step = (d_min_max-d_min_min)/5
  d_min_result = optimize(d_min_min=d_min_min, d_min_max=d_min_max, step=step)
  # fine estimate
  d_min_min=d_min_result-step
  d_min_max=d_min_result+step
  step = (d_min_max-d_min_min)/10
  d_min_result = optimize(d_min_min=d_min_min, d_min_max=d_min_max, step=step)
  return d_min_result

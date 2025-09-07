"""Tools for map analysis and manipulation"""
from __future__ import absolute_import, division, print_function
import cctbx.sgtbx

import boost_adaptbx.boost.python as bp
from six.moves import range
from six.moves import zip
ext = bp.import_ext("cctbx_maptbx_ext")
from cctbx_maptbx_ext import *

from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx import matrix
from scitbx.python_utils import dicts
from libtbx import adopt_init_args
from libtbx.utils import Sorry
import libtbx.load_env
import math
import sys, os
import scitbx.math
from cctbx import adptbx
from libtbx import group_args
from scitbx import fftpack
from libtbx.test_utils import approx_equal
from cctbx import uctbx
import scitbx.math

debug_peak_cluster_analysis = os.environ.get(
  "CCTBX_MAPTBX_DEBUG_PEAK_CLUSTER_ANALYSIS", "")

@bp.inject_into(connectivity)
class _():

  def get_blobs_boundaries_tuples(self):
    """
    get lists of minimum and maximum coordinates for each connected
    region.
    returns 2 lists of tuples: first is minimum, second is maximum coordinates.
    [(x0, y0, z0), (x1, y1, z1), ...] where 0, 1, ... - number of region
    """
    boundaries = self.get_blobs_boundaries()
    regs = self.regions()
    min_boundaries = []
    max_boundaries = []
    for i in range(len(regs)):
      minb = (boundaries[0, i, 0], boundaries[0, i, 1], boundaries[0, i, 2])
      maxb = (boundaries[1, i, 0], boundaries[1, i, 1], boundaries[1, i, 2])
      min_boundaries.append(minb)
      max_boundaries.append(maxb)
    return min_boundaries, max_boundaries

def smooth_map(map, crystal_symmetry, rad_smooth, method = "exp",
     non_negative = True):
  """Smooth a map with radius rad_smooth, using one of these
   methods:  "exp" (Gaussian smoothing), "box_average" (local box average),
   "top_hat" (Top-hat smoothing, spherical average done in reciprocal
   space)
  """

  from cctbx import miller
  assert method in ["exp", "box_average", "top_hat"]
  map_smooth = None
  if (method == "exp"):
    f_map = miller.structure_factor_box_from_map(
      map              = map,
      crystal_symmetry = crystal_symmetry,
      include_000      = True)
    ddd = f_map.d_spacings().data()
    ddd.set_selected(ddd  ==  -1 , 1.e+10)  # d for (0, 0, 0) was set to -1
    ss = 1./flex.pow2(ddd) / 4.
    b_smooth = 8*math.pi**2*rad_smooth**2
    smooth_scale = flex.exp(-b_smooth*ss)
    f_map = f_map.array(data = f_map.data()*smooth_scale)
    from cctbx.maptbx import crystal_gridding
    cg = crystal_gridding(
      unit_cell             = crystal_symmetry.unit_cell(),
      space_group_info      = crystal_symmetry.space_group_info(),
      pre_determined_n_real = map.all())
    fft_map = miller.fft_map(
      crystal_gridding     = cg,
      fourier_coefficients = f_map)
    fft_map.apply_volume_scaling()
    map_smooth = fft_map.real_map_unpadded()
    if non_negative:
      map_smooth = map_smooth.set_selected(map_smooth<0., 0)

  elif (method  == "top_hat"):

    # use same grid as original
    from cctbx.maptbx import crystal_gridding
    cg = crystal_gridding(
        unit_cell = crystal_symmetry.unit_cell(),
        space_group_info = crystal_symmetry.space_group_info(),
        pre_determined_n_real = map.all())

    average_value = map.as_1d().min_max_mean().mean

    # Get structure factors with f_000
    f_map = miller.structure_factor_box_from_map(
      map              = map,
      crystal_symmetry = crystal_symmetry,
      include_000      = True)

    # The d_spacings get a value of -1 for the f000 term. Set it to a big number
    ddd = f_map.d_spacings().data()
    ddd.set_selected(ddd < 0, 1.e+10)
    d_min = f_map.d_min()

    # G-function for top hat (FT of top hat)
    #complete_set = f_map.complete_set(include_f000 = True)
    sphere_reciprocal = f_map.g_function(R=rad_smooth) # top hat function

    # FT (map) * FT(top hat) = FT (convolution of map and top hat)
    fourier_coeff = f_map.array(data=f_map.data()*sphere_reciprocal)

    # Convolution of map and top hat
    fft_map = fourier_coeff.fft_map(d_min = d_min, crystal_gridding = cg)
    fft_map.apply_volume_scaling()
    map_smooth = fft_map.real_map_unpadded()


    if non_negative:
      map_smooth = map_smooth.set_selected(map_smooth<0., 0)

  elif(method == "box_average"): # assume 0/1 binary map
    assert abs(flex.max(map)-1.)<1.e-6
    mmin = flex.min(map)
    assert mmin<1.e-6 and mmin>= 0.0
    map_smooth = map.deep_copy()
    for i in range(3):
      maptbx.map_box_average(
        map_data      = map_smooth,
        index_span    = 1)
    for i in range(3):
      maptbx.map_box_average(
        map_data      = map_smooth,
        cutoff        = 0.99,
        index_span    = 1)
  return map_smooth

class d99(object):
  def __init__(self, map = None, f_map = None, crystal_symmetry = None):
    adopt_init_args(self, locals())
    if(map is not None):
      assert f_map is None
      assert crystal_symmetry is not None
      map = shift_origin_if_needed(map_data = map).map_data
      from cctbx import miller
      self.f_map = miller.structure_factor_box_from_map(
        map = map, crystal_symmetry = crystal_symmetry)
    else:
      assert [map, crystal_symmetry].count(None) == 2
    self.d_spacings = self.f_map.d_spacings().data()
    self.d_max, self.d_min = flex.max(self.d_spacings), flex.min(self.d_spacings)
    o = ext.d99(
      f          = self.f_map.data(),
      d_spacings = self.d_spacings,
      hkl        = self.f_map.indices(),
      cutoff     = 0.99)
    self.result = group_args(
      d99 = o.d_min())

  def show(self, log):
    fmt = "%12.6f %8.6f"
    for d_min, cc in zip(self.result.d_mins, self.result.ccs):
      print(fmt%(d_min, cc), file = log)

def assert_same_gridding(map_1, map_2,
                         Sorry_message = "Maps have different gridding."):
  f1 = map_1.focus() == map_2.focus()
  f2 = map_1.origin() == map_2.origin()
  f3 = map_1.all() == map_2.all()
  if([f1, f2, f3].count(True)!= 3):
    raise Sorry(Sorry_message)

def shift_origin_if_needed(map_data = None,
    sites_cart = None,
    crystal_symmetry = None,
    ncs_object = None,
    origin_grid_units = None,
    n_xyz = None,
    ):

  if not map_data:
    assert origin_grid_units and n_xyz
    shift_needed = True

  else: # usual
    shift_needed = not \
    (map_data.focus_size_1d() > 0 and map_data.nd()  ==  3 and
     map_data.is_0_based())

  shift_frac = None
  shift_cart = None
  if(shift_needed):
    if map_data:
      N = map_data.all()
      O = map_data.origin()
      map_data = map_data.shift_origin()
    else:
      N = n_xyz
      O = origin_grid_units

    original_origin_grid_units = O
    original_origin_cart = (0, 0, 0)
    if crystal_symmetry:
      if(not crystal_symmetry.space_group().type().number() in [0, 1]):
        raise Sorry("Space groups other than P1 are not supported.")
      a, b, c = crystal_symmetry.unit_cell().parameters()[:3]
      fm = crystal_symmetry.unit_cell().fractionalization_matrix()
      sx, sy, sz = O[0]/N[0], O[1]/N[1], O[2]/N[2]
      shift_frac = [-sx, -sy, -sz]
      shift_cart = crystal_symmetry.unit_cell().orthogonalize(shift_frac)
      original_origin_cart = tuple(-matrix.col(shift_cart))
      if(sites_cart is not None):
        sites_cart = sites_cart + flex.vec3_double(sites_cart.size(), shift_cart)
      if ncs_object:
        ncs_object = ncs_object.deep_copy(coordinate_offset = shift_cart)
    else:
      original_origin_grid_units = None
      original_origin_cart = None
  else:
    original_origin_grid_units = (0, 0, 0)
    original_origin_cart = (0, 0, 0)
  return group_args(
    map_data   = map_data,
    ncs_object = ncs_object,
    sites_cart = sites_cart,
    shift_frac = shift_frac,
    shift_cart = shift_cart,
    original_origin_grid_units = original_origin_grid_units,
    original_origin_cart = original_origin_cart)

def value_at_closest_grid_point(map, x_frac):
  return map[closest_grid_point(map.accessor(), x_frac)]

flex.int.value_at_closest_grid_point = value_at_closest_grid_point
flex.double.value_at_closest_grid_point = value_at_closest_grid_point
flex.double.eight_point_interpolation = eight_point_interpolation
flex.double.eight_point_interpolation_with_gradients = \
  eight_point_interpolation_with_gradients
flex.double.quadratic_interpolation_with_gradients = \
  quadratic_interpolation_with_gradients
flex.double.tricubic_interpolation = tricubic_interpolation
flex.double.tricubic_interpolation_with_gradients = tricubic_interpolation_with_gradients

def cc_peak(cutoff, map_1 = None, map_2 = None, map_coeffs_1 = None, map_coeffs_2 = None):
  """
  Compute CCpeak as described in
    Acta Cryst. (2014). D70, 2593-2606
    Metrics for comparison of crystallographic maps
    A. Urzhumtsev, P. V. Afonine, V. Y. Lunin, T. C. Terwilliger and P. D. Adams
  """
  from cctbx import miller
  assert [map_1, map_2].count(None) in [0, 2]
  assert [map_coeffs_1, map_coeffs_2].count(None) in [0, 2]
  if([map_1, map_2].count(None) == 0):
    # Maps are assumed to be quantile rank scaled (HE).
    return ext.cc_peak(map_1 = map_1, map_2 = map_2, cutoff = cutoff)
  elif([map_coeffs_1, map_coeffs_2].count(None) == 0):
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
    return ext.cc_peak(map_1 = m1_he, map_2 = m2_he, cutoff = cutoff)
  else:
    raise Sorry("Combination of inputs not supported.")

def map_accumulator(n_real, use_max_map, smearing_b = 5, max_peak_scale = 2,
                    smearing_span = 10, use_exp_table = True):
  """
  Good defaults for 2mFo-DFc type maps:
    smearing_b = 1, max_peak_scale = 100, smearing_span = 5
  """
  return ext.map_accumulator(n_real = n_real, smearing_b = smearing_b,
    max_peak_scale = max_peak_scale, smearing_span = smearing_span,
      use_exp_table = use_exp_table, use_max_map = use_max_map)

def peak_volume_estimate(map_data, sites_cart, crystal_symmetry, cutoff,
      atom_radius = 1.5):
  v = flex.double()
  sites_frac = crystal_symmetry.unit_cell().fractionalize(sites_cart)
  for sc, sf in zip(sites_cart, sites_frac):
    if(map_data.value_at_closest_grid_point(sf)>= cutoff):
      sel = grid_indices_around_sites(
        unit_cell  = crystal_symmetry.unit_cell(),
        fft_n_real = map_data.focus(),
        fft_m_real = map_data.all(),
        sites_cart = flex.vec3_double([sc]),
        site_radii = flex.double([atom_radius]*1))
      v.append((map_data.select(sel)>= cutoff).count(True))
  r = flex.min_default(v, None)
  if(r == 0): return None
  return r

def truncate(map_data, by_sigma_less_than, scale_by, set_value = 0):
  """
  Trunate map inplace by standard deviation (sigma) while scale it with
  specified scale, such as volume (scale_by = 1/volume) or sigma
  (scale_by = 1/standard_deviation). Input map_data is expected to be unscaled (
  right out of FT).
  """
  sigma = statistics(map_data).sigma()
  if(sigma  ==  0):
    map_data = map_data*scale_by
    return
  ext.truncate(
    map_data           = map_data,
    standard_deviation = sigma,
    by_sigma_less_than = by_sigma_less_than,
    scale_by           = scale_by,
    set_value          = set_value)

def mask(xray_structure,
         n_real,
         mask_value_inside_molecule = 0,
         mask_value_outside_molecule = 1,
         solvent_radius = 0,
         atom_radius = None):
  xrs_p1 = xray_structure.expand_to_p1(sites_mod_positive = True)
  if(atom_radius is None):
    from cctbx.masks import vdw_radii_from_xray_structure
    atom_radii = vdw_radii_from_xray_structure(xray_structure = xrs_p1)
  else:
    atom_radii = flex.double(xrs_p1.scatterers().size(), atom_radius)
  return ext.mask(
    sites_frac                  = xrs_p1.sites_frac(),
    unit_cell                   = xrs_p1.unit_cell(),
    n_real                      = n_real,
    mask_value_inside_molecule  = mask_value_inside_molecule,
    mask_value_outside_molecule = mask_value_outside_molecule,
    radii                       = atom_radii + solvent_radius)

class statistics(ext.statistics):

  def __init__(self, map):
    ext.statistics.__init__(self, map)

@bp.inject_into(ext.statistics)
class _():

  def show_summary(self, f = None, prefix = ""):
    if (f is None): f = sys.stdout
    print(prefix + "max %.6g" % (self.max()), file = f)
    print(prefix + "min %.6g" % (self.min()), file = f)
    print(prefix + "mean %.6g" % (self.mean()), file = f)
    print(prefix + "sigma %.6g" % (self.sigma()), file = f)

use_space_group_symmetry = sgtbx.search_symmetry_flags(
  use_space_group_symmetry = True)

@bp.inject_into(ext.histogram)
class _():

  """
  Injector for extending cctbx.maptbx.histogram
  """
  # XXX make a method of scitbx
  def get_percentile_cutoffs(self, map, vol_cutoff_plus_percent,
      vol_cutoff_minus_percent):
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
    for i_bin, value in enumerate(self.v_values()):
      if ((value*100) <=  vol_cutoff_plus_percent):
        i_bin_plus = i_bin - 1
        break
    assert (i_bin_plus >=  0)
    cutoffp_lower_limit = self.arguments()[i_bin_plus]
    top_values = map_values.select(map_values >=  cutoffp_lower_limit)
    i_upper = min(int(size * (vol_cutoff_plus_percent / 100.)),
                  top_values.size())
    s = flex.sort_permutation(top_values)
    top_values_sorted = top_values.select(s)
    del s
    assert (top_values_sorted.size() >=  i_upper)
    cutoffp = top_values_sorted[-i_upper]
    del top_values
    del top_values_sorted
    # lower limit
    i_bin_minus = -1
    for i_bin, value in enumerate(self.c_values()):
      if ((value*100) > vol_cutoff_minus_percent):
        i_bin_minus = i_bin
        break
    assert (i_bin_minus >=  0)
    cutoffm_upper_limit = self.arguments()[i_bin_minus]
    bottom_values = map_values.select(map_values <=  cutoffm_upper_limit)
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
                     peak_search_level = 1,
                     max_peaks = 0,
                     peak_cutoff = None,
                     interpolate = True):
    if (peak_cutoff is None):
      ext.peak_list.__init__(self,
        data, tags, peak_search_level, max_peaks, interpolate)
    else:
      ext.peak_list.__init__(self,
        data, tags, peak_search_level, peak_cutoff, max_peaks, interpolate)

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling = True):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy(
  to_map = structure_factors_to_map,
  from_map = structure_factors_from_map)

class crystal_gridding(object):

  def __init__(self, unit_cell,
                     d_min = None,
                     resolution_factor = None,
                     step = None,
                     symmetry_flags = None,
                     space_group_info = None,
                     mandatory_factors = None,
                     max_prime = 5,
                     assert_shannon_sampling = True,
                     pre_determined_n_real = None):
    if (pre_determined_n_real is None):
      assert [d_min, step].count(None)  ==  1
      if (step is not None):
        d_min = step*2
        resolution_factor = 0.5
      elif (resolution_factor is None):
        resolution_factor = 1/3
      if (symmetry_flags is not None): assert space_group_info is not None
      if (mandatory_factors is None): mandatory_factors = (1, 1, 1)
      assert len(mandatory_factors)  ==  3
    else:
      assert d_min is None
      assert step is None
      assert mandatory_factors is None
    adopt_init_args(self, locals(), hide = True)
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
             ==  self.n_real())
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
      unit_cell = self.unit_cell(),
      space_group_info = self.space_group_info())

  def n_grid_points(self):
    result = 1
    for n in self.n_real():
      result *=  n
    return result

  def tags(self):
    return crystal_gridding_tags(self)

class crystal_gridding_tags(crystal_gridding):

  def __init__(self, gridding):
    crystal_gridding._copy_constructor(self, gridding)
    assert gridding.symmetry_flags() is not None
    self._tags = grid_tags(dim = self.n_real())
    self._tags.build(
      space_group_type = self.space_group_info().type(),
      symmetry_flags = self.symmetry_flags())
    assert self._tags.n_grid_misses()  ==  0

  def tags(self):
    return self._tags

  def peak_search(self, parameters, map, verify_symmetry = True):
    if (parameters is None):
      parameters = peak_search_parameters()
    if (verify_symmetry and libtbx.env.full_testing):
      assert self._tags.verify(map)
    if (map.accessor().is_padded()):
      map = copy(map, flex.grid(map.focus()))
    grid_peaks = peak_list(
      data = map,
      tags = self._tags.tag_array(),
      peak_search_level = parameters.peak_search_level(),
      max_peaks = parameters.max_peaks(),
      peak_cutoff = parameters.peak_cutoff(),
      interpolate = parameters.interpolate())
    if (parameters.min_distance_sym_equiv() is None):
      return grid_peaks
    return peak_cluster_analysis(
      peak_list = grid_peaks,
      special_position_settings = crystal.special_position_settings(
        crystal_symmetry = self.crystal_symmetry(),
        min_distance_sym_equiv = parameters.min_distance_sym_equiv()),
      general_positions_only = parameters.general_positions_only(),
      effective_resolution = parameters.effective_resolution(),
      significant_height_fraction = parameters.significant_height_fraction(),
      cluster_height_fraction = parameters.cluster_height_fraction(),
      min_cross_distance = parameters.min_cross_distance(),
      max_clusters = parameters.max_clusters(),
      min_cubicle_edge = parameters.min_cubicle_edge())

class boxes_by_dimension(object):
  def __init__(self,
               n_real,
               abc,
               dim,
               log = None,
               prefix = ""):
    self.n_real = n_real
    #
    step_1 = abc[0]/n_real[0] # step size along edge
    step_2 = abc[1]/n_real[1] # step size along edge
    step_3 = abc[2]/n_real[2] # step size along edge
    i_step_1 = int(dim/step_1) # points per box edge
    i_step_2 = int(dim/step_2) # points per box edge
    i_step_3 = int(dim/step_3) # points per box edge
    #
    n_boxes = self._generate_boxes(i_step_1, i_step_2, i_step_3)
    assert n_boxes  ==  len(self.starts)

  def _generate_boxes(self, ba, bb, bc):
    def regroup(be):
      maxe = be[len(be)-1][1]
      step = int(maxe/len(be))
      result = []
      for i in range(len(be)):
        if(i == 0):
          l = 0
          r = step
        elif(i == len(be)-1):
          l = i*step
          r = maxe
        else:
          l = i*step
          r = (i+1)*step
        result.append([l, r])
      return result
    be = []
    for i, b in enumerate([ba, bb, bc]):
      be_ = self._box_edges(n_real_1d = self.n_real[i], step = b)
      be_ = regroup(be_)
      be.append(be_)
    self.starts = []
    self.ends = []
    for i in be[0]:
      for j in be[1]:
        for k in be[2]:
          self.starts.append([i[0], j[0], k[0]])
          self.ends.append([i[1], j[1], k[1]])
    return len(self.starts)

  def _box_edges(self, n_real_1d, step):
    limits = []
    for i in range(0, n_real_1d, step): limits.append(i)
    limits.append(n_real_1d)
    box_1d = []
    for i in range(len(limits)):
      if(i == 0):               box_1d.append([limits[0],  limits[1]])
      elif(i!= len(limits)-1): box_1d.append([limits[i], limits[i+1]])
    return box_1d

class boxes(object):
  """
  Split box defined by n_real into boxes where each box is a fraction of the
  whole box.
  """
  def __init__(self,
               n_real,
               fraction = None,
               log = None,
               max_boxes = 2000,
               prefix = ""):
    self.n_real = n_real
    i = 0
    n_boxes = 1.e+9
    n_boxes_ = []
    while n_boxes>max_boxes:
      ba, bb, bc = \
        min(10+i, max(3, int(n_real[0]*fraction))), \
        min(10+i, max(3, int(n_real[1]*fraction))), \
        min(10+i, max(3, int(n_real[2]*fraction)))
      n_boxes = self._generate_boxes(ba, bb, bc)
      if(n_boxes_.count(n_boxes)>3): break
      n_boxes_.append(n_boxes)
      i +=  1
    assert n_boxes  ==  len(self.starts)
    if(log):
      print(prefix, "n1, n2, n3 (n_real)  :", n_real, file = log)
      print(prefix, "points per box edge:", ba, bb, bc, file = log)
      print(prefix, "number of boxes    :", len(self.starts), file = log)

  def _generate_boxes(self, ba, bb, bc):
    def regroup(be):
      maxe = be[len(be)-1][1]
      step = int(maxe/len(be))
      result = []
      for i in range(len(be)):
        if(i == 0):
          l = 0
          r = step
        elif(i == len(be)-1):
          l = i*step
          r = maxe
        else:
          l = i*step
          r = (i+1)*step
        result.append([l, r])
      return result
    be = []
    for i, b in enumerate([ba, bb, bc]):
      be_ = self._box_edges(n_real_1d = self.n_real[i], step = b)
      be_ = regroup(be_)
      be.append(be_)
    self.starts = []
    self.ends = []
    for i in be[0]:
      for j in be[1]:
        for k in be[2]:
          self.starts.append([i[0], j[0], k[0]])
          self.ends.append([i[1], j[1], k[1]])
    return len(self.starts)

  def _box_edges(self, n_real_1d, step):
    limits = []
    for i in range(0, n_real_1d, step): limits.append(i)
    limits.append(n_real_1d)
    box_1d = []
    for i in range(len(limits)):
      if(i == 0):               box_1d.append([limits[0],  limits[1]])
      elif(i!= len(limits)-1): box_1d.append([limits[i], limits[i+1]])
    return box_1d

class peak_search_parameters(object):

  def __init__(self, peak_search_level = 1,
                     max_peaks = 0,
                     peak_cutoff = None,
                     interpolate = True,
                     min_distance_sym_equiv = None,
                     general_positions_only = False,
                     effective_resolution = None,
                     significant_height_fraction = None,
                     cluster_height_fraction = None,
                     min_cross_distance = None,
                     max_clusters = None,
                     min_cubicle_edge = 5):
    adopt_init_args(self, locals(), hide = True)

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
                     general_positions_only = False,
                     effective_resolution = None,
                     significant_height_fraction = None,
                     cluster_height_fraction = None,
                     min_cross_distance = None,
                     max_clusters = None,
                     min_cubicle_edge = 5):
    if (effective_resolution is not None):
      if (significant_height_fraction is None):
          significant_height_fraction = 1/5
      if (cluster_height_fraction is None):
          cluster_height_fraction = 1/3
    if (min_cross_distance is None):
        min_cross_distance = special_position_settings.min_distance_sym_equiv()
    adopt_init_args(self, locals(), hide = True)
    assert self._min_cross_distance is not None
    self._gridding = peak_list.gridding()
    if (effective_resolution is not None):
      self._is_processed = flex.bool(peak_list.size(), False)
    else:
      self._is_processed = None
    if (   effective_resolution is not None
        or debug_peak_cluster_analysis  ==  "use_old"):
      self._site_cluster_analysis = None
    else:
      self._site_cluster_analysis = \
        self._special_position_settings.site_cluster_analysis(
          min_cross_distance = self._min_cross_distance,
          min_self_distance
             = self._special_position_settings.min_distance_sym_equiv(),
          general_positions_only = self._general_positions_only,
          min_cubicle_edge = self._min_cubicle_edge)
    self._peak_list_indices = flex.size_t()
    self._peak_list_index = 0
    self._sites = flex.vec3_double()
    self._heights = flex.double()
    self._fixed_site_indices = flex.size_t()

  def __next__(self):
    if (self._effective_resolution is not None):
      return self.next_with_effective_resolution()
    else:
      return self.next_site_cluster_analysis()

  next = __next__

  def all(self, max_clusters = None):
    if (self._effective_resolution is not None):
      return self.all_with_effective_resolution(max_clusters = max_clusters)
    else:
      return self.all_site_cluster_analysis(max_clusters = max_clusters)

  def __iter__(self):
    while 1:
      site_info = next(self)
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
    if (self._peak_list.size()  ==  0):
      return None
    return self._peak_list.heights()[0]

  def append_fixed_site(self, site, height = 0):
    if (self._site_cluster_analysis is not None):
      self._site_cluster_analysis.insert_fixed_site_frac(original_site = site)
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
      if (peak_list_index >=  self._peak_list.size()): return None
      self._peak_list_index +=  1
      site_symmetry = self._special_position_settings.site_symmetry(
        site = self._peak_list.sites()[peak_list_index])
      site = site_symmetry.exact_site()
      if (not self._site_cluster_analysis.process_site_frac(
                original_site = site,
                site_symmetry_ops = site_symmetry)): continue
      height = self._peak_list.heights()[peak_list_index]
      self._peak_list_indices.append(peak_list_index)
      self._sites.append(site)
      self._heights.append(height)
      return cluster_site_info(
        peak_list_index = peak_list_index,
        grid_index = self._peak_list.grid_indices(peak_list_index),
        grid_height = self._peak_list.grid_heights()[peak_list_index],
        site = site,
        height = height)

  def all_site_cluster_analysis(self, max_clusters = None):
    if (max_clusters is None):
      max_clusters = self._max_clusters
    assert max_clusters is not None
    while 1:
      if (self._sites.size() >=  max_clusters): break
      if (self.next_site_cluster_analysis() is None): break
    return self

  def next_with_effective_resolution(self):
    while 1:
      peak_list_index = self._peak_list_index
      if (peak_list_index >=  self._peak_list.size()): return None
      self._peak_list_index +=  1
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
          message = "This function should not be used for"
                  " processing a large number of peaks.",
          category = RuntimeWarning)
      for s in self._sites:
        dist = sgtbx.min_sym_equiv_distance_info(equiv_sites, s).dist()
        if (dist < self._min_cross_distance):
          keep = False
          break
      if (keep  ==  True):
        if (    self._effective_resolution is not None
            and (   self._heights.size()  ==  0
                 or height <   self._heights[0]
                             * self._significant_height_fraction)):
            site, height = self._accumulate_significant(
              site, height, site_symmetry, equiv_sites)
        self._peak_list_indices.append(peak_list_index)
        self._sites.append(site)
        self._heights.append(height)
        return cluster_site_info(
          peak_list_index = peak_list_index,
          grid_index = grid_index,
          grid_height = grid_height,
          site = site,
          height = height)

  def _accumulate_significant(self, site, height, site_symmetry, equiv_sites):
    unit_cell = self.special_position_settings().unit_cell()
    orth = unit_cell.orthogonalize
    frac = unit_cell.fractionalize
    sum_w_sites = matrix.col(orth(site)) * height
    sum_w = height
    height_cutoff = height * self._cluster_height_fraction
    for i in range(self._peak_list_index, self._peak_list.size()):
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
        sum_w_sites +=  matrix.col(orth(close_site)) * other_height
        sum_w +=  other_height
    return frac(sum_w_sites / sum_w), height

  def all_with_effective_resolution(self, max_clusters = None):
    if (max_clusters is None):
      max_clusters = self._max_clusters
    assert max_clusters is not None
    while 1:
      if (self._sites.size() >=  max_clusters): break
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
  small_n_real = [0, 0, 0]
  small_origin_in_large_grid = [0, 0, 0]
  small_abc = [0, 0, 0]
  sites_frac_shift = [0, 0, 0]
  for i in range(3):
    grid_step = large_ucp[i] / large_n_real[i]
    buffer = large_d_min / grid_step
    grid_min = ifloor(large_frac_min[i] * large_n_real[i] - buffer)
    grid_max = iceil(large_frac_max[i] * large_n_real[i] + buffer)
    min_grid = grid_max - grid_min + 1
    small_n_real[i] = fftpack.adjust_gridding(min_grid = min_grid, max_prime = 5)
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
    crystal_symmetry = crystal.symmetry(
      unit_cell = tuple(small_abc)+large_ucp[3:],
      space_group_symbol = "P1"),
    scatterers = work_scatterers)
  small_xray_structure.set_sites_cart(sites_cart = sites_cart_small)
  small_f_calc = small_xray_structure.structure_factors(
    d_min = large_d_min).f_calc()
  small_gridding = crystal_gridding(
    unit_cell = small_f_calc.unit_cell(),
    space_group_info = small_f_calc.space_group_info(),
    pre_determined_n_real = small_n_real)
  from cctbx import miller
  small_fft_map = miller.fft_map(
    crystal_gridding = small_gridding,
    fourier_coefficients = small_f_calc)
  small_fft_map.apply_sigma_scaling()
  small_map = small_fft_map.real_map_unpadded()
  grid_indices = grid_indices_around_sites(
    unit_cell = small_xray_structure.unit_cell(),
    fft_n_real = small_n_real,
    fft_m_real = small_n_real,
    sites_cart = sites_cart_small,
    site_radii = site_radii)
  small_copy_from_large_map = copy(
    map_unit_cell = large_density_map,
    first = small_origin_in_large_grid,
    last = matrix.col(small_origin_in_large_grid)
       + matrix.col(small_n_real)
       - matrix.col((1, 1, 1)))
  assert small_copy_from_large_map.all()  ==  small_map.all()
  corr = flex.linear_correlation(
    x = small_map.select(grid_indices),
    y = small_copy_from_large_map.select(grid_indices))
  if (not corr.is_well_defined()):
    return None
  return corr.coefficient()

def ccv(map_1, map_2, modified, centered, cutoff = None, n_bins = 10000):
  if(modified):
    map_1 = volume_scale(map = map_1, n_bins = n_bins).map_data()
    map_2 = volume_scale(map = map_2, n_bins = n_bins).map_data()
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

class spherical_variance_around_point(object):
  def __init__(self,
      real_map,
      unit_cell,
      site_cart,
      radius,
      n_points = 40,
      spline_interpolation = True,
      write_sphere_points_to_pdb_file = None):
    self.site_cart = site_cart
    self.radius = radius
    assert n_points>0
    sphere_points = []
    x, y, z = site_cart
    # reference: "Distributing many points on a sphere" by E.B. Saff and
    #     A.B.J. Kuijlaars, Mathematical Intelligencer 19.1 (1997) 5--11.
    # derived from http://packinon.sourceforge.net/py_progs/pg_saff.html
    for k in range(1, n_points+1):
      h = -1 + 2 * (k - 1) / float(n_points - 1)
      theta = math.acos(h)
      if (k  ==  1) or (k  ==  n_points):
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
      site_frac = unit_cell.fractionalize(site_cart = point)
      value_at_point = real_map.tricubic_interpolation(site_frac)
      map_values.append(value_at_point)
    self.min = flex.min(map_values)
    self.max = flex.max(map_values)
    self.mean = flex.mean(map_values)
    self.standard_deviation = map_values.standard_deviation_of_the_sample()
    if (write_sphere_points_to_pdb_file is not None):
      f = open(write_sphere_points_to_pdb_file, "w")
      for i, point in enumerate(sphere_points):
        f.write(
          "HETATM    1  O   HOH A   1     %7.3f %7.3f %7.3f  1.00 20.00\n"%
          point)
      f.close()

  def show(self, out = None, prefix = ""):
    if (out is None) : out = sys.stdout
    print("%sMap values around point [%g, %g, %g], radius = %g:" % \
      (prefix, self.site_cart[0], self.site_cart[1], self.site_cart[2],
       self.radius), file = out)
    print("%s  min = %.2f  max = %.2f  mean = %.2f  stddev = %.2f" % \
      (prefix, self.min, self.max, self.mean, self.standard_deviation), file = out)

def principal_axes_of_inertia(
    real_map,
    site_cart,
    unit_cell,
    radius):
  st = sphericity_tensor(
    map_data  = real_map,
    unit_cell = unit_cell,
    radius    = radius,
    site_frac = unit_cell.fractionalize(site_cart))
  es = adptbx.eigensystem(st)
  def center_of_mass_():
    return center_of_mass(
    map_data = real_map, unit_cell = unit_cell, cutoff = 0.1)
  def inertia_tensor():
    return st
  def eigensystem():
    return es
  return group_args(
    center_of_mass = center_of_mass_,
    inertia_tensor = inertia_tensor,
    eigensystem    = eigensystem)

class local_scale(object):
  def __init__(
        self,
        crystal_gridding,
        crystal_symmetry,
        f_map = None,
        map_data = None,
        miller_array = None,
        d_min = None): #XXX  = 1: more features and noise
    # process inputs
    assert [f_map, map_data].count(None)  ==  1
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
    for s, e in zip(b.starts, b.ends):
      box = copy(map_data, s, e)
      box.reshape(flex.grid(box.all()))
      mi, ma, me = box.as_1d().min_max_mean().as_tuple()
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
        complete_set = miller_array.complete_set(d_min = d_min)
      self.map_coefficients = complete_set.structure_factors_from_map(
        map            = self.map_result,
        use_scale      = True,
        anomalous_flag = False,
        use_sg         = False)

def sphericity_by_heuristics(
      map_data,
      unit_cell,
      center_cart,
      radius,
      s_angle_sampling_step = 20,
      t_angle_sampling_step = 20):
  points_on_sphere_cart = flex.vec3_double()
  for s in range(0, 360, s_angle_sampling_step):
    for t in range(0, 360, t_angle_sampling_step):
      xc, yc, zc = scitbx.math.point_on_sphere(r = radius, s_deg = s, t_deg = t,
        center = center_cart)
      points_on_sphere_cart.append([xc, yc, zc])
  o = sphericity2(
    map_data              = map_data,
    center_cart           = center_cart,
    points_on_sphere_cart = points_on_sphere_cart,
    unit_cell             = unit_cell)
  return group_args(rho = o.rho_min_max_mean(), ccs = o.ccs_min_max_mean())

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
  for r in range(0, radius, step):
    r = r/100.
    dist.append(r)
    rho = flex.double()
    for s in range(0, 360, s_angle_sampling_step):
      for t in range(0, 360, t_angle_sampling_step):
        xc, yc, zc = scitbx.math.point_on_sphere(r = r, s_deg = s, t_deg = t,
          center = center_cart)
        xf, yf, zf = unit_cell.fractionalize([xc, yc, zc])
        rho.append(map_data.eight_point_interpolation([xf, yf, zf]))
        #rho.append(map_data.tricubic_interpolation([xf, yf, zf]))
    rho_1d.append(flex.mean(rho))
  return dist, rho_1d

class positivity_constrained_density_modification(object):
  def __init__(self, f, f_000, n_cycles = 100, resolution_factor = 0.25, d_min = None,
               crystal_gridding = None, complete_set = None):
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
    for i in range(n_cycles):
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
      self.f_mod = self.f.complete_with(other = self.f_mod, scale = True,
        replace_phases = True)
      #self.assert_equal()

  def assert_equal(self):
    from libtbx.test_utils import approx_equal
    x, y = self.f, self.f_mod
    x, y = x.common_sets(y)
    x = abs(x).data()
    y = abs(y).data()
    assert approx_equal(x, y)

def d_min_corner(map_data, unit_cell):
  max_index = flex.miller_index( [[(i-1)//2 for i in map_data.all()]] )
  return uctbx.d_star_sq_as_d(unit_cell.max_d_star_sq( max_index ))

def d_min_from_map(map_data, unit_cell, resolution_factor = 1./2.):
  a, b, c = unit_cell.parameters()[:3]
  nx, ny, nz = map_data.all()
  d1, d2, d3 = \
    a/nx/resolution_factor, \
    b/ny/resolution_factor, \
    c/nz/resolution_factor
  return max(d1, d2, d3)

def map_coefficients_to_map(map_coeffs, crystal_symmetry, n_real):
  assert isinstance(map_coeffs.data(), flex.complex_double)
  cg = crystal_gridding(
    unit_cell             = crystal_symmetry.unit_cell(),
    space_group_info      = crystal_symmetry.space_group_info(),
    pre_determined_n_real = n_real)
  fft_map = map_coeffs.fft_map(
    crystal_gridding = cg,
    symmetry_flags   = use_space_group_symmetry)
  fft_map.apply_volume_scaling()
  return fft_map.real_map_unpadded()

def map_to_map_coefficients(m, cs, d_min):
  import cctbx.miller
  fft = fftpack.real_to_complex_3d([i for i in m.all()])
  map_box = copy(
    m, flex.grid(fft.m_real()).set_focus(m.focus()))
  map_box.reshape(flex.grid(fft.m_real()).set_focus(fft.n_real()))
  map_box = fft.forward(map_box)
  box_structure_factors = structure_factors.from_map(
    unit_cell = cs.unit_cell(),
    space_group_type = cs.space_group().type(),
    anomalous_flag = False,
    d_min = d_min,
    complex_map = map_box,
    conjugate_flag = True,
    discard_indices_affected_by_aliasing = True)
  n = map_box.all()[0] * map_box.all()[1] * map_box.all()[2]
  map_coeffs = cctbx.miller.set(
    crystal_symmetry = cs,
    anomalous_flag = False,
    indices = box_structure_factors.miller_indices(),
    ).array(data = box_structure_factors.data()/n)
  return map_coeffs

def atom_radius_as_central_peak_width(element, b_iso, d_min, scattering_table):
  """
  Estimate atom radius as half-width of the central peak of Fourier image.
  """
  from cctbx import xray, miller
  dim = 40.
  cs = crystal.symmetry((dim, dim, dim, 90, 90, 90), "P 1")
  sp = crystal.special_position_settings(cs)
  sc = xray.scatterer(
    scattering_type = element,
    site            = (0, 0, 0),
    u               = adptbx.b_as_u(b_iso))
  scatterers = flex.xray_scatterer([sc])
  xrs = xray.structure(sp, scatterers)
  xrs.scattering_type_registry(table = scattering_table)
  cg = crystal_gridding(
    unit_cell         = xrs.unit_cell(),
    space_group_info  = xrs.space_group_info(),
    step              = 0.1)
  fc = xrs.structure_factors(d_min = d_min, algorithm = "direct").f_calc()
  fft_map = miller.fft_map(
    crystal_gridding     = cg,
    fourier_coefficients = fc,
    f_000                = xrs.f_000())
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()
  def search_curve(map_data, dim):
    x = 0.
    step = 0.01
    mv_max = None
    mv_prev = None
    while x<= dim:
      mv = map_data.eight_point_interpolation([x/dim, 0, 0])
      if(mv_prev is not None and mv>mv_prev): return x-step
      if(mv_max is None): mv_max = mv
      if(mv_max/mv>100.): return x-step
      if(mv<0.):          return x-step
      x+= step
      mv_prev = mv
    return None
  radius = search_curve(map_data = map_data, dim = dim)
  assert radius is not None
  return radius

class atom_curves(object):
  """
Class-toolkit to compute various 1-atom 1D curves: exact electron density,
Fourier image of specified resolution, etc.
  """

  def __init__(self, scattering_type, scattering_table = "wk1995",
               scattering_dictionary=None):
    adopt_init_args(self, locals())
    assert [self.scattering_table, self.scattering_dictionary].count(None)==1
    self.scr = self.get_xray_structure(box = 1, b = 0).scattering_type_registry()
    self.uff = self.scr.unique_form_factors_at_d_star_sq

  def get_xray_structure(self, box, b):
    cs = crystal.symmetry((box, box, box, 90, 90, 90), "P 1")
    sp = crystal.special_position_settings(cs)
    from cctbx import xray
    sc = xray.scatterer(
      scattering_type = self.scattering_type,
      site            = (0, 0, 0),
      u               = adptbx.b_as_u(b))
    scatterers = flex.xray_scatterer([sc])
    xrs = xray.structure(sp, scatterers)
    if(self.scattering_table is not None):
      xrs.scattering_type_registry(table = self.scattering_table)
    else:
      xrs.scattering_type_registry(custom_dict = self.scattering_dictionary)
    return xrs

  def exact_density_at_r(self, r, b_iso):
    return self.scr.gaussian(self.scattering_type).electron_density(r, b_iso)

  def exact_gradient_at_r(self, r, t, t0, b_iso):
    return self.scr.gaussian(self.scattering_type).gradient(r = r, t = t, t0 = t0,
      b_iso = b_iso)

  def exact_density(self, b_iso, radius_max = 5., radius_step = 0.001):
    r = 0.0
    density = flex.double()
    radii   = flex.double()
    ed = self.scr.gaussian(self.scattering_type)
    while r < radius_max:
      density.append(ed.electron_density(r, b_iso))
      radii.append(r)
      r+= radius_step
    return group_args(radii = radii, density = density)

  def form_factor(self, ss, b_iso):
    dss = 4*ss
    return self.uff(dss)[0]*math.exp(-b_iso*ss)

  def integrand(self, r, b_iso):
    def compute(s):
      ss = (s/2)**2
      if(abs(r)>1.e-9):
        return 2/r * s * self.form_factor(ss, b_iso) * math.sin(2*math.pi*r*s)
      else:
        return 4*math.pi * s**2 * self.form_factor(ss, b_iso)
    return compute

  def bcr_approx(self,
                 d_min,
                 radius_max,
                 radius_step,
                 mxp,
                 epsc,
                 epsp  = 0.000,
                 edist = 1.0E-13,
                 kpres = 1,
                 kprot = 112,
                 ):
    b_iso = 0 # Must always be 0! All image vals below are for b_iso=0 !!!
    from cctbx.maptbx.bcr import bcr
    im = self.image(
      d_min=d_min, b_iso=0, radius_max=radius_max, radius_step=radius_step)
    bpeak, cpeak, rpeak, _,_,_,_ = bcr.get_BCR(
      dens  = im.image_values,
      dist  = im.radii,
      dmax  = radius_max,
      mxp   = mxp,
      epsc  = epsc,
      epsp  = epsp,
      edist = edist,
      kpres = kpres,
      kprot = kprot,
      )
    #
    bcr_approx_values = flex.double()
    # FILTER
    bpeak_, cpeak_, rpeak_ = [],[],[]
    for bi, ci, ri in zip(bpeak, cpeak, rpeak):
      if(abs(bi)<1.e-6 or abs(ci)<1.e-6): continue
      else:
        bpeak_.append(bi)
        cpeak_.append(ci)
        rpeak_.append(ri)
    bpeak, cpeak, rpeak = bpeak_, cpeak_, rpeak_
    #
    for r in im.radii:
      first = 0
      second = 0
      for B, C, R in zip(bpeak, cpeak, rpeak):
        if(abs(R)<1.e-6):
          first += bcr.gauss(B=B, C=C, r=r, b_iso=0)
        else:
          second += C*bcr.chi(B=B, R=R, r=r, b_iso=0)
      bcr_approx_values.append(first + second)
    return group_args(
      radii             = im.radii,
      image_values      = im.image_values,
      bcr_approx_values = bcr_approx_values,
      B=bpeak, C=cpeak, R=rpeak)

  def image(self,
            d_min,
            b_iso,
            d_max = None,
            radius_min = 0,
            radius_max = 5.,
            radius_step = 0.001,
            n_integration_steps = 2000):
    r = radius_min
    assert d_max !=  0.
    if(d_max is None): s_min = 0
    else:              s_min = 1./d_max
    assert d_min !=  0.
    s_max = 1./d_min
    image_values = flex.double()
    radii        = flex.double()
    while r < radius_max:
      s = scitbx.math.simpson(
        f = self.integrand(r, b_iso), a = s_min, b = s_max, n = n_integration_steps)
      image_values.append(s)
      radii.append(r)
      r+= radius_step
    # Fine first inflection point
    first_inflection_point = None
    i_first_inflection_point = None
    size = image_values.size()
    second_derivatives = flex.double()
    for i in range(size):
      if(i>0 and i<size-1):
        dxx = image_values[i-1]+image_values[i+1]-2*image_values[i]
      elif(i == 0):
        dxx = 2*image_values[i+1]-2*image_values[i]
      else:
        dxx = second_derivatives[i-1]*radius_step**2
      if(first_inflection_point is None and dxx>0):
        first_inflection_point = (radii[i-1]+radii[i])/2.
        i_first_inflection_point = i
      second_derivatives.append(dxx/radius_step**2)
    return group_args(
      radii                    = radii,
      image_values             = image_values,
      first_inflection_point   = first_inflection_point,
      i_first_inflection_point = i_first_inflection_point,
      radius                   = first_inflection_point*2,
      second_derivatives       = second_derivatives)

  def image_from_miller_indices(self, miller_indices, b_iso, uc,
                                radius_max, radius_step):
    p2 = flex.double()
    tmp = flex.double()
    for mi in miller_indices:
      p2.append(self.form_factor(ss = uc.d_star_sq(mi)/4, b_iso = b_iso))
      tmp.append( 2*math.pi*mi[2] )
    mv  = flex.double()
    rad = flex.double()
    z = 0.0
    while z < radius_max:
      result = 0
      for mi, p2i, tmpi in zip(miller_indices, p2, tmp):
        result +=  p2i*math.cos(tmpi*z)
      rad.append(z)
      mv.append(result*2)
      z+= radius_step
    return group_args(radii = rad, image_values = mv/uc.volume())

  def image_from_3d(self, box, b, step, unit_cell, space_group_info,
                          miller_array):
    from cctbx import miller
    xrs = self.get_xray_structure(box = box, b = b)
    fc = miller_array.structure_factors_from_scatterers(
      xray_structure = xrs, algorithm = "direct").f_calc()
    cg = crystal_gridding(
      unit_cell         = unit_cell,
      space_group_info  = space_group_info,
      step              = step,
      symmetry_flags    = use_space_group_symmetry)
    fft_map = miller.fft_map(
      crystal_gridding     = cg,
      fourier_coefficients = fc)
    fft_map.apply_volume_scaling()
    map_data = fft_map.real_map_unpadded()
    mv = flex.double()
    radii = flex.double()
    r = 0
    while r < box:
      mv_ = map_data.eight_point_interpolation([r/box, 0, 0])
      mv.append(mv_)
      radii.append(r)
      r+= step
    return group_args(radii = radii, image_values = mv)

  def one_gaussian_exact(self, r, A0, B0, b = 0):
    cmn = 4*math.pi/(B0+b)
    return A0*cmn**1.5 * math.exp(-math.pi*cmn*r**2)

  def one_gaussian_approximation(self, d_min, b, use_inflection_point = True):
    ib0 = self.image(
      d_min = d_min, b_iso = 0, radius_max = 5, radius_step = 0.01)
    if(use_inflection_point):
      i_cut = ib0.i_first_inflection_point
    else:
      i_cut = None
      for i in range(ib0.radii.size()):
        if(ib0.image_values[i]<= 0):
          rad_cut = ib0.radii[i-1]
          i_cut = i-1
          break
    assert i_cut is not None
    # this gives a*exp(-b*x**2)
    r = scitbx.math.gaussian_fit_1d_analytical(
      x = ib0.radii[:i_cut], y = ib0.image_values[:i_cut])
    B0 = 4*math.pi**2/r.b
    A0 = r.a/(r.b/math.pi)**1.5
    image_approx_values = flex.double()
    for rad in ib0.radii:
      v = self.one_gaussian_exact(r = rad, A0 = A0, B0 = B0, b = b)
      image_approx_values.append(v)
    return group_args(image_b0 = ib0, image_approx_at_b = image_approx_values,
      i_cut = i_cut, n_points = ib0.radii.size())

def sharpen2(map, xray_structure, resolution, file_name_prefix):
  from cctbx import miller
  fo = miller.structure_factor_box_from_map(
    crystal_symmetry = xray_structure.crystal_symmetry(), map = map)
  #
  fc = fo.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  d_fsc_model = fc.d_min_from_fsc(other = fo, fsc_cutoff = 0.).d_min
  print("d_fsc_model:", d_fsc_model)
  #resolution = min(resolution, d_fsc_model)
  #resolution = d_fsc_model
  print(resolution, d_fsc_model)
  #
  xray_structure = xray_structure.set_b_iso(value = 0)
  fc = fo.structure_factors_from_scatterers(
    xray_structure = xray_structure).f_calc()
  d_spacings = fo.d_spacings().data()
  #
  cc = -999
  d_best = None
  data = fc.data().deep_copy()
  for d in [i/10. for i in range(10, 100)]:
    sel = d_spacings<d
    data_ = data.set_selected(sel, 0)
    fc_ = fc.customized_copy(data = data_)
    cc_ = fo.map_correlation(other = fc_)
    if(cc_>cc):
      cc = cc_
      d_best = d
    #print "%8.1f %10.6f"%(d, cc_)
  print("Best d:", d_best)
  #
  fc1 = xray_structure.structure_factors(d_min = resolution).f_calc()
  fc2 = fc1.resolution_filter(d_min = d_best)
  cg = crystal_gridding(
    unit_cell        = xray_structure.crystal_symmetry().unit_cell(),
    space_group_info = xray_structure.crystal_symmetry().space_group_info(),
    d_min            = resolution)
  map2 = fc2.fft_map(crystal_gridding = cg).real_map_unpadded()
  cc = -999
  b = None
  ss = 1./flex.pow2(fc1.d_spacings().data()) / 4.
  data = fc1.data()
  for b_ in range(1, 500, 1):
    xray_structure = xray_structure.set_b_iso(value = b_)
    sc = flex.exp(-b_*ss)
    fc1 = fc1.customized_copy(data = data*sc)
    map1 = fc1.fft_map(crystal_gridding = cg).real_map_unpadded()
    cc_ = flex.linear_correlation(x = map1.as_1d(), y = map2.as_1d()).coefficient()
    if(cc_>cc):
      cc = cc_
      b = b_
    #print "%8.0f %10.6f"%(b_, cc_)
  print("Best B:", b)
  #
  fo_sharp = fo.resolution_filter(d_min = resolution)
  ss = 1./flex.pow2(fo_sharp.d_spacings().data()) / 4.
  #B_sharp = -35.
  B_sharp = -1*b
  sc = flex.exp(-B_sharp*ss)
  fo_sharp = fo_sharp.customized_copy(data = fo_sharp.data()*sc)
  # output
  mtz_dataset = fo_sharp.as_mtz_dataset(column_root_label = "F")
  mtz_object = mtz_dataset.mtz_object()
  mtz_object.write(file_name = "%s.mtz"%file_name_prefix)
  fft_map = fo_sharp.fft_map(crystal_gridding = cg)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  #
  import mmtbx.masks
  mask_object = mmtbx.masks.smooth_mask(
    xray_structure = xray_structure,
    n_real         = map_data.all(),
    rad_smooth     = 2.0)
  map_data = map_data * mask_object.mask_smooth
  #
  from iotbx import mrcfile
  mrcfile.write_ccp4_map(
    file_name = "%s.ccp4"%file_name_prefix,
    unit_cell = cg.unit_cell(),
    space_group = cg.space_group(),
    #gridding_first = (0, 0, 0), # This causes a bug (map gets shifted)
    #gridding_last = n_real,  # This causes a bug (map gets shifted)
    map_data = map_data,
    labels = flex.std_string([""]))
  return fo_sharp, map_data

def loc_res(map_model_manager,
            chunk_size = 10,
            method = "fsc",
            b_min = 0,
            b_max = 500,
            b_step = 5,
            res_min = 1.5,
            res_max = 10.0,
            res_step = 0.1,
            fsc_cutoff = 0.143):
  assert method in ["fsc", "rscc", "rscc_d_min_b", "d99"]
  from cctbx import miller
  mmm = map_model_manager
  mmm.map_manager().set_mean_zero_sd_one()
  mmm.model().setup_scattering_dictionaries(scattering_table = "electron")
  result = flex.double(mmm.model().size())
  chunk_selections = mmm.model().get_hierarchy().chunk_selections(
    residues_per_chunk = chunk_size)
  for chunk_sel in chunk_selections:
    box_mmm = mmm.extract_all_maps_around_model(
       selection = chunk_sel,  box_cushion=3)
    box_mmm.remove_origin_shift_and_unit_cell_crystal_symmetry()

    # <<<<<<<<<<<<<<<
    if 0: ######## WHY THIS:
      box_mmm.mask_all_maps_around_atoms(soft_mask = True,
        mask_atoms_atom_radius=3, soft_mask_radius = 3)
    else: ######## IS NOT THE SAME AS THIS:
      box_mmm.map_manager().create_mask_around_atoms(model = box_mmm.model(),
        mask_atoms_atom_radius = 3.)
      box_mmm.map_manager().soft_mask(soft_mask_radius = 3)
      box_mmm.map_manager().apply_mask()
    # <<<<<<<<<<<<<<<

    # <<<<<<<<<<<<<<<
    if 1: ######## WHY THIS:
      fo = miller.structure_factor_box_from_map(
        crystal_symmetry = box_mmm.model().get_xray_structure().crystal_symmetry(),
        map              = box_mmm.map_manager().map_data())
    else: ######## IS NOT THE SAME AS THIS:
      fo = box_mmm.map_as_fourier_coefficients()
    # <<<<<<<<<<<<<<<

    if method in ["rscc_d_min_b", "rscc"]:
      box_mmm.model().get_xray_structure().set_b_iso(value = 0.0)
      d_spacings = fo.d_spacings().data()
      ss = 1./flex.pow2(d_spacings) / 4.
      def rcast(x): return int(x*10)
      resolutions = flex.double(
        [i/10. for i in range(rcast(res_min), rcast(res_max), rcast(res_step))])
    if method != "d99":
      fc = fo.structure_factors_from_scatterers(
        xray_structure = box_mmm.model().get_xray_structure()).f_calc()
    if(method == "d99"):
      d_min = d99(f_map=fo).result.d99
    elif(method == "fsc"):
      d_min = fc.d_min_from_fsc(other = fo, fsc_cutoff = fsc_cutoff).d_min
    elif(method == "rscc"):
      d_min, cc = cc_complex_complex(
        f_1        = fo.data(),
        f_2        = fc.data(),
        d_spacings = d_spacings,
        ss         = ss,
        d_mins     = flex.double([i/10. for i in range(15, 100)]),
        b_iso      = 0)
    elif(method == "rscc_d_min_b"):
      d_min_best = None
      cc_best = -1
      for b_iso in  [ i for i in range(b_min, b_max, b_step)]:
        d_min, cc = cc_complex_complex(
          f_1        = fo.data(),
          f_2        = fc.data(),
          d_spacings = d_spacings,
          ss         = ss,
          d_mins     = flex.double([i/10. for i in range(15, 100)]),
          b_iso      = b_iso)
        if cc>cc_best:
          cc_best = cc
          d_min_best = d_min
      cc = cc_best
      d_min = d_min_best
    result.set_selected(chunk_sel, d_min)
  return group_args(result = result, selections = chunk_selections)

def is_bounded_by_constant(map_data,
     relative_sd_tol = 0.1):
    ''' Determine if this map is bounded on all sides by values that are
       zero or a constant, within relative tolerance of relative_sd_tol to
       the SD of the map as a whole

       Returns True if map boundary values are nearly constant,
         and False if they vary

       Requires that map is at origin (0,0,0)
    '''

    assert tuple(map_data.origin()) == (0,0,0)

    relative_sd = relative_sd_on_edges(map_data)

    # Determine whether values at boundaries are all about the same:
    if relative_sd  > relative_sd_tol:
      return False  # Not uniform on edges
    else:
      return True  # uniform on edges


def relative_sd_on_edges(map_data,
    skip_if_greater_than = None,
    use_maximum = None):

    '''
     Determine relative SD of values on edges to the map as a whole

     Requires that map is at origin (0,0,0)

    '''
    assert tuple(map_data.origin()) == (0,0,0)

    sd_overall = map_data.as_1d().standard_deviation_of_the_sample()

    all = list(map_data.all())
    boundary_data = flex.double()
    relative_sd_on_edges = 0

    from cctbx.maptbx import copy
    for i in (0, all[0]-1):
      new_map_data = copy(map_data,
         tuple((i, 0, 0)),
         tuple((i, all[1], all[2])))
      boundary_data.extend(new_map_data.as_1d())
      relative_sd_on_edges = max(relative_sd_on_edges,
        new_map_data.as_1d().standard_deviation_of_the_sample() / max(
          1.e-10,sd_overall))
      if (skip_if_greater_than is not None) and (
        relative_sd_on_edges > skip_if_greater_than):
          return relative_sd_on_edges

    for j in (0, all[1]-1):
      new_map_data = copy(map_data,
         tuple((0, j, 0)),
         tuple((all[0], j, all[2])))
      boundary_data.extend(new_map_data.as_1d())
      relative_sd_on_edges = max(relative_sd_on_edges,
        new_map_data.as_1d().standard_deviation_of_the_sample() / max(
          1.e-10,sd_overall))
      if (skip_if_greater_than is not None) and (
        relative_sd_on_edges > skip_if_greater_than):
          return relative_sd_on_edges

    for k in (0, all[2]-1):
      new_map_data = copy(map_data,
         tuple((0, 0, k)),
         tuple((all[0], all[1], k)))
      boundary_data.extend(new_map_data.as_1d())
      relative_sd_on_edges = max(relative_sd_on_edges,
        new_map_data.as_1d().standard_deviation_of_the_sample() / max(
          1.e-10,sd_overall))
      if (skip_if_greater_than is not None) and (
        relative_sd_on_edges > skip_if_greater_than):
          return relative_sd_on_edges

    if use_maximum: # Take maximum for any edge
      return relative_sd_on_edges
    else:  # use overall
      return boundary_data.standard_deviation_of_the_sample(
        ) / max(1.e-10,sd_overall)

def get_resolution_where_significant_data_present(ma,
   minimum_fraction_data_points=0.1):
    # Now filter ma at resolution where there are significant data
    sel = ( ma.amplitudes().data() > 1.e-10)
    ma_with_data = ma.select(sel)
    n_bins = int(0.5+10 * 1/minimum_fraction_data_points)
    ma_with_data.setup_binner(n_bins = n_bins, d_max = 10000.,
      d_min = ma_with_data.d_min())
    dsd = ma_with_data.d_spacings().data()
    ibin_list=list(ma_with_data.binner().range_used())
    ibin_list.reverse()
    total_data = ma_with_data.size()
    minimum_data_points = int(minimum_fraction_data_points * total_data)
    total_found = 0
    for i_bin in ibin_list:
      sel2      = ma_with_data.binner().selection(i_bin)
      dd        = dsd.select(sel2)
      d_max     = dd.min_max_mean().max
      n         = dd.size()
      total_found += n
      if total_found >= minimum_data_points and total_found < total_data//2:
        return d_max
    return None

def get_diff_score_towards_periodic(map_data,
      minimum_fraction_data_points = None):

    '''
      Evaluate consistency of high-pass filtered difference map analysis
      with that expected for a map that is periodic.

      The difference map is difference between the map and the map lacking high-
      resolution terms.  This difference map shows only high-frequency
      information

      A map that is periodic should give a difference map that is more or less
      uniform everywhere.  A non-periodic map should have a discontinuity at the
      borders and have high variation in the difference map at the edges.
    '''

    from cctbx import crystal
    dummy_uc_parameters=tuple(list(map_data.all())+[90.,90.,90.])
    cs= crystal.symmetry( dummy_uc_parameters, 1)

    # Normalize the map data
    sd=max(1.e-20,map_data.as_1d().standard_deviation_of_the_sample())
    mean=map_data.as_1d().min_max_mean().mean
    map_data=(map_data - mean)/sd

    # Test for difference map variation at edges of map

    # Get all structure factors, back transform to get map that can
    #  be represented by FT of all map coefficients in box (may not
    #  be the same as original because gridding may not allow it)

    from cctbx import miller
    ma = miller.structure_factor_box_from_map(
      crystal_symmetry = cs,
      map              = map_data,
      d_min            = None)

    map_data = map_coefficients_to_map(
      map_coeffs       = ma,
      crystal_symmetry = cs,
      n_real           = map_data.all())

    # Now we have map that can be represented by Fourier coefficients.
    # First get the map as Fourier coefficients

    ma = miller.structure_factor_box_from_map(
      crystal_symmetry = cs,
      map              = map_data,
      d_min            = None)

    # Ready with map as Fourier coefficients (FT of ma will give map_data again)

    # Find highest resolution where there are some non-zero data

    d_min_value = get_resolution_where_significant_data_present(ma,
      minimum_fraction_data_points = minimum_fraction_data_points)

    # High-frequency filter at this resolution
    filtered_ma = ma.resolution_filter(d_min = d_min_value)

    filtered_map       = map_coefficients_to_map(
      map_coeffs       = filtered_ma,
      crystal_symmetry = cs,
      n_real           = map_data.all())

    # Make a difference map to look at only high_frequency terms

    diff_map=map_data - filtered_map

    # Get the overall SD of the map and SD on edges:

    diff_sd = diff_map.as_1d().standard_deviation_of_the_sample()
    diff_relative_sd_on_edges = relative_sd_on_edges(diff_map,
       use_maximum = True)

    # Score based on expectation that a periodic map has a value of about 1
    #  and a non-periodic map has a value about 2

    diff_score_towards_aperiodic = max(0,min(1,(
         diff_relative_sd_on_edges - 1)/(2 - 1)))
    diff_score_towards_periodic = 1 - diff_score_towards_aperiodic

    return diff_score_towards_periodic

def get_edge_score_towards_periodic(map_data,
  use_minimum = True):
    '''
      Measure of whether facing edges have correlated data with correlation
     similar to that found for adjacent planes and different than randomly
     chosen points

     If use_minimum is set, take minimum of values on all pairs of faces

    '''

    all = list(map_data.all())
    one_data = flex.double()
    middle_plus_one_data = flex.double()
    middle_data = flex.double()
    boundary_zero_data = flex.double()
    boundary_zero_one_data = flex.double()
    boundary_one_data = flex.double()

    lowest_relative_cc = 1.0

    from cctbx.maptbx import copy
    unique_list=[]
    for i in (0,1, all[0]-1):
      if not i in unique_list: unique_list.append(i)
      new_map_data = copy(map_data,
         tuple((i, 0, 0)),
         tuple((i, all[1], all[2])))
      if i == 0:
        boundary_zero_data_local=new_map_data.as_1d()
        boundary_zero_data.extend(new_map_data.as_1d())
      elif i == 1:
        one_data_local=new_map_data.as_1d()
        one_data.extend(new_map_data.as_1d())
      else:
        boundary_one_data_local=new_map_data.as_1d()
        boundary_one_data.extend(new_map_data.as_1d())
    lowest_relative_cc = min(lowest_relative_cc,get_relative_cc(
      boundary_zero_data=boundary_zero_data_local,
      boundary_one_data=boundary_one_data_local,
      one_data=one_data_local,))

    assert len(unique_list) == 3

    unique_list=[]
    for j in (0,1, all[1]-1):
      if not j in unique_list: unique_list.append(j)
      new_map_data = copy(map_data,
         tuple((0, j, 0)),
         tuple((all[0], j, all[2])))
      if j == 0:
        boundary_zero_data_local=new_map_data.as_1d()
        boundary_zero_data.extend(new_map_data.as_1d())
      elif j == 1:
        one_data_local=new_map_data.as_1d()
        one_data.extend(new_map_data.as_1d())
      else:
        boundary_one_data_local=new_map_data.as_1d()
        boundary_one_data.extend(new_map_data.as_1d())
    assert len(unique_list) == 3
    lowest_relative_cc = min(lowest_relative_cc,get_relative_cc(
      boundary_zero_data=boundary_zero_data_local,
      boundary_one_data=boundary_one_data_local,
      one_data=one_data_local,))

    unique_list=[]
    for k in (0, 1, all[2]-1):
      if not k in unique_list: unique_list.append(k)
      new_map_data = copy(map_data,
         tuple((0, 0, k)),
         tuple((all[0], all[1], k)))
      if k == 0:
        boundary_zero_data_local=new_map_data.as_1d()
        boundary_zero_data.extend(new_map_data.as_1d())
      elif k == 1:
        one_data_local=new_map_data.as_1d()
        one_data.extend(new_map_data.as_1d())
      else:
        boundary_one_data_local=new_map_data.as_1d()
        boundary_one_data.extend(new_map_data.as_1d())
    assert len(unique_list) == 3
    lowest_relative_cc = min(lowest_relative_cc,get_relative_cc(
      boundary_zero_data=boundary_zero_data_local,
      boundary_one_data=boundary_one_data_local,
      one_data=one_data_local,))

    # use lowest value of relative_cc for any pair of faces so
    #  we can detect any faces that are trimmed

    if use_minimum:
      relative_cc=lowest_relative_cc
    else:
      relative_cc = get_relative_cc(
        boundary_zero_data=boundary_zero_data,
        boundary_one_data=boundary_one_data,
        one_data=one_data,)

    edge_score_towards_periodic = max(0,min(1,relative_cc ))

    return edge_score_towards_periodic

def get_relative_cc(
      boundary_zero_data = None,
      boundary_one_data = None,
      one_data = None):

    cc_boundary_zero_one= flex.linear_correlation(boundary_zero_data,
       boundary_one_data).coefficient()
    cc_positive_control= flex.linear_correlation(boundary_zero_data,
      one_data).coefficient()

    # Make negative control with randomized order of data
    one_data_random_perm= one_data.select(
         flex.random_permutation(len(one_data)))
    cc_negative_control = flex.linear_correlation(boundary_zero_data,
       one_data_random_perm).coefficient()


    # Expect that negative controls about zero, positive control high near 1,
    #  then cc_boundary_zero_one like negative_control means planes at
    #  boundaries differ, and cc_boundary_zero_one like positive means
    #  boundaries similar (as in wrapped)

    relative_cc = (cc_boundary_zero_one - cc_negative_control)/max(1.e-10,
         cc_positive_control - cc_negative_control)
    return relative_cc

def is_periodic(map_data,
     minimum_fraction_data_points = 0.1,
     high_confidence_delta = 0.2,
     medium_confidence_delta = 0.25):

    '''
       Determine if this map is periodic.  If values on opposite faces are
       about as similar as values on adjacent planes, it is probably periodic.

       Two tests are used: (1) correlation of facing edges of map and
        (2) test whether difference map between original and map without
        high resolution data shows most variation at edges (due to mismatch
        of edge data at facing edges of map).

       Map edge correlation score:
       Normally values on adjacent planes are very highly correlated (> 0.9)
       and random points in a map have very low correlation (< 0.1). This
       allows a test based on correlation of facing edges of a map and comparison
       to random pairs of points in map.

       Difference map score:
       If a map is boxed then if it is treated as a periodic map, there will
       be a discontinuity at the edges of the map.  This can be detected by
       calculating the Fourier transform of the high-resolution map coefficients
       for the map and detecting if this high-pass filtered map is dominated by
       features at the edge of the map.

       Returns True if periodic, False if not, and None if map gridding is
       too small (too few planes) or sampling is insufficiently fine to tell.

       Requires that map is at origin (0,0,0)
    '''

    assert tuple(map_data.origin()) == (0,0,0)


    # The difference map score is solid if > 0.8 or < 0.2. Otherwise best to
    #   combine it with edge score (correlation of edges) to get sum score

    diff_score_towards_periodic = get_diff_score_towards_periodic(map_data,
      minimum_fraction_data_points = minimum_fraction_data_points)

    if diff_score_towards_periodic >  (1 - high_confidence_delta):
      return True
    elif diff_score_towards_periodic < high_confidence_delta:
      return False

    # Get edge score and sum score now

    edge_score_towards_periodic = get_edge_score_towards_periodic(map_data)

    sum_score_towards_periodic = 0.5* (
         diff_score_towards_periodic + edge_score_towards_periodic)

    # We can be confident if sum_score is < .25 or > 0.75

    if sum_score_towards_periodic > (1 - medium_confidence_delta):
        return True
    elif sum_score_towards_periodic < medium_confidence_delta:
        return False
    else:
        return None # Really do not know

def map_values_along_line_connecting_two_points(map_data, points_cart,
      unit_cell, interpolation, step=None, n_steps=None):
  """
  Calculate interpolated map values along the line connecting two points in
  space.
  """
  assert interpolation in ["eight_point", "tricubic"]
  assert [step, n_steps].count(None)==1
  points_frac = unit_cell.fractionalize(points_cart)
  dist = unit_cell.distance(points_frac[0], points_frac[1])
  #
  def get_points(start, end, step=None, n_steps=None):
    assert [step, n_steps].count(None)==1
    dx = end[0] - start[0]
    dy = end[1] - start[1]
    dz = end[2] - start[2]
    direction = (dx, dy, dz)
    length = math.sqrt(direction[0]**2 + direction[1]**2 + direction[2]**2)
    direction_unit = (direction[0] / length, direction[1] / length, direction[2] / length)
    if n_steps is None:
      n_steps = int(max(abs(dx), abs(dy), abs(dz), length) / step)
    if step is None:
      step = length/n_steps
    points = []
    for i in range(n_steps + 1):
      point = (
        start[0] + i * step * direction_unit[0],
        start[1] + i * step * direction_unit[1],
        start[2] + i * step * direction_unit[2])
      points.append(point)
    # check end point
    end_ = points[-1]
    d = math.sqrt(
      (end_[0] - start[0])**2 +
      (end_[1] - start[1])**2 +
      (end_[2] - start[2])**2)
    if d<length: points.append(end)
    #
    return points
  #
  points = get_points(start=points_cart[0], end=points_cart[1], step=step,
    n_steps=n_steps)

  dist = flex.double()
  vals = flex.double()
  mv_max = None
  point_max = None
  for p in points:
    xp = p[0]
    yp = p[1]
    zp = p[2]
    rp = unit_cell.fractionalize([xp,yp,zp])
    d = unit_cell.distance(points_frac[0], rp)
    dist.append(d)
    pf = unit_cell.fractionalize([xp,yp,zp])
    if(interpolation=="eight_point"):
      mv = map_data.eight_point_interpolation(pf)
    else:
      mv = map_data.tricubic_interpolation(pf)
    if mv_max is None:
      mv_max = mv
      point_max = p[:]
    else:
      if mv > mv_max:
        mv_max = mv
        point_max = p[:]
    vals.append(mv)
  return group_args(dist = dist, vals = vals, point_max = point_max)

import cctbx.sgtbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.maptbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx.python_utils import dicts
from scitbx.python_utils.misc import adopt_init_args
import sys

class statistics(ext.statistics):

  def __init__(self, map):
    ext.statistics.__init__(self, map)

  def show_summary(self, f=sys.stdout):
    print >> f, "max %.6g" % (self.max())
    print >> f, "min %.6g" % (self.min())
    print >> f, "mean %.6g" % (self.mean())
    print >> f, "sigma %.6g" % (self.sigma())

def symmetry_flags(use_space_group_symmetry,
                   use_normalizer_k2l=00000,
                   use_structure_seminvariants=00000):
  return ext.symmetry_flags(use_space_group_symmetry,
                            use_normalizer_k2l,
                            use_structure_seminvariants)

use_space_group_symmetry = symmetry_flags(use_space_group_symmetry=0001)

def peak_list(data,
              tags,
              peak_search_level=1,
              max_peaks=0,
              peak_cutoff=None):
  if (peak_cutoff == None):
    return ext.peak_list(data, tags, peak_search_level, max_peaks)
  return ext.peak_list(data, tags, peak_search_level, peak_cutoff, max_peaks)

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling=0001):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy(
  to_map=structure_factors_to_map,
  from_map=structure_factors_from_map)

class crystal_gridding:

  def __init__(self, unit_cell,
                     d_min,
                     resolution_factor=1./3,
                     symmetry_flags=None,
                     space_group_info=None,
                     mandatory_factors=None,
                     max_prime=5,
                     assert_shannon_sampling=0001):
    adopt_init_args(self, locals(), hide=0001)
    assert symmetry_flags == None or mandatory_factors == None
    if (symmetry_flags != None): assert space_group_info != None
    if (symmetry_flags != None):
      self._n_real = determine_gridding(
        unit_cell, d_min, resolution_factor,
        symmetry_flags, space_group_info.type(),
        max_prime, assert_shannon_sampling)
    else:
      if (mandatory_factors == None): mandatory_factors = (1,1,1)
      assert len(mandatory_factors) == 3
      self._n_real = determine_gridding(
        unit_cell, d_min, resolution_factor,
        mandatory_factors,
        max_prime, assert_shannon_sampling)

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

  def mandatory_factors(self):
    return self._mandatory_factors

  def max_prime(self):
    return self._max_prime

  def n_real(self):
    return self._n_real

  def space_group(self):
    assert self.space_group_info() != None
    return self.space_group_info().group()

  def crystal_symmetry(self):
    assert self.space_group_info() != None
    return crystal.symmetry(
      unit_cell=self.unit_cell(),
      space_group_info=self.space_group_info())

class peak_list_cluster_reduction:

  def __init__(self, peak_list,
                     special_position_settings,
                     general_positions_only=00000,
                     min_cross_distance=None,
                     max_reduced_sites=None):
    adopt_init_args(self, locals(), hide=0001)
    sites = flex.vec3_double()
    gridding = peak_list.gridding()
    for entry in peak_list.entries():
      site = [float(entry.index[i]) / gridding[i] for i in xrange(3)]
      sites.append(special_position_settings.site_symmetry(site).exact_site())
    if (min_cross_distance == None):
      min_cross_distance = special_position_settings.min_distance_sym_equiv()
    self._unreduced_indices = flex.size_t()
    self._reduced_sites = flex.vec3_double()
    for unreduced_index,site in sites.items():
      site_symmetry = special_position_settings.site_symmetry(site)
      if (general_positions_only and not site_symmetry.is_point_group_1()):
        continue
      equiv_sites = sgtbx.sym_equiv_sites(site_symmetry)
      keep = 0001
      for reduced_site in self._reduced_sites:
        dist = sgtbx.min_sym_equiv_distance_info(
          equiv_sites, reduced_site).dist()
        if (dist < min_cross_distance):
          keep = 00000
          break
      if (keep == 0001):
        self._unreduced_indices.append(unreduced_index)
        self._reduced_sites.append(site)
        if (len(self._reduced_sites) == max_reduced_sites): break

  def unreduced_peak_list(self):
    return self._peak_list

  def unreduced_indices(self):
    return self._unreduced_indices

  def reduced_sites(self):
    return self._reduced_sites

  def unreduced_index(self, reduced_index):
    return self._unreduced_indices[reduced_index]

  def reduced_site(self, reduced_index):
    return self._reduced_sites[reduced_index]

  def peak_height(self, reduced_index):
    return self._peak_list.entries()[
      self._unreduced_indices[reduced_index]].value

  def peak_heights(self):
    result = flex.double()
    entries = self._peak_list.entries()
    for i in self._unreduced_indices:
      result.append(entries[i].value)
    return result

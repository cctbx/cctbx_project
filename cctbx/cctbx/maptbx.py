import cctbx.sgtbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.maptbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx.python_utils import dicts
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
              max_peaks=0):
  return ext.peak_list(data, tags, peak_search_level, max_peaks)

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling=0001):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy()
structure_factors.to_map = structure_factors_to_map
structure_factors.from_map = structure_factors_from_map

def determine_gridding(unit_cell,
                       d_min,
                       resolution_factor=1./3,
                       symmetry_flags=None,
                       space_group_info=None,
                       mandatory_factors=None,
                       max_prime=5,
                       assert_shannon_sampling=0001):
  assert symmetry_flags == None or mandatory_factors == None
  if (symmetry_flags != None): assert space_group_info != None
  if (symmetry_flags != None):
    return ext.determine_gridding(
      unit_cell, d_min, resolution_factor,
      symmetry_flags, space_group_info.type(),
      max_prime, assert_shannon_sampling)
  if (mandatory_factors == None): mandatory_factors = (1,1,1)
  assert len(mandatory_factors) == 3
  return ext.determine_gridding(
    unit_cell, d_min, resolution_factor,
    mandatory_factors,
    max_prime, assert_shannon_sampling)

import cctbx.sgtbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.maptbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx import fftpack
from scitbx.python_utils import dicts, list_algebra
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
  if (mandatory_factors == None):
    mandatory_factors = (1,1,1)
  if (symmetry_flags != None):
    if (symmetry_flags.use_structure_seminvariants()):
      mandatory_factors = space_group_info.structure_seminvariant().gridding()
    sub_space_group = symmetry_flags.select_sub_space_group(
      space_group_info.type())
    mandatory_factors = sub_space_group.refine_gridding(mandatory_factors)
  assert len(mandatory_factors) == 3
  assert d_min > 0
  if (assert_shannon_sampling): assert resolution_factor <= 0.5
  grid = unit_cell.max_miller_indices(d_min * 2 * resolution_factor)
  grid = [2 * n + 1 for n in grid]
  grid = fftpack.adjust_gridding_triple(grid, max_prime, mandatory_factors)
  if (symmetry_flags == None): return grid
  best_size = None
  ss = space_group_info.structure_seminvariant()
  g_limit = max(grid) + 1
  for g0 in xrange(grid[0], g_limit, mandatory_factors[0]):
    for g1 in xrange(grid[1], g_limit, mandatory_factors[1]):
      for g2 in xrange(grid[2], g_limit, mandatory_factors[2]):
        trial_grid = fftpack.adjust_gridding_triple(
          (g0, g1, g2), max_prime, mandatory_factors)
        if (symmetry_flags.use_structure_seminvariants()):
          trial_grid = ss.refine_gridding(trial_grid)
        trial_grid = sub_space_group.refine_gridding(trial_grid)
        assert fftpack.adjust_gridding_triple(
          trial_grid, max_prime, mandatory_factors) == trial_grid
        trial_size = list_algebra.product(trial_grid)
        if (best_size == None and trial_grid == grid):
          return grid
        if (best_size == None or trial_size < best_size):
          best_grid = trial_grid
          best_size = trial_size
  return best_grid

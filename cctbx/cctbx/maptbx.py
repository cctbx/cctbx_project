from scitbx.python_utils.misc import import_regular_symbols
from cctbx_boost import maptbx_ext as ext
import_regular_symbols(globals(), ext.__dict__)
del import_regular_symbols

from scitbx.python_utils import dicts

def symmetry_flags(use_space_group_symmetry,
                   use_normalizer_k2l=False,
                   use_structure_seminvariants=False):
  return ext.symmetry_flags(use_space_group_symmetry,
                            use_normalizer_k2l,
                            use_structure_seminvariants)

def determine_grid(unit_cell,
                   d_min,
                   resolution_factor=1./3,
                   max_prime=5,
                   mandatory_factors=(1,1,1)):
  return ext.determine_grid(
           unit_cell, d_min, resolution_factor, max_prime, mandatory_factors)

def peak_list(data,
              tags,
              peak_search_level=1,
              max_peaks=0):
  return ext.peak_list(data, tags, peak_search_level, max_peaks)

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling=True):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy()
structure_factors.to_map = structure_factors_to_map
structure_factors.from_map = structure_factors_from_map

import cctbx.sgtbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.maptbx_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx.python_utils import dicts

def symmetry_flags(use_space_group_symmetry,
                   use_normalizer_k2l=00000,
                   use_structure_seminvariants=00000):
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

def as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling=0001):
  return ext.as_CObjectZYX(map_unit_cell, first, last, apply_sigma_scaling)

structure_factors = dicts.easy()
structure_factors.to_map = structure_factors_to_map
structure_factors.from_map = structure_factors_from_map

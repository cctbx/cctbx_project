import cctbx.maptbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.translation_search_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

def symmetry_flags(is_isotropic_search_model,
                   have_f_part):
  return ext.symmetry_flags(is_isotropic_search_model, have_f_part)

def map_gridding(unit_cell,
                 space_group_type,
                 symmetry_flags,
                 resolution_factor,
                 miller_indices_f_obs,
                 max_prime=5):
  return ext.map_gridding(unit_cell, space_group_type, symmetry_flags,
                          resolution_factor, miller_indices_f_obs, max_prime)

def fast_nv1995(gridding,
                space_group,
                anomalous_flag,
                miller_indices_f_obs,
                f_obs,
                f_part,
                miller_indices_p1_f_calc,
                p1_f_calc):
  return ext.fast_nv1995(gridding, space_group, anomalous_flag,
                         miller_indices_f_obs, f_obs, f_part,
                         miller_indices_p1_f_calc, p1_f_calc)

import cctbx.maptbx

from scitbx.python_utils import misc
ext = misc.import_ext("cctbx_boost.translation_search_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

from scitbx.boost_python_utils import injector

class symmetry_flags(ext.symmetry_flags):

  def __init__(self, is_isotropic_search_model,
                     have_f_part):
    return ext.symmetry_flags.__init__(self,
      is_isotropic_search_model, have_f_part)

class fast_terms(ext.fast_terms):

  def __init__(self, gridding,
                     anomalous_flag,
                     miller_indices_p1_f_calc,
                     p1_f_calc):
    return ext.fast_terms.__init__(self,
      gridding, anomalous_flag, miller_indices_p1_f_calc, p1_f_calc)

class _fast_terms(injector, ext.fast_terms):

  def summation(self, space_group,
                      miller_indices_f_obs,
                      m,
                      f_obs,
                      f_part,
                      squared_flag):
    return self.raw_summation(
      space_group, miller_indices_f_obs, m, f_obs, f_part, squared_flag)

class fast_nv1995(ext.fast_nv1995):

  def __init__(self, gridding,
                     space_group,
                     anomalous_flag,
                     miller_indices_f_obs,
                     f_obs,
                     f_part,
                     miller_indices_p1_f_calc,
                     p1_f_calc):
    return ext.fast_nv1995.__init__(self,
      gridding, space_group, anomalous_flag,
      miller_indices_f_obs, f_obs, f_part,
      miller_indices_p1_f_calc, p1_f_calc)

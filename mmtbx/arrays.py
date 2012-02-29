from scitbx.array_family import flex
from libtbx import adopt_init_args
import boost.python
ext = boost.python.import_ext("mmtbx_f_model_ext")
from cctbx import miller

class init(object):
  def __init__(self,
               f_obs,
               f_calc,
               f_mask,
               r_free_flags,
               k_isotropic,
               k_anisotropic,
               k_mask):
    adopt_init_args(self, locals())
    assert f_obs.indices().all_eq(f_calc.indices())
    assert f_obs.indices().all_eq(r_free_flags.indices())
    assert f_obs.indices().all_eq(f_mask.indices())
    assert k_isotropic.size() == k_anisotropic.size()
    assert k_isotropic.size() == k_mask.size()
    self.ss = 1./flex.pow2(self.f_calc.d_spacings().data()) / 4.
    self.data = ext.core(
      f_calc        = self.f_calc.data(),
      f_mask        = self.f_mask.data(),
      k_isotropic   = self.k_isotropic,
      k_anisotropic = self.k_anisotropic,
      k_mask        = self.k_mask)
    self.f_model = miller.array(miller_set=self.f_obs, data=self.data.f_model)
    self.f_model_no_aniso_scale = miller.array(
      miller_set=self.f_obs,
      data      =self.data.f_model_no_aniso_scale)
    self.selection_work = miller.array(
      miller_set=self.f_obs,
      data      =~self.r_free_flags.data())

  def select(self, selection=None):
    if(selection is None): return self
    assert self.f_obs.indices().size() == selection.size()
    return init(
      f_obs         = self.f_obs.select(selection),
      f_calc        = self.f_calc.select(selection),
      f_mask        = self.f_mask.select(selection),
      r_free_flags  = self.r_free_flags.select(selection),
      k_isotropic   = self.k_isotropic.select(selection),
      k_mask        = self.k_mask.select(selection),
      k_anisotropic = self.k_anisotropic.select(selection))

  def deep_copy(self):
    return self.select(selection=flex.bool(self.f_obs.indices().size(), True))

  def update(self,
             f_obs=None,
             f_calc=None,
             f_mask=None,
             r_free_flags=None,
             k_isotropic=None,
             k_mask=None,
             k_anisotropic=None):
    if(f_obs is None):         f_obs         = self.f_obs
    if(f_calc is None):        f_calc        = self.f_calc
    if(f_mask is None):        f_mask        = self.f_mask
    if(r_free_flags is None):  r_free_flags  = self.r_free_flags
    if(k_isotropic is None):   k_isotropic   = self.k_isotropic
    if(k_mask is None):        k_mask        = self.k_mask
    if(k_anisotropic is None): k_anisotropic = self.k_anisotropic
    return init(
      f_obs         = f_obs,
      f_calc        = f_calc,
      f_mask        = f_mask,
      r_free_flags  = r_free_flags,
      k_isotropic   = k_isotropic,
      k_mask        = k_mask,
      k_anisotropic = k_anisotropic)

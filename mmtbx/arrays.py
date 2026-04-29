"""Holder object for arrays"""
from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from libtbx import adopt_init_args
import boost_adaptbx.boost.python as bp
from six.moves import range
ext = bp.import_ext("mmtbx_f_model_ext")
from cctbx import miller

class init(object): # XXX PVA: Why Fobs, HL and r_free_flags are not here?
  def __init__(self,
               f_calc,
               f_masks=None,
               f_part1=None,
               f_part2=None,
               k_masks=None,
               k_isotropic_exp=None,
               k_isotropic=None,
               k_anisotropic=None):
    adopt_init_args(self, locals())
    if(self.f_masks is None):
      self.f_masks = [self.f_calc.customized_copy(
        data=flex.complex_double(f_calc.data().size(), 0))]
    else:
      if(not (type(self.f_masks) in [list, tuple])):
        self.f_masks = [self.f_masks]
      for fm in self.f_masks:
        assert self.f_calc.indices().all_eq(fm.indices())
    if(self.k_isotropic_exp is not None):
      assert self.k_isotropic_exp.size() == self.f_calc.indices().size()
    else: self.k_isotropic_exp = flex.double(f_calc.data().size(), 1)
    if(self.k_isotropic is not None):
      assert self.k_isotropic.size() == self.f_calc.indices().size()
    else: self.k_isotropic = flex.double(f_calc.data().size(), 1)
    if(self.k_anisotropic is not None):
      assert self.k_anisotropic.size() == self.f_calc.indices().size()
    else: self.k_anisotropic = flex.double(f_calc.data().size(), 1)
    if(self.k_masks is None):
      n=len(self.f_masks)
      self.k_masks = [flex.double(f_calc.data().size(), 0)]*n
    else:
      if(not (type(self.k_masks) in [list, tuple])):
        self.k_masks = [self.k_masks]
    if(self.f_part1 is not None):
      assert self.f_calc.indices().all_eq(self.f_part1.indices())
    else:
      self.f_part1 = self.f_calc.customized_copy(
        data=flex.complex_double(f_calc.data().size(), 0))
    if(self.f_part2 is not None):
      assert self.f_calc.indices().all_eq(self.f_part2.indices())
    else:
      self.f_part2 = self.f_calc.customized_copy(
        data=flex.complex_double(f_calc.data().size(), 0))
    # assemble f_bulk
    f_bulk_data = flex.complex_double(f_calc.data().size(), 0)
    assert len(self.k_masks) == len(self.f_masks)
    for i in range(len(self.k_masks)):
      f_bulk_data += self.k_masks[i]*self.f_masks[i].data()
    #
    self.data = ext.data(
      f_calc          = self.f_calc.data(),
      f_bulk          = f_bulk_data,
      k_isotropic_exp = self.k_isotropic_exp,
      k_isotropic     = self.k_isotropic,
      k_anisotropic   = self.k_anisotropic,
      f_part1         = self.f_part1.data(),
      f_part2         = self.f_part2.data())
    self.f_model = miller.array(miller_set=self.f_calc, data=self.data.f_model)
    self.f_model_no_aniso_scale = miller.array(
      miller_set=self.f_calc,
      data      =self.data.f_model_no_aniso_scale)

  def f_mask(self):
    assert len(self.f_masks)==1
    return self.f_masks[0]

  def k_mask(self):
    assert len(self.k_masks)==1
    return self.k_masks[0]

  def select(self, selection=None):
    if(selection is None): return self
    assert self.f_calc.indices().size() == selection.size()
    f_masks = [fm.select(selection=selection) for fm in self.f_masks]
    k_masks = [km.select(selection) for km in self.k_masks]
    return init(
      f_calc          = self.f_calc.select(selection),
      f_masks         = f_masks,
      k_isotropic_exp = self.k_isotropic_exp.select(selection),
      k_isotropic     = self.k_isotropic.select(selection),
      k_masks         = k_masks,
      k_anisotropic   = self.k_anisotropic.select(selection),
      f_part1         = self.f_part1.select(selection),
      f_part2         = self.f_part2.select(selection))

  def deep_copy(self):
    return self.select(selection=flex.bool(self.f_calc.indices().size(), True))

  def __getstate__(self):
    return {"args": (
      self.f_calc,
      self.f_masks,
      self.f_part1,
      self.f_part2,
      self.k_masks,
      self.k_isotropic_exp,
      self.k_isotropic,
      self.k_anisotropic)}

  def __setstate__(self, state):
    assert len(state) == 1
    self.__init__(*state["args"])

  # XXX PVA: returns new object, not itself updated. This is misleading!
  # See fmodel_kbu where this is done rigth.
  def update(self,
             f_calc=None,
             f_masks=None,
             k_isotropic_exp=None,
             k_isotropic=None,
             k_masks=None,
             k_anisotropic=None,
             f_part1=None,
             f_part2=None):
    if(f_calc is None):          f_calc          = self.f_calc
    if(f_masks is None):         f_masks         = self.f_masks
    if(k_isotropic_exp is None): k_isotropic_exp = self.k_isotropic_exp
    if(k_isotropic is None):     k_isotropic     = self.k_isotropic
    if(k_masks is None):         k_masks         = self.k_masks
    if(k_anisotropic is None):   k_anisotropic   = self.k_anisotropic
    if(f_part1 is None):         f_part1         = self.f_part1
    if(f_part2 is None):         f_part2         = self.f_part2
    self.__init__(
      f_calc          = f_calc,
      f_masks         = f_masks,
      k_isotropic_exp = k_isotropic_exp,
      k_isotropic     = k_isotropic,
      k_masks         = k_masks,
      k_anisotropic   = k_anisotropic,
      f_part1         = f_part1,
      f_part2         = f_part2)
    return self

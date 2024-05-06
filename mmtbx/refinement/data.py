from __future__ import absolute_import, division
from cctbx.array_family import flex
from libtbx import adopt_init_args

class fs(object):
  """
  Fourier (reciprocal) space target and gradinets manager.
  'selection' selects parameters to be refined.
  """
  def __init__(self,
               fmodel,
               sites_cart = False,
               u_iso = False,
               occupancy = False,
               selection = None):
    adopt_init_args(self, locals())
    self.target_functor_xray = fmodel.target_functor(
      alpha_beta = None # XXX Check what that means
      )
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    if selection is None:
      selection = self.fmodel.xray_structure.all_selection().iselection()
    if isinstance(selection, flex.bool):
      selection = selection.iselection()
    sc = self.fmodel.xray_structure.scatterers()
    self.size = sc.size()
    if  (self.sites_cart): sc.flags_set_grad_site(iselection = selection)
    elif(self.u_iso):      sc.flags_set_grad_u_iso(iselection = selection)
    elif(self.occupancy):  sc.flags_set_grad_occupancy(iselection = selection)

  def update(self, x):
    if  (self.sites_cart):
      pass
    elif(self.u_iso):
      self.fmodel.xray_structure.set_u_iso(values = x)
    elif(self.occupancy):
      self.fmodel.xray_structure.set_occupancies(value=x)
    self.fmodel.update_xray_structure(update_f_calc = True)
    self.tg = self.target_functor_xray(compute_gradients = True)

  def target(self):
    return self.tg.target_work()

  def gradients(self):
    if self.selection is None:
      return self.tg.gradients_wrt_atomic_parameters().packed()
    else:
      g = self.tg.gradients_wrt_atomic_parameters().packed()
      if not self.sites_cart:
        result = flex.double(self.size, 0)
        result.set_selected(self.selection, g)
      else:
        assert 0 # not implemented
      return result

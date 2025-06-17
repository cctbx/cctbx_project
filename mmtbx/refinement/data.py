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
               occupancy = False):
    adopt_init_args(self, locals())
    self.scatterers = self.fmodel.xray_structure.scatterers()
    self.size = self.scatterers.size()
    self.target_functor_xray = fmodel.target_functor(
      alpha_beta = None # XXX Check what that means
      )
    self.scatterers.flags_set_grads(state=False)
    all_selection = self.fmodel.xray_structure.all_selection()
    self.selection = all_selection
    self.sites_cart, self.u_iso, self.occupancy = None, None, None

  def get_x(self):
    assert [self.sites_cart, self.u_iso, self.occupancy].count(True) == 1
    xrs = self.fmodel.xray_structure
    if   self.sites_cart: return xrs.sites_cart().as_double()
    elif self.u_iso:      return xrs.extract_u_iso_or_u_equiv()
    elif self.occupancy:  return self.scatterers.extract_occupancies()
    else: assert 0

  def _set_flags(self, scf, selection):
    assert [self.sites_cart, self.u_iso, self.occupancy].count(True) == 1

    #self.scatterers = self.fmodel.xray_structure.scatterers()

    self.scatterers.flags_set_grads(state=False)

    if selection is not None:
      assert isinstance(selection, flex.bool)
      self.selection = selection
    scf(iselection = self.selection.iselection())

  def set_refine_occupancy(self, selection = None):
    self.sites_cart, self.u_iso, self.occupancy = False, False, True
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_occupancy,
      selection = selection)

  def set_refine_u_iso(self, selection = None):
    self.sites_cart, self.u_iso, self.occupancy = False, True, False
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_u_iso,
      selection = selection)

  def set_refine_sites(self, selection = None):
    self.sites_cart, self.u_iso, self.occupancy = True, False, False
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_site,
      selection = selection)

  def update(self, x):
    xrs = self.fmodel.xray_structure
    if  (self.sites_cart): xrs.set_sites_cart(sites_cart = flex.vec3_double(x))
    elif(self.u_iso):      xrs.set_u_iso(values=x, selection=self.selection)
    elif(self.occupancy):  xrs.set_occupancies(value=x, selection=self.selection)
    self.fmodel.update_xray_structure(update_f_calc = True)
    self.tg = self.target_functor_xray(compute_gradients = True)

  def target(self):
    return self.tg.target_work()

  def gradients(self):
    if self.selection is None:
      return self.tg.gradients_wrt_atomic_parameters().packed()
      assert 0
    else:
      if not self.sites_cart:
        g = self.tg.gradients_wrt_atomic_parameters().packed()
        result = flex.double(self.size, 0)
        result.set_selected(self.selection, g)
        return result
      else:
        g = self.tg.d_target_d_site_cart()
        result = flex.vec3_double(self.size, [0,0,0])
        result.set_selected(self.selection, g)
        return result.as_double()

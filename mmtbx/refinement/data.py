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
    # joint mode: x = (sites_cart of all atoms | u_iso of all atoms)
    self.sites_and_u_iso = None

  def _check_mode(self):
    assert [self.sites_cart, self.u_iso, self.occupancy,
            self.sites_and_u_iso].count(True) == 1

  def get_x(self):
    self._check_mode()
    xrs = self.fmodel.xray_structure
    if   self.sites_cart: return xrs.sites_cart().as_double()
    elif self.u_iso:      return xrs.extract_u_iso_or_u_equiv()
    elif self.occupancy:  return self.scatterers.extract_occupancies()
    elif self.sites_and_u_iso:
      return xrs.sites_cart().as_double().concatenate(
        xrs.extract_u_iso_or_u_equiv())
    else: assert 0

  def _set_flags(self, scf, selection):
    self._check_mode()

    #self.scatterers = self.fmodel.xray_structure.scatterers()

    self.scatterers.flags_set_grads(state=False)

    if selection is not None:
      assert isinstance(selection, flex.bool)
      self.selection = selection
    scf(iselection = self.selection.iselection())

  def _set_mode(self, sites_cart=False, u_iso=False, occupancy=False,
                sites_and_u_iso=False):
    self.sites_cart      = sites_cart
    self.u_iso           = u_iso
    self.occupancy       = occupancy
    self.sites_and_u_iso = sites_and_u_iso

  def set_refine_occupancy(self, selection = None):
    self._set_mode(occupancy=True)
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_occupancy,
      selection = selection)

  def set_refine_u_iso(self, selection = None):
    self._set_mode(u_iso=True)
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_u_iso,
      selection = selection)

  def set_refine_sites(self, selection = None):
    self._set_mode(sites_cart=True)
    self._set_flags(
      scf       = self.scatterers.flags_set_grad_site,
      selection = selection)

  def set_refine_sites_and_u_iso(self, selection = None):
    """
    Refine Cartesian sites and u_iso of the selected atoms together.

    The parameter vector is in block order, the sites_cart of all atoms
    (3 per atom) followed by the u_iso of all atoms; gradients of
    unselected atoms are zero. One structure-factor gradient calculation
    serves both blocks.
    """
    self._set_mode(sites_and_u_iso=True)
    def set_both(iselection):
      self.scatterers.flags_set_grad_site(iselection = iselection)
      self.scatterers.flags_set_grad_u_iso(iselection = iselection)
    self._set_flags(scf = set_both, selection = selection)

  def update(self, x):
    xrs = self.fmodel.xray_structure
    if  (self.sites_cart): xrs.set_sites_cart(sites_cart = flex.vec3_double(x))
    elif(self.u_iso):      xrs.set_u_iso(values=x, selection=self.selection)
    elif(self.occupancy):  xrs.set_occupancies(value=x, selection=self.selection)
    elif(self.sites_and_u_iso):
      n = self.size
      xrs.set_sites_cart(sites_cart = flex.vec3_double(x[:3*n]))
      xrs.set_u_iso(values = x[3*n:], selection = self.selection)
    self.fmodel.update_xray_structure(update_f_calc = True)
    self.tg = self.target_functor_xray(compute_gradients = True)

  def target(self):
    return self.tg.target_work()

  def gradients(self):
    if self.selection is None:
      return self.tg.gradients_wrt_atomic_parameters().packed()
      assert 0
    else:
      if self.sites_and_u_iso:
        # packed gradients are per selected scatterer: Cartesian site (3)
        # then u_iso (cctbx packing order convention 2)
        g = self.tg.gradients_wrt_atomic_parameters().packed()
        m = self.selection.count(True)
        assert g.size() == 4*m
        g.reshape(flex.grid(m, 4))
        sites = flex.vec3_double(self.size, [0,0,0])
        sites.set_selected(self.selection,
          flex.vec3_double(g.matrix_copy_block(0, 0, m, 3).as_1d()))
        u = flex.double(self.size, 0)
        u.set_selected(self.selection, g.matrix_copy_column(3))
        return sites.as_double().concatenate(u)
      elif not self.sites_cart:
        g = self.tg.gradients_wrt_atomic_parameters().packed()
        result = flex.double(self.size, 0)
        result.set_selected(self.selection, g)
        return result
      else:
        g = self.tg.d_target_d_site_cart()
        result = flex.vec3_double(self.size, [0,0,0])
        result.set_selected(self.selection, g)
        return result.as_double()

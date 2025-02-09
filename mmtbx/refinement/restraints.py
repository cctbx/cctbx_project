from __future__ import absolute_import, division
from cctbx.array_family import flex
from libtbx import adopt_init_args
from cctbx import adptbx

class manager(object):
  """
  Refinement restraints: restraints manager + refinable parameter limits

  XXX LIMITED. Restraint targets not used. Selection is not used or used
      incorrectly.

  """
  def __init__(self,
               model,
               use_target=False):
    adopt_init_args(self, locals())
    self._restraints = None
    self._lower_bound = None
    self._upper_bound = None
    self._bound_flags = None
    self._use_xyz     = None
    self._use_adp     = None
    self._use_occ     = None

  def set_flags(self, use_xyz=False, use_adp=False, use_occ=False):
    self._use_xyz = use_xyz
    self._use_adp = use_adp
    self._use_occ = use_occ

  def check_flags(self):
    assert [self._use_xyz, self._use_adp, self._use_occ].count(True) == 1

  def check_selection(self, selection):
    if selection is not None:
      assert isinstance(selection, flex.bool)
      assert self.model.size() == selection.size()

  def set_use_xyz(self, selection = None, max_shift = 3.0):
    self.set_flags(use_xyz=True)
    self.check_selection(selection)
    x = flex.vec3_double( self.get_x() )
    max_shift = flex.vec3_double(x.size(), [max_shift, max_shift, max_shift])
    self._lower_bound = x.deep_copy()
    self._upper_bound = x.deep_copy()
    if selection is not None:
      self._lower_bound.set_selected(selection, x-max_shift)
      self._upper_bound.set_selected(selection, x+max_shift)
    else:
      self._lower_bound = x-max_shift
      self._upper_bound = x+max_shift
    self.x            = x.as_double()
    self._lower_bound = self._lower_bound.as_double()
    self._upper_bound = self._upper_bound.as_double()
    self._bound_flags = flex.int(x.size()*3, 2)
    if self.use_target:
      self._restraints = self.model.restraints_manager_energies_sites(
        compute_gradients=True)

  def get_x(self):
    self.check_flags()
    xrs = self.model.get_xray_structure()
    if   self._use_xyz: return xrs.sites_cart().as_double()
    elif self._use_adp: return xrs.extract_u_iso_or_u_equiv()
    elif self._use_occ: return xrs.scatterers().extract_occupancies()
    else: assert 0

  def set_use_adp(self, selection=None, b_min=1, b_max=200):
    self.set_flags(use_adp=True)
    self.check_selection(selection)
    x = self.get_x()
    u_min, u_max = adptbx.b_as_u(b_min), adptbx.b_as_u(b_max)
    self._lower_bound = x.deep_copy()
    self._upper_bound = x.deep_copy()
    if selection is not None:
      self._lower_bound.set_selected(selection, u_min)
      self._upper_bound.set_selected(selection, u_max)
    else:
      self._lower_bound = flex.double(x.size(), u_min)
      self._upper_bound = flex.double(x.size(), u_max)
    self._bound_flags = flex.int(x.size(), 2)
    if self.use_target:
      self._restraints = self.model.energies_adp(
        iso_restraints    = None,
        use_hd            = self.model.is_neutron(),
        compute_gradients = True)

  def set_use_occ(self, selection=None, q_min=0.004, q_max=1.0):
    self.set_flags(use_occ=True)
    self.check_selection(selection)
    x = self.get_x()
    self._lower_bound = x.deep_copy()
    self._upper_bound = x.deep_copy()
    if selection is not None:
      self._lower_bound.set_selected(selection, q_min)
      self._upper_bound.set_selected(selection, q_max)
    else:
      self._lower_bound = flex.double(x.size(), q_min)
      self._upper_bound = flex.double(x.size(), q_max)
    self._bound_flags = flex.int(x.size(), 2)
    self._restraints = None

  def lower_bound(self):
    return self._lower_bound

  def upper_bound(self):
    return self._upper_bound

  def bound_flags(self):
    return self._bound_flags

  def update(self, x):
    if not self.use_target: return
    if self._use_xyz:
      self.model.set_sites_cart(flex.vec3_double(x))
      self._restraints = self.model.restraints_manager_energies_sites(
        compute_gradients=True)
    elif self._use_adp:
      b_iso = x*adptbx.u_as_b(1.)
      self.model.set_b_iso(values = b_iso)
      self._restraints = self.model.energies_adp(
        iso_restraints    = None,
        use_hd            = self.model.is_neutron(),
        compute_gradients = True)


  def target(self):
    if self._restraints is None: return None
    if   self._use_xyz: return self._restraints.target
    elif self._use_adp: return self._restraints.target
    elif self._use_occ: return None
    else: assert 0

  def gradients(self):
    if self._restraints is None: return None
    if   self._use_xyz: return self._restraints.gradients.as_double()
    elif self._use_adp: return self._restraints.u_iso_gradients
    elif self._use_occ: return None
    else: assert 0

from cctbx.xray.structure_factors.misc import from_scatterers_common
import cctbx.xray.structure_factors.gradient_flags
from cctbx.xray import ext
from cctbx import miller
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time

class from_scatterers_direct(from_scatterers_common):

  def __init__(self, manager=None,
                     xray_structure=None,
                     miller_set=None,
                     d_target_d_f_calc=None,
                     gradient_flags=None):
    from_scatterers_common.__init__(self, manager, xray_structure, miller_set)
    self._d_target_d_f_calc = d_target_d_f_calc
    if (d_target_d_f_calc is None):
      d_target_d_f_calc = flex.complex_double()
    if (gradient_flags is None):
      gradient_flags = cctbx.xray.structure_factors.gradient_flags()
    timer = user_plus_sys_time()
    self._results = ext.structure_factors_direct_with_first_derivatives(
      self._miller_set.unit_cell(),
      self._miller_set.space_group(),
      self._miller_set.indices(),
      self._xray_structure.scatterers(),
      d_target_d_f_calc,
      gradient_flags)
    if (manager is not None):
      manager.estimate_time_direct.register(
        xray_structure.scatterers().size() * miller_set.indices().size(),
        timer.elapsed())

  def f_calc(self):
    return miller.array(self._miller_set, self._results.f_calc())

  def d_target_d_f_calc(self):
    return self._d_target_d_f_calc

  def d_target_d_site(self):
    d_target_d_site = self._results.d_target_d_site()
    xray_structure = self.xray_structure()
    assert d_target_d_site.size() == xray_structure.scatterers().size()
    return d_target_d_site

  def d_target_d_u_iso(self):
    d_target_d_u_iso = self._results.d_target_d_u_iso()
    xray_structure = self.xray_structure()
    assert d_target_d_u_iso.size() == xray_structure.scatterers().size()
    return d_target_d_u_iso

  def d_target_d_u_star(self):
    d_target_d_u_star = self._results.d_target_d_u_star()
    xray_structure = self.xray_structure()
    assert d_target_d_u_star.size() == xray_structure.scatterers().size()
    return d_target_d_u_star

  def d_target_d_occupancy(self):
    d_target_d_occupancy = self._results.d_target_d_occupancy()
    xray_structure = self.xray_structure()
    assert d_target_d_occupancy.size() == xray_structure.scatterers().size()
    return d_target_d_occupancy

  def d_target_d_fp(self):
    d_target_d_fp = self._results.d_target_d_fp()
    xray_structure = self.xray_structure()
    assert d_target_d_fp.size() == xray_structure.scatterers().size()
    return d_target_d_fp

  def d_target_d_fdp(self):
    d_target_d_fdp = self._results.d_target_d_fdp()
    xray_structure = self.xray_structure()
    assert d_target_d_fdp.size() == xray_structure.scatterers().size()
    return d_target_d_fdp

  def d_target_d_site_inplace_frac_as_cart(self, d_target_d_site):
    ext.structure_factors_d_target_d_site_in_place_frac_as_cart(
      self.miller_set().unit_cell(), d_target_d_site)

from cctbx.xray.structure_factors.manager import managed_calculation_base
import cctbx.xray.structure_factors.gradient_flags
from cctbx.xray import ext
from cctbx import miller
from cctbx.array_family import flex
from scitbx.python_utils.misc import user_plus_sys_time

class from_scatterers_direct(managed_calculation_base):

  def __init__(self, xray_structure,
                     miller_set,
                     manager=None):
    managed_calculation_base.__init__(self,manager, xray_structure,miller_set)
    timer = user_plus_sys_time()
    self._results = ext.structure_factors_direct_with_first_derivatives(
      self._miller_set.unit_cell(),
      self._miller_set.space_group(),
      self._miller_set.indices(),
      self._xray_structure.scatterers(),
      None,
      cctbx.xray.structure_factors.gradient_flags())
    if (manager is not None):
      manager.estimate_time_direct.register(
        xray_structure.scatterers().size() * miller_set.indices().size(),
        timer.elapsed())

  def f_calc(self):
    return miller.array(self._miller_set, self._results.f_calc())

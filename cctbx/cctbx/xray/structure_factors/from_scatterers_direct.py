from cctbx.xray.structure_factors.manager import managed_calculation_base
from cctbx.xray.structure_factors.manager import default_cos_sin_table
from cctbx.xray.structure_factors import global_counters
from cctbx.xray import ext
from cctbx import miller
from scitbx.python_utils.misc import user_plus_sys_time

class from_scatterers_direct(managed_calculation_base):

  def __init__(self, xray_structure,
                     miller_set,
                     manager=None,
                     cos_sin_table=False):
    time_all = user_plus_sys_time()
    managed_calculation_base.__init__(self,
      manager, xray_structure, miller_set, algorithm="direct")
    timer = user_plus_sys_time()
    if (manager is not None):
      cos_sin_table = manager.cos_sin_table()
    if (cos_sin_table == True):
      cos_sin_table = default_cos_sin_table
    elif (cos_sin_table == False):
      cos_sin_table = None
    if (cos_sin_table is None):
      self._results = ext.structure_factors_direct(
        self._miller_set.unit_cell(),
        self._miller_set.space_group(),
        self._miller_set.indices(),
        self._xray_structure.scatterers(),
        self._xray_structure.scattering_type_registry())
    else:
      self._results = ext.structure_factors_direct(
        cos_sin_table,
        self._miller_set.unit_cell(),
        self._miller_set.space_group(),
        self._miller_set.indices(),
        self._xray_structure.scatterers(),
        self._xray_structure.scattering_type_registry())
    if (manager is not None):
      manager.estimate_time_direct.register(
        xray_structure.scatterers().size() * miller_set.indices().size(),
        timer.elapsed())
    global_counters.calls_from_scatterers_direct += 1
    global_counters.time_from_scatterers_direct += time_all.elapsed()

  def f_calc(self):
    return miller.array(self._miller_set, self._results.f_calc())

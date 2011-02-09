from cctbx.xray.structure_factors.manager import managed_calculation_base
from cctbx.xray.structure_factors import global_counters
from cctbx.xray import ext
from cctbx import miller
from cctbx import maptbx
from libtbx.utils import user_plus_sys_time
from libtbx import introspection

class from_scatterers_fft(managed_calculation_base):

  def __init__(self, manager,
                     xray_structure,
                     miller_set,
                     algorithm="fft"):
    time_all = user_plus_sys_time()
    managed_calculation_base.__init__(self,
      manager, xray_structure, miller_set, algorithm="fft")
    assert miller_set.d_min() > manager.d_min() * (1-1e-6)
    manager.setup_fft() # before timing
    time_sampling = user_plus_sys_time()
    sampled_density = ext.sampled_model_density(
      unit_cell=xray_structure.unit_cell(),
      scatterers=xray_structure.scatterers(),
      scattering_type_registry=xray_structure.scattering_type_registry(),
      fft_n_real=manager.rfft().n_real(),
      fft_m_real=manager.rfft().m_real(),
      u_base=manager.u_base(),
      wing_cutoff=manager.wing_cutoff(),
      exp_table_one_over_step_size=manager.exp_table_one_over_step_size(),
      force_complex=manager.force_complex(),
      sampled_density_must_be_positive=
        manager.sampled_density_must_be_positive(),
      tolerance_positive_definite=manager.tolerance_positive_definite())
    time_sampling = time_sampling.elapsed()
    time_fft = user_plus_sys_time()
    if (not sampled_density.anomalous_flag()):
      sf_map = manager.rfft().forward(sampled_density.real_map())
      collect_conj = True
    else:
      sf_map = manager.cfft().backward(sampled_density.complex_map())
      collect_conj = False
    time_fft = time_fft.elapsed()
    time_from_map = user_plus_sys_time()
    self._f_calc_data = maptbx.structure_factors.from_map(
      space_group=miller_set.space_group(),
      anomalous_flag=sampled_density.anomalous_flag(),
      miller_indices=miller_set.indices(),
      complex_map=sf_map,
      conjugate_flag=collect_conj).data()
    time_from_map = time_from_map.elapsed()
    time_apply_u_extra = user_plus_sys_time()
    sampled_density.eliminate_u_extra_and_normalize(
      miller_set.indices(),
      self._f_calc_data)
    time_apply_u_extra = time_apply_u_extra.elapsed()
    introspection.virtual_memory_info().update_max()
    manager.estimate_time_fft.register(
      n_scatterers=xray_structure.scatterers().size(),
      n_miller_indices=miller_set.indices().size(),
      time_sampling=time_sampling,
      time_fft=time_fft,
      time_from_or_to_map=time_from_map,
      time_apply_u_extra=time_apply_u_extra)
    global_counters.from_scatterers_fft.process(time_all.elapsed())

  def f_calc(self):
    return miller.array(self.miller_set(), self._f_calc_data)

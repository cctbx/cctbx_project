from cctbx.xray.structure_factors.manager import managed_calculation_base
from cctbx.xray import ext
from cctbx import miller
from cctbx import maptbx
from scitbx.python_utils.misc import user_plus_sys_time

class from_scatterers_fft(managed_calculation_base):

  def __init__(self, manager,
                     xray_structure,
                     miller_set):
    managed_calculation_base.__init__(self, manager,xray_structure,miller_set)
    assert miller_set.d_min() >= manager.d_min()
    manager.setup_fft() # before timing
    time_sampling = user_plus_sys_time()
    sampled_density = ext.sampled_model_density(
      unit_cell=xray_structure.unit_cell(),
      scatterers=xray_structure.scatterers(),
      scattering_dict=xray_structure.scattering_dict(),
      fft_n_real=manager.rfft().n_real(),
      fft_m_real=manager.rfft().m_real(),
      u_base=manager.u_base(),
      wing_cutoff=manager.wing_cutoff(),
      exp_table_one_over_step_size=manager.exp_table_one_over_step_size(),
      force_complex=manager.force_complex(),
      electron_density_must_be_positive=
        manager.electron_density_must_be_positive(),
      tolerance_positive_definite=manager.tolerance_positive_definite())
    time_sampling = time_sampling.elapsed()
    time_fft = user_plus_sys_time()
    if (not sampled_density.anomalous_flag()):
      sf_map = manager.rfft().forward(sampled_density.real_map())
      collect_conj = 1
    else:
      sf_map = manager.cfft().backward(sampled_density.complex_map())
      collect_conj = 0
    time_fft = time_fft.elapsed()
    time_from_map = user_plus_sys_time()
    self._f_calc_data = maptbx.structure_factors.from_map(
      miller_set.space_group(),
      sampled_density.anomalous_flag(),
      miller_set.indices(),
      sf_map,
      collect_conj).data()
    time_from_map = time_from_map.elapsed()
    time_apply_u_extra = user_plus_sys_time()
    sampled_density.eliminate_u_extra_and_normalize(
      miller_set.indices(),
      self._f_calc_data)
    time_apply_u_extra = time_apply_u_extra.elapsed()
    manager.estimate_time_fft.register(
      n_scatterers=xray_structure.scatterers().size(),
      n_miller_indices=miller_set.indices().size(),
      time_sampling=time_sampling,
      time_fft=time_fft,
      time_from_or_to_map=time_from_map,
      time_apply_u_extra=time_apply_u_extra)

  def f_calc(self):
    return miller.array(self.miller_set(), self._f_calc_data)

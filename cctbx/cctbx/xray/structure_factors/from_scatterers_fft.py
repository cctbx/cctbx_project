from cctbx.xray.structure_factors.misc import from_scatterers_common
from cctbx.xray import ext
from cctbx import miller
from cctbx import maptbx
from scitbx import fftpack
from scitbx.python_utils.misc import user_plus_sys_time

class from_scatterers_fft(from_scatterers_common):

  def __init__(self, manager,
                     xray_structure,
                     miller_set,
                     d_target_d_f_calc=None,
                     gradient_flags=None,
                     force_complex=00000,
                     electron_density_must_be_positive=0001):
    from_scatterers_common.__init__(self, manager, xray_structure, miller_set)
    assert manager.symmetry_flags().use_space_group_symmetry()
    assert miller_set.d_min() >= manager.d_min()
    assert d_target_d_f_calc is None, "FFT derivatives not implemented."
    assert gradient_flags is None, "FFT derivatives not implemented."
    manager.setup_fft() # before timing
    time_sampling = user_plus_sys_time()
    sampled_density = ext.sampled_model_density(
      xray_structure.unit_cell(),
      xray_structure.scatterers(),
      manager.rfft().n_real(),
      manager.rfft().m_real(),
      manager.u_extra(),
      manager.wing_cutoff(),
      manager.exp_table_one_over_step_size(),
      force_complex,
      electron_density_must_be_positive)
    time_sampling = time_sampling.elapsed()
    time_symmetry_mapping = user_plus_sys_time()
    sampled_density.apply_symmetry(manager.crystal_gridding_tags().tags())
    time_symmetry_mapping = time_symmetry_mapping.elapsed()
    time_fft = user_plus_sys_time()
    if (not sampled_density.anomalous_flag()):
      map = sampled_density.real_map()
      sf_map = manager.rfft().forward(map)
      collect_conj = 1
    else:
      cfft = fftpack.complex_to_complex_3d(manager.rfft().n_real())
      map = sampled_density.complex_map()
      sf_map = cfft.backward(map)
      collect_conj = 0
    time_fft = time_fft.elapsed()
    time_collect = user_plus_sys_time()
    self._f_calc_data = maptbx.structure_factors.from_map(
      sampled_density.anomalous_flag(),
      miller_set.indices(),
      sf_map,
      collect_conj).data()
    sampled_density.eliminate_u_extra_and_normalize(
      miller_set.indices(),
      self._f_calc_data)
    time_collect = time_collect.elapsed()
    manager.estimate_time_fft.register(
      xray_structure.scatterers().size(),
      miller_set.indices().size(),
      time_sampling, time_symmetry_mapping, time_fft, time_collect)

  def f_calc(self):
    return miller.array(self.miller_set(), self._f_calc_data)

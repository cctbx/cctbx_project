from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.gradients_direct \
  import gradients_direct
from cctbx.xray.structure_factors.gradients_fft \
  import gradients_fft
from cctbx import maptbx

class gradients(manager):

  def __call__(self, xray_structure,
                     miller_set,
                     d_target_d_f_calc,
                     gradient_flags,
                     direct=00000,
                     fft=00000):
    assert direct == 00000 or fft == 00000
    if (direct == 00000 and fft == 00000):
      n_scatterers = xray_structure.scatterers().size()
      n_miller_indices = miller_set.indices().size()
      if (not self.have_good_timing_estimates()):
        # rough estimate
        if (  n_scatterers * self.space_group().order_z() * n_miller_indices
            < self.crystal_gridding().n_grid_points()):
          direct = 0001
      else:
        if (   self.estimate_time_direct(n_scatterers * n_miller_indices)
            <= self.estimate_time_fft(n_scatterers, n_miller_indices)):
          direct = 0001
    if (direct): f = gradients_direct
    else:        f = gradients_fft
    return f(
      manager=self,
      xray_structure=xray_structure,
      miller_set=miller_set,
      d_target_d_f_calc=d_target_d_f_calc,
      gradient_flags=gradient_flags)

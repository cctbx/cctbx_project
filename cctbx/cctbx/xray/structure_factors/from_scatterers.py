from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.from_scatterers_direct \
  import from_scatterers_direct
from cctbx.xray.structure_factors.from_scatterers_fft \
  import from_scatterers_fft
from cctbx import maptbx

class from_scatterers(manager):

  def __call__(self, xray_structure,
                     miller_set,
                     direct=00000,
                     fft=00000):
    assert direct == 00000 or fft == 00000
    if (direct == 00000 and fft == 00000):
      n_scatterers = xray_structure.scatterers().size()
      n_miller_indices = miller_set.indices().size()
      if (not self.have_good_timing_estimates()):
        # rough estimate
        if (  4 * n_scatterers * self.space_group().order_p() * n_miller_indices
            < self.crystal_gridding().n_grid_points()):
          direct = 0001
      else:
        if (   self.estimate_time_direct(n_scatterers * n_miller_indices)
            <= self.estimate_time_fft(n_scatterers, n_miller_indices)):
          direct = 0001
    if (direct): f = from_scatterers_direct
    else:        f = from_scatterers_fft
    return f(
      manager=self,
      xray_structure=xray_structure,
      miller_set=miller_set)

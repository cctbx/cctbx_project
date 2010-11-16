from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.from_scatterers_direct \
  import from_scatterers_direct
from cctbx.xray.structure_factors.from_scatterers_fft \
  import from_scatterers_fft

class from_scatterers(manager):
  """ Factory class for structure factor evaluations """

  def __call__(self, xray_structure,
                     miller_set,
                     algorithm=None):
    """
    Evaluate structure factors and return the result

    :Parameters:

      xray_structure : `cctbx.xray.structure`
        the X-ray structure to evaluate the structure factors of
      miller_set : `cctbx.miller.set`
        the set of miller indicies to evaluate the structure factors at
      algorithm : string
        the name of the evaluation method, either "direct", "fft", or None

    :return:
      an instance of
      `cctbx.xray.structure_factors.from_scatterers_direct` or
      `cctbx.xray.structure_factors.from_scatterers_fft`
      when C{algorithm} is respectively equal to "direct" or "fft",
      or the best suited of the two of them when C{algorithm} is None
    :rtype:
      a type s.t. an instance e provides the evaluated structure factors
      as C{e.f_calc()}
    """
    assert algorithm in ("direct", "fft", None)
    if (algorithm is None):
      n_scatterers = xray_structure.scatterers().size()
      n_miller_indices = miller_set.indices().size()
      if (not self.have_good_timing_estimates()):
        # rough estimate
        if (  4*n_scatterers*self.space_group().order_p()*n_miller_indices
            < self.crystal_gridding().n_grid_points()):
          algorithm = "direct"
      else:
        if (   self.estimate_time_direct(n_scatterers * n_miller_indices)
            <= self.estimate_time_fft(n_scatterers, n_miller_indices)):
          algorithm = "direct"
    if (algorithm == "direct"): f = from_scatterers_direct
    else:                       f = from_scatterers_fft
    return f(
      manager=self,
      xray_structure=xray_structure,
      miller_set=miller_set)

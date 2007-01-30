from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.gradients_direct \
  import gradients_direct
from cctbx.xray.structure_factors.gradients_fft \
  import gradients_fft
from cctbx import maptbx

class gradients(manager):
  """ Factory class for structure factor derivatives evaluation """

  def __call__(self, xray_structure,
                     u_iso_refinable_params,
                     miller_set,
                     d_target_d_f_calc,
                     n_parameters,
                     algorithm=None):
    """
    Evaluate structure factors derivatives and return the result
    @type xray_structure: cctbx.xray.structure
    @param xray_structure: the structure to differentiate the structure
    factors of
    @type u_iso_reinable_params: ??
    @param u_iso_reinable_params: ??
    @type miller_set: L{cctbx.miller.set}
    @param miller_set: the set of miller indicies to evaluate the structure
    factors at
    @type d_target_d_f_calc: ??
    @param d_target_d_f_calc: ??
    @type n_parameters: ??
    @param n_parameters: ??
    @type algorithm: string
    @param algorithm: the name of the evaluation method, either "direct",
    "fft", or None
    @rtype: a type s.t. an instance e provides the evaluated structure factors
    gradients as C{e.packed()}
    @return: an instance of
    L{cctbx.xray.structure_factors.gradients_direct} or
    L{cctbx.xray.structure_factors.gradients_fft} when C{algorithm} is
    respectively equal to "direct" or "fft",
    or the best suited of the two of them when C{algorithm} is None
    """
    assert algorithm in ("direct", "fft", None)
    if (algorithm is None):
      n_scatterers = xray_structure.scatterers().size()
      n_miller_indices = miller_set.indices().size()
      if (not self.have_good_timing_estimates()):
        if (  4*n_scatterers*self.space_group().order_p()*n_miller_indices
            < self.crystal_gridding().n_grid_points()):
          algorithm = "direct"
      else:
        if (   self.estimate_time_direct(n_scatterers * n_miller_indices)
            <= self.estimate_time_fft(n_scatterers, n_miller_indices)):
          algorithm = "direct"
    if (algorithm == "direct"): f = gradients_direct
    else:                       f = gradients_fft
    return f(
      manager=self,
      xray_structure=xray_structure,
      u_iso_refinable_params=u_iso_refinable_params,
      miller_set=miller_set,
      d_target_d_f_calc=d_target_d_f_calc,
      n_parameters=n_parameters)

from __future__ import absolute_import, division, print_function
from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.algorithm import algorithms

class gradients(manager):
  """ Factory class for structure factor derivatives evaluation """

  def __call__(self, xray_structure,
                     u_iso_refinable_params,
                     miller_set,
                     d_target_d_f_calc,
                     n_parameters,
                     algorithm=None,
                     extra_params=None):
    """
    Evaluate structure factors derivatives and return the result

    :Parameters:

      xray_structure : `cctbx.xray.structure`
        the structure to differentiate the structure factors of
      u_iso_refinable_params
        TODO
      miller_set : `cctbx.miller.set`
        the set of miller indices to evaluate the structure factors at
      d_target_d_f_calc : `scitbx.array_family.flex.complex_double`
        the derivative of the target wrt f_calc for each miller index
      n_parameters
        TODO
      algorithm : string
        the name of the evaluation method, either "direct", "fft", or None

      :return:
        an instance of
        `cctbx.xray.structure_factors.gradients_direct` or
        `cctbx.xray.structure_factors.gradients_fft`
        when C{algorithm} is respectively equal to "direct" or "fft",
        or the best suited of the two of them when C{algorithm} is None
      :rtype:
        a type s.t. an instance e provides the evaluated structure factors
        gradients as C{e.packed()}

    """
    assert algorithm in algorithms.keys() or algorithm is None
    if (algorithm is None):
      algorithm = "fft"
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
    f = algorithms[algorithm].gradients
    return f(
      manager=self,
      xray_structure=xray_structure,
      u_iso_refinable_params=u_iso_refinable_params,
      miller_set=miller_set,
      d_target_d_f_calc=d_target_d_f_calc,
      n_parameters=n_parameters,
      extra_params=extra_params)

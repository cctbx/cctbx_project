from __future__ import absolute_import, division, print_function
# -*- coding: utf-8 -*-
from cctbx.xray.structure_factors.manager import manager
from cctbx.xray.structure_factors.algorithm import algorithms

class from_scatterers(manager):
  """ Factory class for structure factor evaluations """

  def __call__(self, xray_structure,
                     miller_set,
                     algorithm=None,
                     extra_params=None):
    """Evaluate structure factors and return the result

    :type xray_structure: cctbx.xray.structure
    :param xray_structure: the X-ray structure to evaluate the structure factors of
    :type miller_set: cctbx.miller.set
    :param miller_set: the set of miller indices to evaluate the structure factors at
    :type algorithm: string
    :param algorithm: the name of the evaluation method, either "direct", "fft", or None

    :rtype: an instance of
      `cctbx.xray.structure_factors.from_scatterers_direct` or
      `cctbx.xray.structure_factors.from_scatterers_fft`
    :retruns: an instance e of
      `cctbx.xray.structure_factors.from_scatterers_direct` or
      `cctbx.xray.structure_factors.from_scatterers_fft`
      when C{algorithm} is respectively equal to "direct" or "fft",
      or the best suited of the two of them when C{algorithm} is None
      providing the evaluated structure factors as C{e.f_calc()}
    """
    assert algorithm in algorithms.keys() or algorithm is None
    if (algorithm is None):
      algorithm = "fft"
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
    f = algorithms[algorithm].from_scatterers
    return f(
      manager=self,
      xray_structure=xray_structure,
      miller_set=miller_set,
      algorithm=algorithm, # passing algorithm allows f to decide on CPU/GPU implementation
      extra_params=extra_params)

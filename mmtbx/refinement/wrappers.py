from __future__ import division, print_function
from libtbx import adopt_init_args
from scitbx import minimizers
import mmtbx.refinement.data
from mmtbx.refinement import calculators
from cctbx import adptbx

class unrestrained_qb_fsr(object):
  """
  Unrestrained sequential occupancy (q) and isotropic ADP (b) refinement.
  Use example: refinement of ocupancy and B-factors of newly placed water only.
  Input fmodel is changed in-place.
  """
  def __init__(self,
               fmodel,
               selection,
               refine_q = True,
               refine_b = True,
               refine_xyz = True,
               macro_cycles=3,
               max_iterations=50,
               q_min = 0.,
               q_max = 1,
               b_min = 5,
               b_max = 45,
               max_xyz_shift = 0.5,
               log   = None):
    adopt_init_args(self, locals())
    data = mmtbx.refinement.data.fs(fmodel = fmodel)
    for micro_cycle in range(macro_cycles):
      # Occupancy
      if refine_q:
        calculator = calculators.occ(
          data      = data,
          selection = selection,
          q_min     = q_min,
          q_max     = q_max).calculator()
        minimized = minimizers.lbfgsb(
            calculator     = calculator,
            max_iterations = max_iterations)
        if log is not None:
          print("occ: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # ADP
      if refine_b:
        calculator = calculators.adp(
          data      = data,
          selection = selection,
          u_min     = adptbx.b_as_u(b_min),
          u_max     = adptbx.b_as_u(b_max)).calculator()
        minimized = minimizers.lbfgsb(
            calculator     = calculator,
            max_iterations = max_iterations)
        if log is not None:
          print("adp: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # XYZ
      if refine_xyz:
        calculator = calculators.xyz(
          data      = data,
          selection = selection,
          max_shift = max_xyz_shift).calculator()
        minimized = minimizers.lbfgsb(
          calculator     = calculator,
          max_iterations = max_iterations)
        if log is not None:
          print("xyz: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)

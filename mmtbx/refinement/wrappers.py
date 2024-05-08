from __future__ import division, print_function
from libtbx import adopt_init_args
from scitbx import minimizers
from scitbx.array_family import flex
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
               macro_cycles=3,
               max_iterations=50,
               q_min = 0.,
               q_max = 1,
               b_min = 5,
               b_max = 45,
               log   = None):
    adopt_init_args(self, locals())
    for micro_cycle in range(macro_cycles):
      # Occupancy
      if refine_q:
        data = mmtbx.refinement.data.fs(
          fmodel = fmodel, occupancy=True, selection=selection)
        x = fmodel.xray_structure.scatterers().extract_occupancies()
        lower_bound = x.deep_copy()
        lower_bound.set_selected(selection, q_min)
        upper_bound = x.deep_copy()
        upper_bound.set_selected(selection, q_max)
        calculator = calculators.individual(
          data        = data,
          x           = x,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          bound_flags = flex.int(x.size(), 2))
        minimized = minimizers.lbfgsb(
            calculator     = calculator,
            max_iterations = max_iterations)
        if log is not None:
          print("occ: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)
      # ADP
      if refine_b:
        data = mmtbx.refinement.data.fs(
          fmodel = fmodel, u_iso=True, selection=selection)
        x = fmodel.xray_structure.extract_u_iso_or_u_equiv()
        lower_bound = x.deep_copy()
        lower_bound.set_selected(selection, adptbx.b_as_u(b_min))
        upper_bound = x.deep_copy()
        upper_bound.set_selected(selection, adptbx.b_as_u(b_max))
        calculator = calculators.individual(
          data        = data,
          x           = x,
          lower_bound = lower_bound,
          upper_bound = upper_bound,
          bound_flags = flex.int(x.size(), 2))
        minimized = minimizers.lbfgsb(
            calculator     = calculator,
            max_iterations = max_iterations)
        if log is not None:
          print("adp: r_work=%6.4f r_free=%6.4f"%(
            fmodel.r_work(), fmodel.r_free()), file = log)

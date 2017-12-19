from __future__ import division

from xfel.cxi.postrefinement_hybrid_rs import per_frame_helper
from scitbx.lstbx import normal_eqns_solving
from scitbx.array_family import flex
import math

from xfel.merging.algorithms.error_model.sdfac_refine_lbfgs import \
  sdfac_refinery, sdfac_refine_refltable_lbfgs

class sdfac_helper(per_frame_helper):
  def build_up(pfh, objective_only=False):
    values = pfh.parameterization(pfh.x)
    residuals = pfh.refinery.fvec_callable(values)
    pfh.reset()
    if objective_only:
      pfh.add_residuals(residuals, weights=pfh.refinery.weights)
    else:
      grad_r = pfh.refinery.jacobian_callable(values)
      jacobian = flex.double(
        flex.grid(len(pfh.refinery.bins), pfh.n_parameters))
      for j, der_r in enumerate(grad_r):
        jacobian.matrix_paste_column_in_place(der_r,j)
        #print >> pfh.out, "COL",j, list(der_r)
      pfh.add_equations(residuals, jacobian, weights=pfh.refinery.weights)
    #print >> pfh.out, "rms %10.3f"%math.sqrt(flex.mean(pfh.refinery.weights*residuals*residuals)),
    print >> pfh.out, "functional value % 20.3f"%pfh.refinery.functional(residuals),
    values.show(pfh.out)

class sdfac_refine_refltable_levmar(sdfac_refine_refltable_lbfgs):
  def run_minimzer(self, values, sels, **kwargs):
    # base class uses non-squared values, but lbfgs version refines the squares.
    values = self.parameterization(values.reference**2)
    self.refinery = sdfac_refinery(self.scaler.ISIGI, self.scaler.miller_set.indices(), sels, self.log)
    self.helper = sdfac_helper(current_x = values.reference,
                               parameterization = self.parameterization, refinery = self.refinery,
                               out = self.log )
    self.iterations = normal_eqns_solving.naive_iterations(non_linear_ls = self.helper,
                                                           step_threshold = 0.0001,
                                                           gradient_threshold = 1.E-10)
    return self

  def get_refined_params(self):
    return self.parameterization(self.helper.x)

  def apply_sd_error_params(self, data, values):
    self.refinery.apply_sd_error_params(data, values)

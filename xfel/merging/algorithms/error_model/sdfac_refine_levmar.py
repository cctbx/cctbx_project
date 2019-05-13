from __future__ import division
from __future__ import print_function

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
      # note the minus sign on the jacobian
      pfh.add_equations(residuals, -jacobian, weights=pfh.refinery.weights)
    #print >> pfh.out, "rms %10.3f"%math.sqrt(flex.mean(pfh.refinery.weights*residuals*residuals)),
    print("functional value % 20.3f"%pfh.refinery.functional(residuals), end=' ', file=pfh.out)
    values.show(pfh.out)

class sdfac_refine_refltable_levmar(sdfac_refine_refltable_lbfgs):
  def run_minimzer(self, values, sels, **kwargs):
    self.refinery = sdfac_refinery(self.scaler, self, self.scaler.miller_set.indices(), sels, self.log)
    self.helper = sdfac_helper(current_x = values.reference,
                               parameterization = self.parameterization, refinery = self.refinery,
                               out = self.log )
    self.iterations = normal_eqns_solving.levenberg_marquardt_iterations(
      non_linear_ls = self.helper,
      track_all=True,
      gradient_threshold=1e-08,
      step_threshold=1e-08,
      tau=1e-08,
      n_max_iterations=200)
    return self

  def get_refined_params(self):
    return self.parameterization(self.helper.x)

  def apply_sd_error_params(self, data, values):
    self.refinery.apply_sd_error_params(data, values)

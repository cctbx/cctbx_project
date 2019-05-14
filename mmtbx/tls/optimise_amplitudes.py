from __future__ import absolute_import, division, print_function

# Import dependencies
import scitbx.array_family.flex
import scitbx.lbfgs

import boost.python
ext = boost.python.import_ext("mmtbx_tls_optimise_amplitudes_ext")
from mmtbx_tls_optimise_amplitudes_ext import *

from libtbx import adopt_init_args


class OptimiseAmplitudes:


  def __init__(self,
      target_uijs,
      target_weights,
      base_amplitudes,
      base_uijs,
      base_atom_indices,
      dataset_hash,
      residual_uijs,
      residual_amplitudes=None,
      residual_mask=None,
      convergence_tolerance=1e-5):

    # Compatibility of parameters is checked in the optimiser so not checking here

    adopt_init_args(self, locals())

    self.initial = None
    self.result = None

  def run(self):
    self._optimise()
    return self

  def _optimise(self):

    # Setup
    self.f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
        target_uijs = self.target_uijs,
        target_weights = self.target_weights,
        base_amplitudes = self.base_amplitudes,
        base_uijs = self.base_uijs,
        base_atom_indices = self.base_atom_indices,
        dataset_hash = self.dataset_hash,
        residual_uijs = self.residual_uijs,
        )
    # Which datasets should be used for residual optimisation
    if self.residual_mask is not None:
      self.f_g_calculator.set_residual_mask(self.residual_mask)

    # Apply residual amplitudes if given
    if self.residual_amplitudes is not None:
      current_amps = self.f_g_calculator.x
      isel = scitbx.array_family.flex.size_t_range(current_amps.size()-self.residual_uijs.size(), current_amps.size())
      assert len(isel) == len(self.residual_uijs)
      assert len(isel) == len(self.residual_amplitudes)
      assert max(isel) == current_amps.size() - 1
      current_vals = current_amps.select(isel)
      assert current_vals.all_eq(1.0)
      self.f_g_calculator.x.set_selected(isel, self.residual_amplitudes)

    self.initial = self.f_g_calculator.x.deep_copy()

    # Optimise
    t_params = scitbx.lbfgs.termination_parameters(
        traditional_convergence_test_eps=self.convergence_tolerance,
        )
    minimizer = scitbx.lbfgs.run(
        target_evaluator=self.f_g_calculator,
        termination_params=t_params,
        )

    # Extract
    self.n_iter = minimizer.iter()
    self.result = self.f_g_calculator.x.deep_copy()
    self.result.set_selected(self.result < 0.0, 0.0) # Report this!?!

    # Tidy up
    del self.f_g_calculator # Delete after usage as not pickle-enabled

    return self


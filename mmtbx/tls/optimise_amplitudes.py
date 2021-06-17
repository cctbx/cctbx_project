from __future__ import absolute_import, division, print_function

# Import dependencies
import scitbx.array_family.flex
import scitbx.lbfgs

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_tls_optimise_amplitudes_ext")
from mmtbx_tls_optimise_amplitudes_ext import *

from libtbx import adopt_init_args

import collections
OptimisationWeights = collections.namedtuple(
  typename = 'OptimisationWeights',
  field_names = ('sum_of_amplitudes', 'sum_of_squared_amplitudes', 'sum_of_amplitudes_squared'),
)

class OptimiseAmplitudes:
  """
  Fit an arbitrary collection of sets of Uijs to a set of target Uijs (for multiple datasets).

  Parameters
  ----------
  target_uijs : flex.sym_mat3_double
    Array of target symmetric Uij matrices for each dataset.
    (shape: n_datasets, n_atoms)
  target_weights : flex.double
    Array of target weights during optimisation.
    (shape: n_datasets, n_atoms)
  base_amplitudes : flex.double
    Array of starting amplitudes for each of the base Uij sets.
    (length n_base)
  base_uijs : list of flex.sym_mat3_double objects
    List of sets of symmetric Uij matrices that are the "components" for which the amplitudes are being optimised.
    (length n_base)
  base_atom_indices : list of flex.size_t objects
    Array of indices which indicate which atom each of the components in base_uijs corresponds to.
    (length n_base; component lengths must match corresponding elements in base_uijs)
  base_dataset_hash : flex.size_t
    Array of indices which indicate which dataset each of the components in base_uijs corresponds to.
    (length n_base)
  atomic_uijs : flex.sym_mat3_double
    Array of symmetric Uij matrices that describe the individual Uij for each atom (in all datasets).
    (length n_atoms)
  atomic_amplitudes : flex.double or None
    Array of starting amplitudes for each of the atomic uijs. If None, starting amplitudes are set to 1.0.
    (length n_atoms)
  atomic_optimisation_mask : flex.bool or None
    Array of boolean values for which datasets to use for optimising the amplitudes of the residual level.
    (length n_datasets)
  optimisation_weights : OptimisationWeights or None
    Weights object for lasso, ridge-regression and other penalty terms weights.
  convergence_tolerance : double
    convergence cutoff for optimisation
  """

  def __init__(self,
      target_uijs,
      target_weights,
      base_amplitudes,
      base_uijs,
      base_atom_indices,
      base_dataset_hash,
      atomic_uijs,
      atomic_amplitudes=None,
      atomic_optimisation_mask=None,
      optimisation_weights=None,
      convergence_tolerance=1e-5,
      ):

    # Compatibility of parameters is checked in the optimiser so not checking here

    # Initialise weights object
    if optimisation_weights is None:
      optimisation_weights = OptimisationWeights(0.0, 0.0, 0.0)

    adopt_init_args(self, locals())

    self.initial = None
    self.result = None

  def get_f_g_calculator(self):
    """
    Initialise and return the functional and gradient calculator object
    """

    f_g_calculator = MultiGroupMultiDatasetUijAmplitudeFunctionalAndGradientCalculator(
        target_uijs = self.target_uijs,
        target_weights = self.target_weights,
        base_amplitudes = self.base_amplitudes,
        base_uijs = self.base_uijs,
        base_atom_indices = self.base_atom_indices,
        base_dataset_hash = self.base_dataset_hash,
        atomic_uijs = self.atomic_uijs,
        weight_sum_of_amplitudes = self.optimisation_weights.sum_of_amplitudes,
        weight_sum_of_squared_amplitudes = self.optimisation_weights.sum_of_squared_amplitudes,
        weight_sum_of_amplitudes_squared = self.optimisation_weights.sum_of_amplitudes_squared,
        )
    # Which datasets should be used for atomic optimisation
    if self.atomic_optimisation_mask is not None:
      f_g_calculator.set_atomic_optimisation_mask(self.atomic_optimisation_mask)

    # Apply atomic amplitudes if given
    if self.atomic_amplitudes is not None:
      current_amps = f_g_calculator.x
      isel = scitbx.array_family.flex.size_t_range(current_amps.size()-self.atomic_uijs.size(), current_amps.size())
      assert len(isel) == len(self.atomic_uijs)
      assert len(isel) == len(self.atomic_amplitudes)
      assert max(isel) == current_amps.size() - 1
      current_vals = current_amps.select(isel)
      assert current_vals.all_eq(1.0)
      f_g_calculator.x.set_selected(isel, self.atomic_amplitudes)

    return f_g_calculator

  def run(self):
    """
    Run the optimisation and return self. Optimised amplitudes are stored in self.result.
    """

    # Setup
    f_g_calculator = self.get_f_g_calculator()

    self.initial = f_g_calculator.x.deep_copy()

    # Optimise
    t_params = scitbx.lbfgs.termination_parameters(
        traditional_convergence_test_eps=self.convergence_tolerance,
        )
    e_params = scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound=True,
        )
    minimizer = scitbx.lbfgs.run(
        target_evaluator=f_g_calculator,
        termination_params=t_params,
        exception_handling_params=e_params,
        )

    # Extract
    self.n_iter = minimizer.iter()
    self.result = f_g_calculator.x.deep_copy()
    self.result.set_selected(self.result < 0.0, 0.0) # Report this!?!

    return self

from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from cctbx import miller
from cctbx.crystal import symmetry
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
import numpy as np

class scaling_result(object):
  '''Stores results of scaling of an experiment'''
  err_low_signal = 1
  err_low_correlation = 2 # low correlation between the observed and reference, or model, intensities

  def __init__(self):
    self.error = None
    self.data_count = None
    self.slope = None
    self.slope_error = None
    self.offset = None
    self.offset_error = None
    self.correlation = None

class experiment_scaler(worker):
  '''Scales experiment reflection intensities to the reference, or model, intensities'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(experiment_scaler, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Scaling; cross-correlation'

  def run(self, experiments, reflections):
    self.logger.log_step_time("SCALE_FRAMES")
    if self.params.scaling.algorithm != "mark0": # mark1 implies no scaling/post-refinement
      self.logger.log("No scaling was done")
      if self.mpi_helper.rank == 0:
        self.logger.main_log("No scaling was done")
      return experiments, reflections

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    # scale experiments, one at a time. Reject experiments that do not correlate with the reference or fail to scale.
    results = []
    slopes = []
    correlations = []
    high_res_experiments = 0
    experiments_rejected_because_of_low_signal = 0
    experiments_rejected_because_of_low_correlation_with_reference = 0

    target_symm = symmetry(unit_cell = self.params.scaling.unit_cell, space_group_info = self.params.scaling.space_group)
    for expt_id, experiment in enumerate(experiments):
      exp_reflections = reflections.select(reflections['id'] == expt_id)

      # Build a miller array for the experiment reflections
      exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      exp_intensities = miller.array(exp_miller_indices, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))

      model_intensities = self.params.scaling.i_model

      # Extract an array of HKLs from the model to match the experiment HKLs
      matching_indices = miller.match_multi_indices(miller_indices_unique = model_intensities.indices(), miller_indices = exp_intensities.indices())

      # Least squares
      weights = self.params.scaling.weights
      result = self.fit_experiment_to_reference(model_intensities, exp_intensities, matching_indices, weights)

      if result.error == scaling_result.err_low_signal:
        experiments_rejected_because_of_low_signal += 1
        continue
      elif result.error == scaling_result.err_low_correlation:
        experiments_rejected_because_of_low_correlation_with_reference += 1
        continue

      slopes.append(result.slope)
      correlations.append(result.correlation)

      if self.params.output.log_level == 0:
        self.logger.log("Experiment ID: %s; Slope: %f; Correlation %f"%(experiment.identifier, result.slope, result.correlation))

      # count high resolution experiments
      if exp_intensities.d_min() <= self.params.merging.d_min:
        high_res_experiments += 1

      # apply scale factors
      if (
          not self.params.postrefinement.enable or
          'postrefine' not in self.params.dispatch.step_list
      ):
        exp_reflections['intensity.sum.value'] *= result.slope
        exp_reflections['intensity.sum.variance'] *= (result.slope**2)
      exp_reflections['correlation'] = flex.double(len(exp_reflections), result.correlation)
      new_experiments.append(experiment)
      new_reflections.extend(exp_reflections)

    new_reflections.reset_ids()
    rejected_experiments = len(experiments) - len(new_experiments)
    assert rejected_experiments == experiments_rejected_because_of_low_signal + \
                                    experiments_rejected_because_of_low_correlation_with_reference

    reflections_removed_because_of_rejected_experiments = reflections.size() - new_reflections.size()

    self.logger.log("Experiments rejected because of low signal: %d"%experiments_rejected_because_of_low_signal)
    self.logger.log("Experiments rejected because of low correlation with reference: %d"%experiments_rejected_because_of_low_correlation_with_reference)
    self.logger.log("Reflections rejected because of rejected experiments: %d"%reflections_removed_because_of_rejected_experiments)
    self.logger.log("High resolution experiments: %d"%high_res_experiments)
    if self.params.postrefinement.enable and 'postrefine' in self.params.dispatch.step_list:
      self.logger.log("Note: scale factors were not applied, because postrefinement is enabled")

    # MPI-reduce all counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    total_experiments_rejected_because_of_low_signal                            = comm.reduce(experiments_rejected_because_of_low_signal, MPI.SUM, 0)
    total_experiments_rejected_because_of_low_correlation_with_reference        = comm.reduce(experiments_rejected_because_of_low_correlation_with_reference, MPI.SUM, 0)
    total_reflections_removed_because_of_rejected_experiments                   = comm.reduce(reflections_removed_because_of_rejected_experiments, MPI.SUM, 0)
    total_high_res_experiments                                                  = comm.reduce(high_res_experiments, MPI.SUM, 0)
    all_slopes                                                                  = comm.reduce(slopes, MPI.SUM, 0)
    all_correlations                                                            = comm.reduce(correlations, MPI.SUM, 0)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0:
      self.logger.main_log('Experiments rejected because of low signal: %d'%total_experiments_rejected_because_of_low_signal)
      self.logger.main_log('Experiments rejected because of low correlation with reference: %d'%total_experiments_rejected_because_of_low_correlation_with_reference)
      self.logger.main_log('Reflections rejected because of rejected experiments: %d'%total_reflections_removed_because_of_rejected_experiments)
      self.logger.main_log('Experiments with high resolution of %5.2f Angstrom or better: %d'%(self.params.merging.d_min, total_high_res_experiments))

      if len(all_slopes) > 0:
        stats_slope = flex.mean_and_variance(flex.double(all_slopes))
        self.logger.main_log('Average experiment scale factor wrt reference: %f'%(stats_slope.mean()))
      if len(all_correlations) > 1:
        stats_correlation = flex.mean_and_variance(flex.double(all_correlations))
        self.logger.main_log('Average experiment correlation with reference: %f +/- %f'%(
            stats_correlation.mean(), stats_correlation.unweighted_sample_standard_deviation()))

      if self.params.postrefinement.enable and 'postrefine' in self.params.dispatch.step_list:
        self.logger.main_log("Note: scale factors were not applied, because postrefinement is enabled")

    self.logger.log_step_time("SCALE_FRAMES", True)

    # Do we have any data left?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(new_experiments, new_reflections)

    return new_experiments, new_reflections

  def fit_experiment_to_reference(self,
      model_intensities,
      experiment_intensities,
      matching_indices,
      weights='unit'
  ):
    'Scale the observed intensities to the reference, or model, using a linear least squares fit.'
     # Y = offset + slope * X, where Y is I_r and X is I_o

    result = scaling_result()
    result.data_count = matching_indices.pairs().size()
    if result.data_count < 3:
      result.error = scaling_result.err_low_signal
      return result

    model_subset = []
    exp_subset = []
    exp_sigmas = []
    for pair in matching_indices.pairs():
      model_subset.append(model_intensities.data()[pair[0]])
      exp_subset.append(experiment_intensities.data()[pair[1]])
      exp_sigmas.append(experiment_intensities.sigmas()[pair[1]])
    model_subset = np.array(model_subset)
    exp_subset = np.array(exp_subset)
    exp_sigmas = np.array(exp_sigmas)

    correlation = pearsonr(exp_subset, model_subset)[0]
    if correlation < self.params.filter.outlier.min_corr:
      result.error = scaling_result.err_low_correlation
      return result
    result.correlation = correlation

    def linfunc(x, m): return x*m
    # For weighting, we need to put the Icalc values on the scale of Iobs, so
    # we do a first-pass scale with unit weights and apply it to model_subset.
    slope_unwt = curve_fit(linfunc, exp_subset, model_subset)[0][0]
    model_subset_scaled = model_subset / slope_unwt
    if weights == 'unit':
      slope = slope_unwt
    elif weights == 'icalc_sigma':
      sigma = (model_subset_scaled**2 + exp_sigmas**2)**.5
      slope = curve_fit(linfunc, exp_subset, model_subset, sigma=sigma)[0][0]
    elif weights == 'icalc':
      sigma = model_subset_scaled**.5
      slope = curve_fit(linfunc, exp_subset, model_subset, sigma=sigma)[0][0]

    result.slope = slope
    return result

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(experiment_scaler)

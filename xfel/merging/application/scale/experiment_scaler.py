from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
import math
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from cctbx import miller
from cctbx.crystal import symmetry
import sys

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
    for experiment in experiments:
      exp_reflections = reflections.select(reflections['exp_id'] == experiment.identifier)

      # Build a miller array for the experiment reflections
      exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      exp_intensities = miller.array(exp_miller_indices, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))

      model_intensities = self.params.scaling.i_model

      # Extract an array of HKLs from the model to match the experiment HKLs
      matching_indices = miller.match_multi_indices(miller_indices_unique = model_intensities.indices(), miller_indices = exp_intensities.indices())

      # Least squares
      result = self.fit_experiment_to_reference(model_intensities, exp_intensities, matching_indices)

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

      new_experiments.append(experiment)
      new_reflections.extend(exp_reflections)

    rejected_experiments = len(experiments) - len(new_experiments)
    assert rejected_experiments == experiments_rejected_because_of_low_signal + \
                                    experiments_rejected_because_of_low_correlation_with_reference

    reflections_removed_because_of_rejected_experiments = reflections.size() - new_reflections.size()

    self.logger.log("Experiments rejected because of low signal: %d"%experiments_rejected_because_of_low_signal)
    self.logger.log("Experiments rejected because of low correlation with reference: %d"%experiments_rejected_because_of_low_correlation_with_reference)
    self.logger.log("Reflections rejected because of rejected experiments: %d"%reflections_removed_because_of_rejected_experiments)
    self.logger.log("High resolution experiments: %d"%high_res_experiments)
    if self.params.postrefinement.enable:
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

      if self.params.postrefinement.enable:
        self.logger.main_log("Note: scale factors were not applied, because postrefinement is enabled")

    self.logger.log_step_time("SCALE_FRAMES", True)

    # Do we have any data left?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(new_experiments, new_reflections)

    return new_experiments, new_reflections

  def fit_experiment_to_reference(self, model_intensities, experiment_intensities, matching_indices):
    'Scale the observed intensities to the reference, or model, using a linear least squares fit.'
     # Y = offset + slope * X, where Y is I_r and X is I_o

    result = scaling_result()
    result.data_count = matching_indices.pairs().size()

    if result.data_count == 0:
      result.error = scaling_result.err_low_signal
      return result

    # Do various auxilliary summations
    sum_xx = 0.
    sum_xy = 0.
    sum_yy = 0.
    sum_x = 0.
    sum_y = 0.
    sum_w = 0.
    for pair in matching_indices.pairs():
      I_w = 1. # Use unit weights for starters
      I_r = model_intensities.data()[pair[0]]
      I_o = experiment_intensities.data()[pair[1]]

      sum_xx += I_w * I_o**2
      sum_yy += I_w * I_r**2
      sum_xy += I_w * I_o * I_r
      sum_x += I_w * I_o
      sum_y += I_w * I_r
      sum_w += I_w

    # calculate Pearson correlation coefficient between X and Y and test it
    DELTA_1 = result.data_count * sum_xx - sum_x**2
    DELTA_2 = result.data_count * sum_yy - sum_y**2
    if (abs(DELTA_1) < sys.float_info.epsilon) or (abs(DELTA_2) < sys.float_info.epsilon):
      result.error = scaling_result.err_low_signal
      return result
    result.correlation = (result.data_count * sum_xy - sum_x * sum_y) / (math.sqrt(DELTA_1) * math.sqrt(DELTA_2))
    if result.correlation < self.params.filter.outlier.min_corr:
      result.error = scaling_result.err_low_correlation
      return result

    DELTA = sum_w * sum_xx
    if abs(DELTA) < sys.float_info.epsilon:
      result.error = scaling_result.err_low_signal
      return result
    result.slope = sum_w * sum_xy / DELTA

    return result

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(experiment_scaler)

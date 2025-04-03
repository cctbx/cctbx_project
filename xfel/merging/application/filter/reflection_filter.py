from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from rstbx.dials_core.integration_core import show_observations
from cctbx import miller
from cctbx.crystal import symmetry
from six.moves import cStringIO as StringIO
import numpy as np


class reflection_filter(worker):
  '''Reject individual reflections based on various criteria'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(reflection_filter, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Filter reflections'

  def validate(self):
    filter_by_significance = 'significance_filter' in self.params.select.algorithm
    filter_by_isolation_forest = 'isolation_forest' in self.params.select.algorithm
    if filter_by_isolation_forest:
      check0 = self.params.select.reflection_filter.tail_percentile > 0
      check1 = self.params.select.reflection_filter.tail_percentile < 1
      assert check0 and check1, \
        'tail_percentile must be between 0 and 1'

      check0 = self.params.select.reflection_filter.contamination_lower > 0
      check1 = self.params.select.reflection_filter.contamination_lower < 1
      check2 = self.params.select.reflection_filter.contamination_upper > 0
      check3 = self.params.select.reflection_filter.contamination_upper < 1
      assert check0 and check1 and check2 and check3, \
        'contamination must be between 0 and 1'
    if filter_by_isolation_forest:
      check0 = self.params.select.reflection_filter.sampling_fraction > 0
      check1 = self.params.select.reflection_filter.sampling_fraction < 1
      assert check0 and check1, \
        'sampling_fraction must be between 0 and 1'

  def plot_reflections(self, experiments, reflections, tag):
    correct_info = 'miller_index' in reflections.keys()
    if correct_info == False:
      reflections['miller_index'] = reflections['miller_index_asymmetric']
    q2_rank = 1 / reflections.compute_d(experiments).as_numpy_array()**2
    if correct_info == False:
      del reflections['miller_index']
    intensity_rank = reflections['intensity.sum.value'].as_numpy_array()
    q2 = self.mpi_helper.comm.gather(q2_rank, root=0)
    intensity = self.mpi_helper.comm.gather(intensity_rank, root=0)
    if self.mpi_helper.rank == 0:
      import matplotlib.pyplot as plt
      import os
      q2 = np.concatenate(q2)
      intensity = np.concatenate(intensity)
      fig, axes = plt.subplots(1, 1, figsize=(8, 3))
      axes.scatter(q2, intensity, s=1, color=[0, 0, 0], marker='.')
      axes.set_ylabel('Intensity')
      axes.set_title(tag)
      axes.set_xlabel(r'Resolution ($\mathrm{\AA}$)')
      xticks = axes.get_xticks()
      xticks = xticks[xticks > 0]
      xticklabels = [f'{l:0.2f}' for l in 1 / np.sqrt(xticks)]
      axes.set_xticks(xticks)
      axes.set_xticklabels(xticklabels)
      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + f'_Iobs_{tag}.png'
        ))
      plt.close()

  def run(self, experiments, reflections):
    filter_by_significance = 'significance_filter' in self.params.select.algorithm
    filter_by_isolation_forest = 'isolation_forest' in self.params.select.algorithm
    # only "unit_cell" "n_obs" and "resolution" algorithms are supported
    if (not filter_by_significance) and (not filter_by_isolation_forest):
      return experiments, reflections

    n_reflections_initial = len(reflections)
    n_experiments_initial = len(experiments)
    if self.params.select.reflection_filter.do_diagnostics:
      self.plot_reflections(experiments, reflections, 'Initial')
    if filter_by_significance:
      experiments, reflections = self.apply_significance_filter(experiments, reflections)
      removed_reflections_significance_rank = n_reflections_initial - len(reflections)
      removed_experiments_significance_rank = n_experiments_initial - len(experiments)
      removed_reflections_significance = self.mpi_helper.comm.reduce(
        removed_reflections_significance_rank, self.mpi_helper.MPI.SUM, root=0
        )
      removed_experiments_significance = self.mpi_helper.comm.reduce(
        removed_experiments_significance_rank, self.mpi_helper.MPI.SUM, root=0
        )
      self.logger.log(f"Reflections rejected because of significant: {removed_reflections_significance_rank}")
      self.logger.log(f"Experiments rejected because of significant: {removed_experiments_significance_rank}")
      if self.mpi_helper.rank == 0:
        self.logger.main_log(f"Total reflections rejected because of significant: {removed_reflections_significance}")
        self.logger.main_log(f"Total experiments rejected because of significant: {removed_experiments_significance}")
      if self.params.select.reflection_filter.do_diagnostics:
        self.plot_reflections(experiments, reflections, 'After Significance Filter')
    if filter_by_isolation_forest:
      experiments, reflections = self.apply_isolation_forest(experiments, reflections)
      filter_type = 'Isolation Forest'
      if self.params.select.reflection_filter.do_diagnostics:
        self.plot_reflections(experiments, reflections, 'After Isolation Forest')

    if filter_by_isolation_forest:
      removed_reflections_filter_rank = n_reflections_initial - len(reflections)
      removed_experiments_filter_rank = n_experiments_initial - len(experiments)
      if filter_by_significance:
        removed_reflections_filter_rank -= removed_reflections_significance_rank
        removed_experiments_filter_rank -= removed_experiments_significance_rank
      removed_reflections_filter = self.mpi_helper.comm.reduce(
        removed_reflections_filter_rank, self.mpi_helper.MPI.SUM, root=0
        )
      removed_experiments_filter = self.mpi_helper.comm.reduce(
        removed_experiments_filter_rank, self.mpi_helper.MPI.SUM, root=0
        )
      self.logger.log(f"Reflections rejected because of {filter_type}: {removed_reflections_filter_rank}")
      self.logger.log(f"Experiments rejected because of {filter_type}: {removed_experiments_filter_rank}")
      if self.mpi_helper.rank == 0:
        self.logger.main_log(f"Total reflections rejected because of {filter_type}: {removed_reflections_filter}")
        self.logger.main_log(f"Total experiments rejected because of {filter_type}: {removed_experiments_filter}")

    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(experiments, reflections)
    return experiments, reflections

  def apply_significance_filter(self, experiments, reflections):
    self.logger.log_step_time("SIGNIFICANCE_FILTER")

    # Apply an I/sigma filter ... accept resolution bins only if they
    #   have significant signal; tends to screen out higher resolution observations
    #   if the integration model doesn't quite fit
    unit_cell = self.params.scaling.unit_cell
    if unit_cell is None:
      try:
        unit_cell = self.params.statistics.average_unit_cell
      except AttributeError:
        pass
    target_symm = symmetry(unit_cell = unit_cell, space_group_info = self.params.scaling.space_group)

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    kap = 'kapton_absorption_correction' in reflections
    for expt_id, experiment in enumerate(experiments):
      exp_reflections = reflections.select(reflections['id'] == expt_id)
      if not len(exp_reflections): continue

      N_obs_pre_filter = exp_reflections.size()

      N_bins_small_set = N_obs_pre_filter // self.params.select.significance_filter.min_ct
      N_bins_large_set = N_obs_pre_filter // self.params.select.significance_filter.max_ct

      # Ensure there is at least one bin.
      N_bins = max([min([self.params.select.significance_filter.n_bins,N_bins_small_set]), N_bins_large_set, 1])
      if kap:
        unattenuated = exp_reflections['kapton_absorption_correction'] == 1.0
        iterable = [exp_reflections.select(unattenuated), exp_reflections.select(~unattenuated)]
        exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index'], True)
        exp_observations = miller.array(exp_miller_indices, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))
        binner = exp_observations.setup_binner(n_bins = N_bins)
        N_bins = None
      else:
        iterable = [exp_reflections]

      #print ("\nN_obs_pre_filter %d"%N_obs_pre_filter)
      #print >> out, "Total obs %d Choose n bins = %d"%(N_obs_pre_filter,N_bins)
      #if indices_to_edge is not None:
      #  print >> out, "Total preds %d to edge of detector"%indices_to_edge.size()

      new_exp_reflections = flex.reflection_table()
      for refls in iterable:
        # Build a miller array for the experiment reflections
        exp_miller_indices = miller.set(target_symm, refls['miller_index'], True)
        exp_observations = miller.array(exp_miller_indices, refls['intensity.sum.value'], flex.sqrt(refls['intensity.sum.variance']))
        if kap:
          exp_observations.use_binning(binner)

        assert exp_observations.size() == refls.size()

        out = StringIO()
        ### !!! CRITICAL BOTTLENECK !!! ###
        bin_results = show_observations(exp_observations, out=out, n_bins=N_bins)

        if self.params.output.log_level == 0:
          self.logger.log(out.getvalue())

        acceptable_resolution_bins = [bin.mean_I_sigI > self.params.select.significance_filter.sigma for bin in bin_results]

        acceptable_nested_bin_sequences = [i for i in range(len(acceptable_resolution_bins)) if False not in acceptable_resolution_bins[:i+1]]

        if len(acceptable_nested_bin_sequences) == 0:
          continue
        else:
          N_acceptable_bins = max(acceptable_nested_bin_sequences) + 1

          imposed_res_filter = float(bin_results[N_acceptable_bins-1].d_range.split()[2])
          if self.params.output.log_level == 0:
            ident = experiment.identifier
            self.logger.log(
              "Experiment id %d, resolution cutoff %f, experiment identifier %s\n"
              %(expt_id, imposed_res_filter, ident)
            )
          else:
            self.logger.log(
              "Experiment id %d, resolution cutoff %f\n"
              %(expt_id, imposed_res_filter)
            )

          if self.params.select.significance_filter.d_min and imposed_res_filter > self.params.select.significance_filter.d_min:
            self.logger.log("Resolution below %f, rejecting"%self.params.select.significance_filter.d_min)
            continue
          imposed_res_sel = exp_observations.resolution_filter_selection(d_min=imposed_res_filter)

          assert imposed_res_sel.size() == refls.size()

          new_exp_reflections.extend(refls.select(imposed_res_sel))

      if new_exp_reflections.size() > 0:
        new_experiments.append(experiment)
        new_reflections.extend(new_exp_reflections)

      #self.logger.log("N acceptable bins %d"%N_acceptable_bins)
      #self.logger.log("Old n_obs: %d, new n_obs: %d"%(N_obs_pre_filter, exp_observations.size()))
      #if indices_to_edge is not None:
      #  print >> out, "Total preds %d to edge of detector"%indices_to_edge.size()

    self.logger.log_step_time("SIGNIFICANCE_FILTER", True)
    new_reflections.reset_ids()
    return new_experiments, new_reflections

  def _common_initial(self, experiments, reflections):
    correct_info = 'miller_index' in reflections.keys()
    if correct_info == False:
      reflections['miller_index'] = reflections['miller_index_asymmetric']
    resolution = reflections.compute_d(experiments)
    if correct_info == False:
      del reflections['miller_index']
    reflections['q2'] = 1 / resolution**2
    q2_rank = reflections['q2'].as_numpy_array()
    intensity_rank = reflections['intensity.sum.value'].as_numpy_array()

    # get bin edges in q2
    n_bins = self.params.select.reflection_filter.n_bins
    q2_min = self.mpi_helper.comm.reduce(q2_rank.min(), op=self.mpi_helper.MPI.MIN, root=0)
    q2_max = self.mpi_helper.comm.reduce(q2_rank.max(), op=self.mpi_helper.MPI.MIN, root=0)
    if self.mpi_helper.rank == 0:
      q2_bins = np.linspace(q2_min, q2_max, n_bins + 1)
    else:
      q2_bins = np.zeros(n_bins + 1)
    self.mpi_helper.comm.Bcast(q2_bins, root=0)

    # Get the mean intensity binned in q2 using a histogram method
    intensity_summation_rank, _ = np.histogram(q2_rank, bins=q2_bins, weights=intensity_rank)
    intensity_counts_rank, _ = np.histogram(q2_rank, bins=q2_bins)
    intensity_summation = np.zeros(n_bins, dtype=intensity_summation_rank.dtype)
    intensity_counts = np.zeros(n_bins, dtype=intensity_counts_rank.dtype)
    self.mpi_helper.comm.Reduce(
      intensity_summation_rank, intensity_summation, op=self.mpi_helper.MPI.SUM, root=0
      )
    self.mpi_helper.comm.Reduce(
      intensity_counts_rank, intensity_counts, op=self.mpi_helper.MPI.SUM, root=0
      )
    if self.mpi_helper.rank == 0:
      binned_mean = intensity_summation / intensity_counts
    else:
      binned_mean = np.zeros(n_bins)
    self.mpi_helper.comm.Bcast(binned_mean, root=0)

    # Normalize intensities on each rank
    intensity_normalized_rank = np.zeros(len(reflections))
    for bin_index in range(n_bins):
      indices = np.logical_and(
        q2_rank >= q2_bins[bin_index],
        q2_rank < q2_bins[bin_index + 1]
        )
      intensity_normalized_rank[indices] = intensity_rank[indices] / binned_mean[bin_index]
    reflections['intensity_normalized'] = flex.double(intensity_normalized_rank)

    # Find the lower and upper thresholds for the data's tails
    percentile = self.params.select.reflection_filter.tail_percentile
    # These flags are used to identify which reflections are in the tails. They evently will be
    # added to the reflection table to simplify logistical management
    lower_tail_flag_rank = np.zeros(len(reflections), dtype=bool)
    upper_tail_flag_rank = np.zeros(len(reflections), dtype=bool)
    lower_tail = []
    upper_tail = []
    for bin_index in range(n_bins):
      indices = np.logical_and(
        q2_rank >= q2_bins[bin_index],
        q2_rank < q2_bins[bin_index + 1]
        )
      intensity_normalized_rank_bin = intensity_normalized_rank[indices]
      q2_rank_bin = q2_rank[indices]

      bin_sizes = self.mpi_helper.comm.gather(q2_rank_bin.size, root=0)
      if self.mpi_helper.rank == 0:
        q2_bin = np.zeros(intensity_counts[bin_index])
        intensity_normalized_bin = np.zeros(intensity_counts[bin_index])
      else:
        q2_bin = None
        intensity_normalized_bin = None
      self.mpi_helper.comm.Gatherv(
        sendbuf=q2_rank_bin,
        recvbuf=[q2_bin, bin_sizes],
        root=0
        )
      self.mpi_helper.comm.Gatherv(
        sendbuf=intensity_normalized_rank_bin,
        recvbuf=[intensity_normalized_bin, bin_sizes],
        root=0
        )

      if self.mpi_helper.rank == 0:
        sort_indices = np.argsort(intensity_normalized_bin)
        intensity_normalized_bin = intensity_normalized_bin[sort_indices]
        q2_bin = q2_bin[sort_indices]

        lower = int(percentile * intensity_normalized_bin.size)
        upper = int((1 - percentile) * intensity_normalized_bin.size)

        thresholds = np.array([intensity_normalized_bin[lower], intensity_normalized_bin[upper]])
        lower_tail.append(np.column_stack((
          intensity_normalized_bin[:lower],
          q2_bin[:lower]
          )))
        upper_tail.append(np.column_stack((
          intensity_normalized_bin[upper:],
          q2_bin[upper:]
          )))
      else:
        thresholds = np.zeros(2)
      self.mpi_helper.comm.Bcast(thresholds, root=0)

      lower_tail_flag_rank[indices] = intensity_normalized_rank_bin <= thresholds[0]
      upper_tail_flag_rank[indices] = intensity_normalized_rank_bin >= thresholds[1]
    reflections['lower_tail_flag'] = flex.bool(lower_tail_flag_rank)
    reflections['upper_tail_flag'] = flex.bool(upper_tail_flag_rank)
    if self.mpi_helper.rank == 0:
      lower_tail = np.concatenate(lower_tail, axis=0)
      upper_tail = np.concatenate(upper_tail, axis=0)
    else:
      lower_tail = None
      upper_tail = None

    return reflections, upper_tail, lower_tail

  def do_diagnostics(self, reflections, model_upper, upper_tail, model_lower, lower_tail):
    def plot_outliers(I_normalized, q2, Y, tag, model):
      inlier_indices = Y == 1
      outlier_indices = Y == -1
      fig, axes = plt.subplots(1, 1, figsize=(8, 3), sharex=True)
      axes.scatter(
        q2[inlier_indices], I_normalized[inlier_indices],
        s=1, color=[0, 0, 0], marker='.', alpha=0.5, label='Inliers'
        )
      axes.scatter(
        q2[outlier_indices], I_normalized[outlier_indices],
        s=20, color=[0.8, 0, 0], marker='.', alpha=1, label='Outliers'
        )
      xx, yy = np.meshgrid(
        np.linspace(I_normalized.min(), I_normalized.max(), 150),
        np.linspace(q2.min(), q2.max(), 150)
        )
      Z = model.predict(np.c_[xx.ravel(), yy.ravel()])
      Z = Z.reshape(xx.shape)
      # https://github.com/matplotlib/matplotlib/issues/23303
      # The label for the contour does not appear in the legend. Must add manually.
      contour = axes.contour(yy, xx, Z, levels=[0], linewidths=2, colors="green")
      contour_handle, _ = contour.legend_elements()
      handles, labels = axes.get_legend_handles_labels()
      handles += contour_handle
      labels += ['Decision Boundary']
      axes.set_xlabel(r'Resolution ($\mathrm{\AA}$)')
      xticks = axes.get_xticks()
      xticks = xticks[xticks > 0]
      xticklabels = [f'{l:0.2f}' for l in 1 / np.sqrt(xticks)]
      axes.set_xticks(xticks)
      axes.set_xticklabels(xticklabels)
      axes.set_ylabel('Normalized Intensity')
      if tag == 'upper':
        loc = 'upper right'
      elif tag == 'lower':
        loc = 'lower right'
      axes.legend(handles, labels, loc=loc, frameon=False)
      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + f'_model_pred_{tag}.png'
        ))
      plt.close()

    intensity_normalized = self.mpi_helper.comm.gather(
      reflections['intensity_normalized'].as_numpy_array(), root=0
      )
    q2 = self.mpi_helper.comm.gather(reflections['q2'].as_numpy_array(), root=0)
    if self.mpi_helper.rank == 0:
      import matplotlib.pyplot as plt
      import os
      plot_outliers(upper_tail[:, 0], upper_tail[:, 1], model_upper.predict(upper_tail), 'upper', model_upper)
      plot_outliers(lower_tail[:, 0], lower_tail[:, 1], model_lower.predict(lower_tail), 'lower', model_lower)

      fig, axes = plt.subplots(3, 1, figsize=(8, 6), sharex=True)
      axes[0].scatter(
        np.concatenate(q2), np.concatenate(intensity_normalized),
        s=1, color=[0, 0, 0], marker='.', alpha=0.5
        )
      axes[1].scatter(
        upper_tail[:, 1], upper_tail[:, 0],
        s=1, color=[0, 0, 0], marker='.', alpha=0.5
        )
      axes[2].scatter(
        lower_tail[:, 1], lower_tail[:, 0],
        s=1, color=[0, 0, 0], marker='.', alpha=0.5
        )
      axes[2].set_xlabel(r'$q^2$ = 1/$d^2$ (1/$\mathrm{\AA^2}$)')
      for i in range(3):
        axes[i].set_ylabel('Normalized Intensity')
      fig.tight_layout()
      fig.savefig(os.path.join(
        self.params.output.output_dir,
        self.params.output.prefix + '_normalized_I.png'
        ))
      plt.close()

  def apply_isolation_forest(self, experiments, reflections):
    self.logger.log_step_time("ISOLATION_FOREST")
    from sklearn.ensemble import IsolationForest

    reflections, upper_tail, lower_tail = self._common_initial(experiments, reflections)
    if self.mpi_helper.rank == 0:
      sampling_fraction = self.params.select.reflection_filter.sampling_fraction
      model_lower = IsolationForest(
        n_estimators=self.params.select.reflection_filter.n_estimators,
        contamination=self.params.select.reflection_filter.contamination_lower,
        max_features=2,
        max_samples=int(sampling_fraction*lower_tail.shape[0]),
        random_state=self.params.select.reflection_filter.random_seed
        )
      model_lower.fit(lower_tail)
      model_upper = IsolationForest(
        n_estimators=self.params.select.reflection_filter.n_estimators,
        contamination=self.params.select.reflection_filter.contamination_upper,
        max_features=2,
        max_samples=int(sampling_fraction*upper_tail.shape[0]),
        random_state=self.params.select.reflection_filter.random_seed
        )
      model_upper.fit(upper_tail)
    else:
      model_lower = None
      model_upper = None

    if self.params.select.reflection_filter.do_diagnostics:
      self.do_diagnostics(reflections, model_upper, upper_tail, model_lower, lower_tail)
    new_experiments, new_reflections = self._common_final(
      experiments, reflections, model_lower, model_upper, 'isolation forest'
      )
    self.logger.log_step_time("ISOLATION_FOREST", True)
    return new_experiments, new_reflections

  def _common_final(self, experiments, reflections, model_lower, model_upper, filter_type):
    model_lower = self.mpi_helper.comm.bcast(model_lower, root=0)
    model_upper = self.mpi_helper.comm.bcast(model_upper, root=0)

    lower_tail_indices = reflections['lower_tail_flag'].as_numpy_array()
    upper_tail_indices = reflections['upper_tail_flag'].as_numpy_array()
    lower_tail_reflections = reflections.select(reflections['lower_tail_flag'])
    upper_tail_reflections = reflections.select(reflections['upper_tail_flag'])
    lower_outliers = model_lower.predict(np.column_stack((
      lower_tail_reflections['intensity_normalized'].as_numpy_array(),
      lower_tail_reflections['q2'].as_numpy_array()
      )))
    upper_outliers = model_upper.predict(np.column_stack((
      upper_tail_reflections['intensity_normalized'].as_numpy_array(),
      upper_tail_reflections['q2'].as_numpy_array()
      )))

    inlier = np.ones(len(reflections), dtype=bool)
    if self.params.select.reflection_filter.apply_lower:
      if np.sum(lower_outliers == -1) > 0:
        indices = np.argwhere(lower_tail_indices)[:, 0]
        inlier[indices[lower_outliers == -1]] = False
    if self.params.select.reflection_filter.apply_upper:
      if np.sum(upper_outliers == -1) > 0:
        indices = np.argwhere(upper_tail_indices)[:, 0]
        inlier[indices[upper_outliers == -1]] = False
    reflections['inlier'] = flex.bool(inlier)

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()
    for expt_id, experiment in enumerate(experiments):
      exp_reflections = reflections.select(reflections['id'] == expt_id)
      if not len(exp_reflections): continue

      new_exp_reflections = flex.reflection_table()
      new_exp_reflections.extend(exp_reflections.select(exp_reflections['inlier']))
      if new_exp_reflections.size() > 0:
        new_experiments.append(experiment)
        new_reflections.extend(new_exp_reflections)

    new_reflections.reset_ids()
    del new_reflections['q2']
    del new_reflections['intensity_normalized']
    del new_reflections['lower_tail_flag']
    del new_reflections['upper_tail_flag']
    del new_reflections['inlier']
    return new_experiments, new_reflections


if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reflection_filter)

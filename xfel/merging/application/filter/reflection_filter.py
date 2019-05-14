from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from rstbx.dials_core.integration_core import show_observations
from cctbx import miller
from cctbx.crystal import symmetry
from six.moves import cStringIO as StringIO

class reflection_filter(worker):
  '''Reject individual reflections based on various criteria'''

  def __repr__(self):
    return 'Filter reflections'

  def run(self, experiments, reflections):
    if 'significance_filter' in self.params.select.algorithm:
      experiments, reflections = self.apply_significance_filter(experiments, reflections)
    return experiments, reflections

  def apply_significance_filter(self, experiments, reflections):

    self.logger.log_step_time("SIGNIFICANCE_FILTER")

    # Apply an I/sigma filter ... accept resolution bins only if they
    #   have significant signal; tends to screen out higher resolution observations
    #   if the integration model doesn't quite fit
    target_symm = symmetry(unit_cell = self.params.scaling.unit_cell, space_group_info = self.params.scaling.space_group)

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    for experiment in experiments:
      exp_reflections = reflections.select(reflections['exp_id'] == experiment.identifier)

      N_obs_pre_filter = exp_reflections.size()

      N_bins_small_set = N_obs_pre_filter // self.params.select.significance_filter.min_ct
      N_bins_large_set = N_obs_pre_filter // self.params.select.significance_filter.max_ct

      # Ensure there is at least one bin.
      N_bins = max([min([self.params.select.significance_filter.n_bins,N_bins_small_set]), N_bins_large_set, 1])

      #print ("\nN_obs_pre_filter %d"%N_obs_pre_filter)
      #print >> out, "Total obs %d Choose n bins = %d"%(N_obs_pre_filter,N_bins)
      #if indices_to_edge is not None:
      #  print >> out, "Total preds %d to edge of detector"%indices_to_edge.size()

      # Build a miller array for the experiment reflections
      exp_miller_indices = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      exp_observations = miller.array(exp_miller_indices, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))

      assert exp_observations.size() == exp_reflections.size()

      out = StringIO()
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

        imposed_res_sel = exp_observations.resolution_filter_selection(d_min=imposed_res_filter)

        assert imposed_res_sel.size() == exp_reflections.size()

        new_exp_reflections = exp_reflections.select(imposed_res_sel)

        if new_exp_reflections.size() > 0:
          new_experiments.append(experiment)
          new_reflections.extend(new_exp_reflections)

        #self.logger.log("N acceptable bins %d"%N_acceptable_bins)
        #self.logger.log("Old n_obs: %d, new n_obs: %d"%(N_obs_pre_filter, exp_observations.size()))
        #if indices_to_edge is not None:
        #  print >> out, "Total preds %d to edge of detector"%indices_to_edge.size()

    removed_reflections = len(reflections) - len(new_reflections)
    removed_experiments = len(experiments) - len(new_experiments)

    self.logger.log("Reflections rejected because of significance filter: %d"%removed_reflections)
    self.logger.log("Experiments rejected because of significance filter: %d"%removed_experiments)

    # MPI-reduce total counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    total_removed_reflections  = comm.reduce(removed_reflections, MPI.SUM, 0)
    total_removed_experiments  = comm.reduce(removed_experiments, MPI.SUM, 0)

    # rank 0: log total counts
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Total reflections rejected because of significance filter: %d"%total_removed_reflections)
      self.logger.main_log("Total experiments rejected because of significance filter: %d"%total_removed_experiments)

    self.logger.log_step_time("SIGNIFICANCE_FILTER", True)

    return new_experiments, new_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reflection_filter)

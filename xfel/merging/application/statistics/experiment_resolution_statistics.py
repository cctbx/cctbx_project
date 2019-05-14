from __future__ import absolute_import, division, print_function
from six.moves import range
from dials.array_family import flex
from libtbx import table_utils
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils

class experiment_resolution_statistics(worker):
  '''Calculates experiments accepted vs resolution bins'''

  def __repr__(self):
    return 'Lattices resolution'

  def run(self, experiments, reflections):

    self.logger.log_step_time("EXPERIMENT_RESOLUTION_STATS")

    # Get pre-created resolution binning objects from the parameters
    self.resolution_binner = self.params.statistics.resolution_binner
    self.hkl_resolution_bins = self.params.statistics.hkl_resolution_bins

    # How many bins do we have?
    self.n_bins = self.resolution_binner.n_bins_all() # (self.params.statistics.n_bins + 2), 2 - to account for the hkls outside of the binner resolution range

    # To enable MPI all-rank reduction, every rank must initialize statistics array(s), even if the rank doesn't have any reflections.
    self.experiment_count_per_resolution_bins = flex.int(self.n_bins, 0)

    # Calculate, format and output statistics for each rank
    if reflections.size() > 0:
      self.count_experiments_per_resolution_bins(reflections)
      Experiment_Table_text = self.get_formatted_table(self.experiment_count_per_resolution_bins, len(experiments))
      self.logger.log(Experiment_Table_text)

    # Accumulate statistics from all ranks
    all_ranks_experiment_count_per_resolution_bins = self.mpi_helper.cumulative_flex(self.experiment_count_per_resolution_bins, flex.int)
    all_ranks_total_experiment_count = self.mpi_helper.sum(len(experiments))

    # Format and output all-rank total statistics
    if self.mpi_helper.rank == 0:
      Experiment_Table_text = self.get_formatted_table(all_ranks_experiment_count_per_resolution_bins, all_ranks_total_experiment_count)
      self.logger.main_log(Experiment_Table_text)

    self.logger.log_step_time("EXPERIMENT_RESOLUTION_STATS", True)

    return experiments, reflections

  def get_formatted_table(self, experiment_count_per_bin, total_experiment_count):
    '''Produce a table with experiment count over resolution bins'''

    table_data = [["Bin", "Resolution Range", "Lattices", "Accepted (%)"]]

    for i_bin in self.resolution_binner.range_used():
      col_legend = '%-13s' % self.resolution_binner.bin_legend(
                                                               i_bin=i_bin,
                                                               show_bin_number=False,
                                                               show_bin_range=False,
                                                               show_d_range=True,
                                                               show_counts=False)
      exp_count_abs = '%8d' % experiment_count_per_bin[i_bin]
      exp_count_percent = '%5.2f'% (100. * experiment_count_per_bin[i_bin] / total_experiment_count)
      table_data.append(['%3d' % i_bin, col_legend, exp_count_abs, exp_count_percent])

    table_data.append([""] * len(table_data[0]))
    table_data.append(["All", "", '%8d' % total_experiment_count])

    return "\n          Image Statistics\n" + table_utils.format(table_data, has_header=1, justify='center', delim=' ')

  def count_experiments_per_resolution_bins(self, reflections):
    '''For each resolution bin, count experiments that contributed reflections to that bin'''

    # Sort all reflections on asu hkls
    self.logger.log_step_time("SORT")
    self.logger.log("Sorting reflection table...")
    reflections.sort('miller_index_asymmetric')
    self.logger.log_step_time("SORT", True)

    # Initialize a dictionary to store unique experiment ids in resolution bins
    experiments_per_resolution_bins = {}
    for i_bin in range(self.n_bins):
      experiments_per_resolution_bins[i_bin] = set()

    # Accumulate experiment ids in the resolution bins where those experiments contributed reflections
    for refls in reflection_table_utils.get_next_hkl_reflection_table(reflections=reflections):
      assert refls.size() > 0
      hkl = refls[0]['miller_index_asymmetric']
      if hkl in self.hkl_resolution_bins:
        i_bin = self.hkl_resolution_bins[hkl]
        for refl in refls:
          experiments_per_resolution_bins[i_bin].add(refl['exp_id'])

    # For each bin, reduce the sets of unique experiment ids to their count
    for i_bin in range(self.resolution_binner.n_bins_all()):
      self.experiment_count_per_resolution_bins[i_bin] = len(experiments_per_resolution_bins[i_bin])

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(experiment_resolution_statistics)

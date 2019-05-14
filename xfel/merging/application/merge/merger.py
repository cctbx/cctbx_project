from __future__ import absolute_import, division, print_function
from dials.array_family import flex
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from dxtbx.model.experiment_list import ExperimentList

try:
  import resource
  import platform
  def get_memory_usage():
    # getrusage returns kb on linux, bytes on mac
    units_per_mb = 1024
    if platform.system() == "Darwin":
      units_per_mb = 1024*1024
    return ('Memory usage: %.1f MB' % (int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb))
except ImportError:
  def debug_memory_usage():
    pass

class merger(worker):
  """
  Merges multiple measurements of symmetry-reduced HKLs.
  """

  def __repr__(self):
    return "Merge multiple measurements of symmetry-reduced HKLs"

  def merging_reflection_table(self):
    '''Create a reflection table for storing merged HKLs'''
    table = flex.reflection_table()
    table['miller_index'] = flex.miller_index()
    table['intensity']    = flex.double()
    table['esd']          = flex.double()
    table['rmsd']         = flex.double()
    table['multiplicity'] = flex.int()
    return table

  def calc_reflection_intensity_stats(self, reflections):
    '''Calculate intensity statistics for reflection table'''
    multiplicity = len(reflections)
    assert multiplicity != 0

    stats = flex.mean_and_variance(reflections['intensity.sum.value'])
    propagated_esd = (flex.sum(reflections['intensity.sum.variance']) ** 0.5)/ multiplicity

    rmsd = 0.0
    if multiplicity > 1:
      rmsd = stats.unweighted_sample_standard_deviation()

    return {'intensity'     : stats.mean(),
            'esd'           : propagated_esd,
            'rmsd'          : rmsd,
            'multiplicity'  : multiplicity}

  def run(self, experiments, reflections):

    # merge reflection intensities: calculate the average and other statistics
    self.logger.log_step_time("AVERAGE")
    self.logger.log("Averaging intensities...")
    all_rank_merged_reflections = self.merging_reflection_table()

    if len(reflections) > 0:
      for hkl_reflection_table in reflection_table_utils.get_next_hkl_reflection_table(reflections):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=hkl_reflection_table)
        intensity_stats['miller_index'] = hkl_reflection_table[0].get('miller_index_asymmetric')
        all_rank_merged_reflections.append(intensity_stats)

    self.logger.log("Merged intensities for %d HKLs"%(all_rank_merged_reflections.size()))
    self.logger.log_step_time("AVERAGE", True)

    # gather all merged intensities at rank 0
    self.logger.log_step_time("GATHER")
    if self.mpi_helper.rank != 0:
      self.logger.log("Executing MPI gathering of all reflection tables at rank 0...")
    all_merged_reflection_tables = self.mpi_helper.comm.gather(all_rank_merged_reflections, root = 0)
    all_experiment_lists = self.mpi_helper.comm.gather(experiments, root = 0)
    #results = self.mpi_helper.comm.gather((experiments, all_rank_merged_reflections), root = 0)
    self.logger.log_step_time("GATHER", True)

    # rank 0: concatenate all merged intensities into the final table
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("MERGE")
      all_experiments = ExperimentList()
      final_merged_reflection_table = self.merging_reflection_table()
      self.logger.log("Performing final merging of reflection tables and experiments received from all ranks...")
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      for expts in all_experiment_lists:
        all_experiments.extend(expts)
      self.logger.main_log("Total merged HKLs: {}".format(final_merged_reflection_table.size()))
      self.logger.log_step_time("MERGE", True)

      return all_experiments, final_merged_reflection_table
    else:
      return None, None

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(merge)

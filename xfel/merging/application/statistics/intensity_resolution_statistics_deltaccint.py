from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex
import math
from libtbx import adopt_init_args
from libtbx.str_utils import format_value
from libtbx import table_utils
from xfel.merging.application.reflection_table_utils import reflection_table_utils
from six.moves import cStringIO as StringIO
from cctbx.crystal import symmetry
from cctbx import miller
from scitbx.math import basic_statistics
from xfel.merging.application.statistics.intensity_resolution_statistics import intensity_resolution_statistics
import copy

class intensity_resolution_statistics_deltaccint(intensity_resolution_statistics):
  '''Calculates hkl intensity statistics for resolution bins'''

  def barrier(self):
    print(f'rank {self.mpi_helper.rank} at barrier')
    self.mpi_helper.comm.barrier()
  def select_fast_(self, refl, mask):
    return refl.select_fast(mask)
  def select_(self, refl, mask):
    return refl.select(mask)

  def __repr__(self):
    return 'Intensity resolution statistics deltaccint'

  def run(self, experiments, reflections):
    self.mpi_helper.comm.barrier()
    if self.mpi_helper.rank==0:
      self.log_cache = []
    self.logger.log_step_time("INTENSITY_STATISTICS_DELTACCINT")

    assert experiments == None, "Must be run after group"

    refl_copy = copy.deepcopy(reflections)
    refl_pruned = reflection_table_utils.prune_reflection_table_keys(
        reflections=refl_copy,
        keys_to_keep=['id', 'intensity.sum.value', 'intensity.sum.variance', 'miller_index_asymmetric'],
        keys_to_ignore=[])

    expt_ids = self.mpi_helper.comm.gather(refl_pruned.experiment_identifiers().values(), 0)
    if self.mpi_helper.rank == 0:
      all_expt_ids = flex.std_string()
      for expt_ids_ in expt_ids:
        all_expt_ids.extend(expt_ids_)
      expt_ids = list(set(all_expt_ids))
    else:
      expt_ids = None
    all_expt_ids = self.mpi_helper.comm.bcast(expt_ids, 0)

    if self.mpi_helper.rank == 0:
      self.logger.main_log("N experiments to run delta cc int on: %d"%len(all_expt_ids))

    mapping = refl_pruned.experiment_identifiers().values().i_seqs_by_value()

    reflections_odd = reflection_table_utils.select_odd_experiment_reflections(refl_pruned)
    reflections_even = reflection_table_utils.select_even_experiment_reflections(refl_pruned)
    assert len(reflections_even) + len(reflections_odd) == len(refl_pruned)

    for i, expt_id in enumerate(all_expt_ids):
      if expt_id in mapping:
        assert len(mapping[expt_id]) == 1
        expt_idx = mapping[expt_id][0]
        mask_odd = (reflections_odd['id'] != expt_idx)
        mask_even = (reflections_even['id'] != expt_idx)
        subset_odd = self.select_fast_(reflections_odd, mask_odd)
        subset_even = self.select_fast_(reflections_even, mask_even)
      else:
        subset_odd = reflections_odd
        subset_even = reflections_even
#      if self.mpi_helper.rank in [0,1]:
#        import line_profiler
#        lp = line_profiler.LineProfiler(self.run_single_delta)
#        lp.enable()
      self.run_single_delta(subset_odd, subset_even, expt_id, i)
#      if self.mpi_helper.rank in [0,1]:
#        lp.disable()
#        lp.print_stats()

    if self.mpi_helper.rank == 0:
      self.logger.main_log('\n'.join(self.log_cache))
    self.logger.log_step_time("INTENSITY_STATISTICS_DELTACCINT", True)

    self.barrier()
    return experiments, reflections

  def run_single_delta(self, subset_odd, subset_even, delta_expt_id, i_expt):
    self.last_bin_incomplete = False
    self.suggested_resolution_scalar = -1.0


#    if self.mpi_helper.rank in [0,1]:
#      import line_profiler
#      lp = line_profiler.LineProfiler(self.calculate_cc_int)
#      lp.enable()
    self.calculate_cc_int(subset_odd, subset_even)
#    if self.mpi_helper.rank in [0,1]:
#      lp.disable()
#      lp.print_stats()

    if self.mpi_helper.rank == 0:
      if i_expt%100==0:
        self.logger.main_log("%d done"%(i_expt))
      Table = self.Total_CC_OneHalf_Table
      if Table is not None:
        self.log_cache.append("%s %d/%d\t\t%f"%(delta_expt_id,
                                                 Table.cumulative_observed_matching_asu_count,
                                                 Table.cumulative_theor_asu_count,
                                                 Table.cumulative_cross_correlation))

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(intensity_resolution_statistics_deltaccint)

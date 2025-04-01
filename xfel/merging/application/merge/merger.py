from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.reflection_table_utils import \
    reflection_table_utils as rt_util
from cctbx.crystal import symmetry
from cctbx import miller
import os
from six.moves import cStringIO as StringIO
import numpy as np
from scitbx.array_family import flex
from dials.util.version import dials_version

class merger(worker):
  """
  Merges multiple measurements of symmetry-reduced HKLs.
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(merger, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Merge multiple measurements of symmetry-reduced HKLs"

  def run(self, experiments, reflections):

    sel_col = "id"
    if self.params.statistics.shuffle_ids:
      refl_ids = list(set(reflections["id"]))
      new_ids = np.random.permutation(refl_ids).astype(np.int32)
      new_id_map = {old_id: new_id for old_id, new_id in zip(refl_ids, new_ids)}
      new_id_col = [new_id_map[i] for i in reflections["id"]]
      reflections["shuffled_id"] = flex.int(new_id_col)
      sel_col = "shuffled_id"

    # select, merge and output odd reflections
    odd_reflections = rt_util.select_odd_experiment_reflections(reflections, sel_col)
    odd_reflections_merged = rt_util.merge_reflections(
        odd_reflections,
        self.params.merging.minimum_multiplicity,
        thresh=self.params.filter.outlier.mad_thresh
    )
    self.gather_and_output_reflections(odd_reflections_merged, 'odd')

    # select, merge and output even reflections
    even_reflections = rt_util.select_even_experiment_reflections(reflections, sel_col)
    even_reflections_merged = rt_util.merge_reflections(
        even_reflections,
        self.params.merging.minimum_multiplicity,
        thresh=self.params.filter.outlier.mad_thresh
    )
    self.gather_and_output_reflections(even_reflections_merged, 'even')

    # merge and output all reflections
    name = "merged_good_refls2/rank%d" % self.mpi_helper.comm.rank
    all_reflections_merged = rt_util.merge_reflections(
        reflections,
        self.params.merging.minimum_multiplicity,
        nameprefix=name,
        thresh=self.params.filter.outlier.mad_thresh
    )
    self.gather_and_output_reflections(all_reflections_merged, 'all')

    return None, reflections

  def gather_and_output_reflections(self, reflections, selection_name):
    # gather merged HKLs from all ranks
    self.logger.log_step_time("GATHER")
    self.logger.log("Running MPI-gather on merged %s HKLs..."%selection_name)
    all_merged_reflection_tables = self.mpi_helper.comm.gather(reflections, root = 0)
    self.logger.log_step_time("GATHER", True)

    # concatenate all merged HKLs
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("MERGE")
      final_merged_reflection_table = rt_util.merged_reflection_table()
      self.logger.log("Concatenating merged %s HKLs at rank 0..."%selection_name)
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      self.logger.main_log("Total (not limited by resolution) merged %s HKLs: %d"%(selection_name, final_merged_reflection_table.size()))
      self.logger.log_step_time("MERGE", True)

      # output as mtz
      if len(final_merged_reflection_table) > 0:
        self.output_reflections_mtz(final_merged_reflection_table, selection_name)

      # free the memory
      del all_merged_reflection_tables
      del final_merged_reflection_table

  def output_reflections_mtz(self, reflections, filename_postfix):
    if self.params.merging.set_average_unit_cell:
      assert 'average_unit_cell' in (self.params.statistics).__dict__
      unit_cell = self.params.statistics.__phil_get__('average_unit_cell')
    else:
      unit_cell = self.params.scaling.unit_cell

    final_symm = symmetry(
                          unit_cell=unit_cell,
                          space_group_info = self.params.scaling.space_group)

    assert 'average_wavelength' in (self.params.statistics).__dict__
    wavelength = self.params.statistics.__phil_get__('average_wavelength')

    if self.params.merging.d_max is None:
      self.logger.main_log("Output merged HKLs limited by %f A resolution"%(self.params.merging.d_min))
    else:
      self.logger.main_log("Output merged HKLs limited by (%f - %f) A resolution range"%(self.params.merging.d_max, self.params.merging.d_min))

    all_obs = miller.array(
      miller_set=miller.set(final_symm, reflections['miller_index'], not self.params.merging.merge_anomalous),
      data=reflections['intensity'],
      sigmas=reflections['sigma']).resolution_filter(
                                   d_min=self.params.merging.d_min,
                                   d_max=self.params.merging.d_max).set_observation_type_xray_intensity()

    mtz_file = os.path.join(self.params.output.output_dir, "%s_%s.mtz"%(self.params.output.prefix, filename_postfix))

    mtz_out = all_obs.as_mtz_dataset(
      column_root_label="Iobs",
      title=self.params.output.title,
      wavelength=wavelength)
    mtz_out.add_miller_array(
      miller_array=all_obs.average_bijvoet_mates(),
      column_root_label="IMEAN")

    if self.params.merging.include_multiplicity_column:
      all_mult = miller.array(
        miller_set=miller.set(final_symm, reflections['miller_index'], not self.params.merging.merge_anomalous),
        data=reflections['multiplicity']).resolution_filter(
                                          d_min=self.params.merging.d_min,
                                          d_max=self.params.merging.d_max).set_observation_type_xray_intensity()
      mtz_out.add_miller_array(
        miller_array=all_mult,
        column_root_label='multiplicity')

    mtz_obj = mtz_out.mtz_object()
    mtz_obj.write(mtz_file)
    self.logger.main_log("Output anomalous and mean data:\n    %s" %os.path.abspath(mtz_file))
    self.logger.main_log("Output data summary:")
    out = StringIO()
    all_obs.show_summary(prefix="  ", f=out)
    self.logger.main_log(out.getvalue())

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(merge)

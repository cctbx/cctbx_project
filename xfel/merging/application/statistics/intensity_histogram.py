from __future__ import division
from xfel.merging.application.worker import worker
from dials.array_family import flex
from six.moves import cStringIO as StringIO
from xfel.merging.application.reflection_table_utils import \
    reflection_table_utils as rt_util
from cctbx.crystal import symmetry
from cctbx import miller


class intensity_histogram(worker):
  '''Builds a histogram of reflection intensities'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(intensity_histogram, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Intensity histogram'

  def run(self, experiments, reflections):
    if self.params.statistics.histogram.mode == "legacy":
      return self.run_legacy(experiments, reflections)
    elif self.params.statistics.histogram.mode == "individual":
      return self.run_individual(experiments, reflections)

  def run_legacy(self, experiments, reflections):
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
#    if self.mpi_helper.rank==0:
#      import IPython;IPython.embed()
#    else:
#      import time;time.sleep(10000)
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("INTENSITY_HISTOGRAM")
      self.histogram(reflections['intensity.sum.value'])
      self.logger.log_step_time("INTENSITY_HISTOGRAM", True)

    return experiments, reflections

  def run_individual(self, experiments, reflections):
    all_reflections_merged = rt_util.merge_reflections(
        reflections,
        self.params.merging.minimum_multiplicity
    )
    all_merged_reflection_tables = self.mpi_helper.comm.gather(
        all_reflections_merged, root=0
    )
    if self.mpi_helper.rank == 0:
      final_merged_reflection_table = rt_util.merged_reflection_table()
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)

      unit_cell = self.params.scaling.unit_cell
      final_symm = symmetry(
          unit_cell=unit_cell,
          space_group_info=self.params.scaling.space_group
      )
      all_obs = miller.array(
          miller_set=miller.set(final_symm, final_merged_reflection_table['miller_index'], False),
          data=final_merged_reflection_table['intensity'],
          sigmas=final_merged_reflection_table['sigma']
      ).resolution_filter(
          d_min=self.params.merging.d_min,
          d_max=self.params.merging.d_max
      ).set_observation_type_xray_intensity()
      print(list(all_obs[0:2].indices()))
      print(list(all_obs[0:2].data()))

#    if self.mpi_helper.rank==0:
#      import IPython;IPython.embed()
#    else:
#      import time;time.sleep(10000)
    return experiments, reflections

    


  def histogram(self, data):
    from matplotlib import pyplot as plt
    nslots = 100
    histogram = flex.histogram(
                               data=data,
                               n_slots=nslots)
    out = StringIO()
    histogram.show(f=out, prefix="    ", format_cutoffs="%6.2f")
    self.logger.main_log(out.getvalue() + '\n' + "Total: %d"%data.size() + '\n')

    if False:
      fig = plt.figure()
      plt.bar(histogram.slot_centers(), histogram.slots(), align="center", width=histogram.slot_width())
      plt.show()

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(intensity_histogram)

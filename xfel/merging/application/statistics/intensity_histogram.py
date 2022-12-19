from __future__ import division
from xfel.merging.application.worker import worker
from dials.array_family import flex
from six.moves import cStringIO as StringIO
from xfel.merging.application.reflection_table_utils import \
    reflection_table_utils as rt_util
from cctbx.crystal import symmetry
from cctbx import miller
import matplotlib.pyplot as plt
import numpy as np


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


      d_vals = [x[1] for x in all_obs.d_spacings()] # dspacings are ((h,k,l),d)
      n_unique = len(d_vals)
      d_percentiles = [0, 5, 10, 40, 80]
      to_plot = [[] for _ in d_percentiles]
      d_slice_starts = [int(x*n_unique/100) for x in d_percentiles]
      d_slice_size = int(n_unique/100)
      d_slices = [slice(x, x+d_slice_size) for x in d_slice_starts]

      # This has to be a list because we're going to slice it... vvv
      d_i_miller = list(zip(d_vals, all_obs.data(), all_obs.indices()))
      d_i_miller = sorted(d_i_miller, reverse=True)
      outer_slices = [d_i_miller[s] for s in d_slices]

      for i_outer, outer in enumerate(outer_slices):
        # now sort the d-slices by decreasing intensity
        outer_sorted = sorted(outer, key=lambda x:x[1], reverse=True)
        n_inner = len(outer_sorted)
        iobs_percentiles = [0,5,10,40,80]
        indices_inner = [int(x*n_inner/100) for x in iobs_percentiles]
        for i_inner, i in enumerate(indices_inner):
          item = outer_sorted[i]
          to_plot[i_outer].append(item)
          #self.single_histogram(reflections, m_i, ax[i_outer][i_inner])

    if self.mpi_helper.rank != 0: to_plot = None
    to_plot = self.mpi_helper.comm.bcast(to_plot, root=0)
    if self.mpi_helper.rank == 0:
      fig, axes = plt.subplots(5,5)
      fig.set_size_inches((12,8))
    else:
      fig, axes = None, None
    for i1, l in enumerate(to_plot):
      for i2, m_i in enumerate(l):
        self.single_histogram(reflections, m_i, axes, i1, i2)
    if self.mpi_helper.rank == 0:
      plt.savefig('plot.png')




#    if self.mpi_helper.rank==0:
#      import IPython;IPython.embed()
#    else:
#      import time;time.sleep(10000)
    return experiments, reflections

  def single_histogram(self, reflections, item, axes, i1, i2):
    d, iobs, miller = item
    sel = reflections['miller_index_asymmetric'] == miller
    matching_refl = reflections.select(sel)
    matching_intensities = matching_refl['intensity.sum.value']
    all_matching_intensities = self.mpi_helper.comm.gather(
        matching_intensities, root=0
    )

    if self.mpi_helper.rank == 0:
      final_intensities = flex.double()
      for x in all_matching_intensities:
        final_intensities.extend(x)
      i_arr = final_intensities.as_numpy_array()

      q1 = np.percentile(i_arr, 25)
      q3 = np.percentile(i_arr, 75)
      med = np.percentile(i_arr, 50)
      iqr = q3-q1
      xmin = min(0, q1) - 2*iqr
      xmin = max(xmin, min(i_arr))
      xmax = iobs + 2*iqr
      xmin = np.percentile(i_arr, 2.5) - iqr/5
      xmax = np.percentile(i_arr, 97.5) + iqr/5


      ax = axes[i1][i2]
      hist = ax.hist(i_arr, bins=np.linspace(xmin, xmax, 50))
      ax.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
      title = "%.2f, %.0f, %d" %(d, iobs, len(i_arr))
      ax.set_title(title)

      ymax = max(hist[0])
      ax.plot((0,0), (0,ymax), 'k-')
      ax.plot((iobs,iobs), (0,ymax), 'r-')



    


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

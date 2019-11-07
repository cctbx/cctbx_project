from __future__ import division
from xfel.merging.application.worker import worker
from dials.array_family import flex
from six.moves import cStringIO as StringIO

class intensity_histogram(worker):
  '''Builds a histogram of reflection intensities'''

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(intensity_histogram, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Intensity histogram'

  def run(self, experiments, reflections):
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("INTENSITY_HISTOGRAM")
      self.histogram(reflections['intensity'])
      self.logger.log_step_time("INTENSITY_HISTOGRAM", True)

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

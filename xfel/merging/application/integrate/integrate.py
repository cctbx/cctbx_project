from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from dxtbx.imageset import ImageSetFactory

from dials.command_line.stills_process import Processor
class integrate_only_processor(Processor):
  def __init__(self, params):
    self.params = params

class integrate(worker):
  """
  Calls the stills process version of dials.integrate
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(integrate, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Integrate reflections'

  def run(self, experiments, reflections):
    self.logger.log_step_time("INTEGRATE")

    # Re-generate the image sets using their format classes so we can read the raw data
    for expt in experiments:
      expt.imageset = ImageSetFactory.make_imageset(expt.imageset.paths(), single_file_indices=expt.imageset.indices())

    processor = integrate_only_processor(self.params)
    reflections = processor.integrate(experiments, reflections)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(integrate)

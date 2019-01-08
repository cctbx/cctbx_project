from __future__ import print_function, division
from xfel.merging.application.worker import worker

class frame_scaler(worker):
  """
  Scale frames of data
  """

  def run(self, experiments, reflections):

    self.logger.log_step_time("SCALE_FRAMES")

    self.logger.log_step_time("SCALE_FRAMES", True)

    return experiments, reflections

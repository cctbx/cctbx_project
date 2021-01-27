from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker

class lunus(worker):
  """
  Calls into the Lunus library to do diffuse scatter integration

  See DOI: 10.1007/978-1-59745-483-4_17
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(lunus, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Compute diffuse scatter using Lunus'

  def run(self, experiments, reflections):

    self.logger.log_step_time("LUNUS")

    # do work

    self.logger.log_step_time("LUNUS", True)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(lunus)

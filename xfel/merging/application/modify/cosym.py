from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage

class cosym(worker):
  """
  Resolve indexing ambiguity using dials.cosym
  """
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(cosym, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Resolve indexing ambiguity using dials.cosym'

  def run(self, experiments, reflections):

    self.logger.log_step_time("COSYM")

    self.logger.log_step_time("COSYM", True)
    self.logger.log("Memory usage: %d MB"%get_memory_usage())

    return experiments, result

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(reindex_to_reference)

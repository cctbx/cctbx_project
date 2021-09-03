from __future__ import division
from xfel.merging.application.worker import worker


class annulus_statistics(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(publisher, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return "Calculate reflection statistics in a resolution annulus"

  def run(self, experiments, reflections):
    pass

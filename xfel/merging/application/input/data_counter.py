from __future__ import absolute_import, division, print_function

from xfel.merging.application.mpi_helper import mpi_helper
from xfel.merging.application.mpi_logger import mpi_logger

class data_counter(object):
  def __init__(self, params):
    # create logger
    self.logger = mpi_logger(params)
    # create MPI helper
    self.mpi_helper = mpi_helper()

  def count(self, experiments, reflections):

    self.logger.log_step_time("CALC_LOAD_STATISTICS")

    # count experiments and reflections
    experiment_count = len(experiments)
    reflection_count = len(reflections)

    # count images
    image_count = sum(len(iset) for iset in experiments.imagesets())


    # MPI-reduce all counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    total_experiment_count  = comm.reduce(experiment_count, MPI.SUM, 0)
    total_image_count       = comm.reduce(image_count, MPI.SUM, 0)
    total_reflection_count  = comm.reduce(reflection_count, MPI.SUM, 0)
    max_reflection_count    = comm.reduce(reflection_count, MPI.MAX, 0)
    min_reflection_count    = comm.reduce(reflection_count, MPI.MIN, 0)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0:
      self.logger.main_log('All ranks have read %d experiments'%total_experiment_count)
      self.logger.main_log('All ranks have read %d images'%total_image_count)
      self.logger.main_log('All ranks have read %d reflections'%total_reflection_count)
      self.logger.main_log('The maximum number of reflections loaded per rank is: %d reflections'%max_reflection_count)
      self.logger.main_log('The minimum number of reflections loaded per rank is: %d reflections'%min_reflection_count)

    self.logger.log_step_time("CALC_LOAD_STATISTICS", True)

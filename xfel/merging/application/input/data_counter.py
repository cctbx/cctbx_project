from __future__ import division
from six.moves import range

from xfel.merging.application.mpi_helper import mpi_helper
from xfel.merging.application.mpi_logger import mpi_logger

#from libtbx.mpi4py import MPI


class data_counter(object):
  def __init__(self, params):
    # create logger
    self.logger = mpi_logger(params)
    # create MPI helper
    self.mpi_helper = mpi_helper()

  def count(self, experiments, reflections):
    # count experiments and reflections
    experiment_count = len(experiments)
    reflection_count = len(reflections)
    
    # count images
    image_count = 0
    all_imgs = []
    for iset in experiments.imagesets():
      all_imgs.extend(iset.paths())
      image_count = len(set(all_imgs))
      
    #from IPython import embed; embed()
      
    self.logger.log_step_time("LOAD_STATISTICS")
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
        
    # MPI-reduce all counts
    total_experiment_count  = comm.reduce(experiment_count, MPI.SUM, 0)
    total_image_count       = comm.reduce(image_count, MPI.SUM, 0)
    total_reflection_count  = comm.reduce(reflection_count, MPI.SUM, 0)
    max_reflection_count    = comm.reduce(reflection_count, MPI.MAX, 0)
    min_reflection_count    = comm.reduce(reflection_count, MPI.MIN, 0)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0:
      self.logger.log('All ranks have read %d experiments'%total_experiment_count)
      self.logger.log('All ranks have read %d images'%total_image_count)
      self.logger.log('All ranks have read %d reflections'%total_reflection_count)
      self.logger.log('The maximum number of reflections loaded per rank is: %d reflections'%max_reflection_count)
      self.logger.log('The minimum number of reflections loaded per rank is: %d reflections'%min_reflection_count)
    self.logger.log_step_time("LOAD_STATISTICS", True)

    


from __future__ import absolute_import, division, print_function

class data_counter(object):
  def __init__(self, params, mpi_helper=None, mpi_logger=None):

    self.mpi_helper = mpi_helper
    if self.mpi_helper == None:
      from xfel.merging.application.mpi_helper import mpi_helper
      self.mpi_helper = mpi_helper()

    self.logger = mpi_logger
    if self.logger == None:
      from xfel.merging.application.mpi_logger import mpi_logger
      self.logger = mpi_logger(params)

  def count_each(self, experiments, reflections, verbose=False):
    self.logger.log_step_time("CALC_LOAD_STATISTICS")

    # MPI-gather individual counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    # count experiments and reflections
    experiment_count = len(experiments) if experiments != None else 0
    reflection_count = len(reflections) if reflections != None else 0

    # count images
    if experiments != None:
      image_count = sum(len(iset) for iset in experiments.imagesets())
    else:
      image_count = 0

    experiment_count_list   = comm.gather(experiment_count, root=0)
    image_count_list        = comm.gather(image_count, root=0)
    reflection_count_list   = comm.gather(reflection_count, root=0)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0 and verbose:
      self.logger.main_log('Experiments by rank: ' + ', '.join(experiment_count_list))
      self.logger.main_log('Images by rank: ' + ', '.join(image_count_list))
      self.logger.main_log('Reflections by rank: ' + ', '.join(reflection_count_list))

    self.logger.log_step_time("CALC_LOAD_STATISTICS", True)

    return (experiment_count_list, image_count_list, reflection_count_list) # for load balancing purposes

  def count(self, experiments, reflections):
    self.logger.log_step_time("CALC_LOAD_STATISTICS")

    # count experiments and reflections
    experiment_count = len(experiments) if experiments != None else 0
    reflection_count = len(reflections) if reflections != None else 0

    # count images
    if experiments != None:
      image_count = sum(len(iset) for iset in experiments.imagesets())
    else:
      image_count = 0

    # MPI-reduce all counts
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    total_experiment_count  = comm.reduce(experiment_count, MPI.SUM)
    min_experiment_count    = comm.reduce(experiment_count, MPI.MIN)
    max_experiment_count    = comm.reduce(experiment_count, MPI.MAX)

    total_image_count       = comm.reduce(image_count, MPI.SUM)
    min_image_count         = comm.reduce(image_count, MPI.MIN)
    max_image_count         = comm.reduce(image_count, MPI.MAX)

    total_reflection_count  = comm.reduce(reflection_count, MPI.SUM)
    min_reflection_count    = comm.reduce(reflection_count, MPI.MIN)
    max_reflection_count    = comm.reduce(reflection_count, MPI.MAX)

    # rank 0: log data statistics
    if self.mpi_helper.rank == 0:
      self.logger.main_log('Experiments: (total,min,max): %d, %d, %d'%(total_experiment_count, min_experiment_count, max_experiment_count))
      self.logger.main_log('Images:      (total,min,max): %d, %d, %d'%(total_image_count, min_image_count, max_image_count))
      self.logger.main_log('Reflections: (total,min,max): %d, %d, %d'%(total_reflection_count, min_reflection_count, max_reflection_count))

      if total_reflection_count == 0:
        self.mpi_helper.set_error("Zero reflection count.")

    self.mpi_helper.check_errors()

    self.logger.log_step_time("CALC_LOAD_STATISTICS", True)

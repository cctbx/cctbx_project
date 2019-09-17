from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker

try:
  import resource
  import platform
  def get_memory_usage():
    # getrusage returns kb on linux, bytes on mac
    units_per_mb = 1024
    if platform.system() == "Darwin":
      units_per_mb = 1024*1024
    return ('Memory usage: %.1f MB' % (int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss) / units_per_mb))
except ImportError:
  def debug_memory_usage():
    pass

class load_balancer(worker):
  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(load_balancer, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Balance input data load'

  def divide_list_into_chunks(self, list_, n): #n: number of chunks
    return [list_[start::n] for start in range(n)]

  def run(self, experiments, reflections):
    self.logger.log("Rebalancing input load...")
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Rebalancing input load...")

    if self.params.input.parallel_file_load.balance == "global":
      new_experiments, new_reflections = self.redistribute(experiments, reflections, self.mpi_helper.comm, self.mpi_helper.size)
    elif self.params.input.parallel_file_load.balance == "per_node":
      mpi_color = int(self.mpi_helper.rank / self.params.input.parallel_file_load.ranks_per_node)
      mpi_new_rank = self.mpi_helper.rank % self.params.input.parallel_file_load.ranks_per_node
      mpi_split_comm = self.mpi_helper.comm.Split(mpi_color, mpi_new_rank)
      new_experiments, new_reflections = self.redistribute(experiments, reflections, mpi_split_comm, self.params.input.parallel_file_load.ranks_per_node)

    # Do we have any data?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(new_experiments, new_reflections)

    return new_experiments, new_reflections

  def redistribute(self, experiments, reflections, mpi_communicator, number_of_mpi_ranks):
    from dials.array_family import flex
    from dxtbx.model.experiment_list import ExperimentList

    split_experiments = self.divide_list_into_chunks(experiments, number_of_mpi_ranks)

    # if some (but not all!) chunks are empty, we want those empty chunks to be randomly distributed
    number_of_empty_chunks = [len(split_experiments[i]) for i in range(len(split_experiments))].count(0)
    if number_of_empty_chunks > 0 and len(experiments) != 0:
      import random
      #random.seed(8)
      random.shuffle(split_experiments)

    '''
    self.logger.log("Split experiment list into %d chunks"%len(split_experiments))
    for chunk in split_experiments:
      self.logger.log(len(chunk))
    '''

    split_reflections = []
    for i in range(number_of_mpi_ranks):
      split_reflections.append(flex.reflection_table())
      for experiment_id, experiment in enumerate(split_experiments[i]):
        refls = reflections.select(reflections['exp_id'] == experiment.identifier)
        split_reflections[i].extend(refls)

    self.logger.log("Split experiments and reflections")
    self.logger.log(get_memory_usage())

    del experiments
    del reflections

    self.logger.log("Deleted experiments and reflections")
    self.logger.log(get_memory_usage())

    new_split_experiments = mpi_communicator.alltoall(split_experiments)
    new_split_reflections = mpi_communicator.alltoall(split_reflections)

    self.logger.log("Consolidating data...")

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()
    for entry in new_split_experiments:
      new_experiments.extend(entry)
    for entry in new_split_reflections:
      new_reflections.extend(entry)

    self.logger.log(get_memory_usage())

    return new_experiments, new_reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(load_balancer)

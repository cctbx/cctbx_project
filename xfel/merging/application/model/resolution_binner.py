from __future__ import absolute_import, division, print_function
from xfel.merging.application.worker import worker
from xfel.merging.application.utils.memory_usage import get_memory_usage
from dials.array_family import flex

class resolution_binner(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(resolution_binner, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Set up resolution bins'

  def run(self, experiments, reflections):
    '''Set up resolution bins; assign bin number to each hkl in the full miller set'''

    '''
    # for debugging purposes on Cori
    memory_MB = get_memory_usage()
    self.logger.log("Starting resolution_binner with memory usage: %d"%memory_MB)
    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    max_memory_usage_MB = comm.reduce(memory_MB, MPI.MAX)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Maximum memory usage per process: %d MB"%max_memory_usage_MB)
    '''

    # Create resolution binner
    full_miller_set = self.params.scaling.miller_set
    full_miller_set.setup_binner(d_max=100000, d_min=self.params.merging.d_min, n_bins=self.params.statistics.n_bins)
    resolution_binner = full_miller_set.binner()

    # Save resolution binner to the parameters
    if not 'resolution_binner' in (self.params.statistics).__dict__:
      self.params.statistics.__inject__('resolution_binner', resolution_binner)
    else:
      self.params.statistics.__setattr__('resolution_binner', resolution_binner)

    # Provide resolution bin number for each asu hkl in the full miller set
    hkl_resolution_bins = {} # hkl vs resolution bin number
    hkls_with_assigned_bin = 0
    for i_bin in resolution_binner.range_used():
      bin_hkl_selection = resolution_binner.selection(i_bin)
      bin_hkls = full_miller_set.select(bin_hkl_selection)
      for hkl in bin_hkls.indices():
        assert not hkl in hkl_resolution_bins # each hkl should be assigned a bin number only once
        hkl_resolution_bins[hkl] = i_bin
        hkls_with_assigned_bin += 1
    self.logger.log("Provided resolution bin number for %d asu hkls"%(hkls_with_assigned_bin))

    from xfel import hkl_resolution_bins_cpp
    bins_cpp = hkl_resolution_bins_cpp(
          flex.miller_index(list(hkl_resolution_bins.keys())),
          flex.int(list(hkl_resolution_bins.values()))
        )
    # Save hkl bin asignments to the parameters
    if not 'hkl_resolution_bins' in (self.params.statistics).__dict__:
      self.params.statistics.__inject__('hkl_resolution_bins', hkl_resolution_bins)
    else:
      self.params.statistics.__setattr__('hkl_resolution_bins', hkl_resolution_bins)
    if not 'hkl_resolution_bins_cpp' in (self.params.statistics).__dict__:
      self.params.statistics.__inject__('hkl_resolution_bins_cpp', bins_cpp)
    else:
      self.params.statistics.__setattr__('hkl_resolution_bins_cpp', bins_cpp)

    return experiments, reflections

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(resolution_binner)

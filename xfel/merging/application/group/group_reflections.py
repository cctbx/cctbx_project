from __future__ import absolute_import, division, print_function
from six.moves import range
from xfel.merging.application.worker import worker
from dials.array_family import flex

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

class hkl_group(worker):
  '''For each asu hkl, gather all of its measurements from all ranks at a single rank, while trying to evenly distribute asu HKLs over the ranks.'''

  def __repr__(self):
    return "Group symmetry-reduced HKLs"

  def distribute_reflection_table(self):
    '''Create a reflection table for storing reflections distributed over hkl chunks'''
    table = flex.reflection_table()
    table['miller_index_asymmetric']  = flex.miller_index()
    table['intensity.sum.value']      = flex.double()
    table['intensity.sum.variance']   = flex.double()
    table['exp_id']                   = flex.std_string()
    return table

  def run(self, experiments, reflections):

    self.logger.log_step_time("GROUP")

    # set up hkl chunks to be used for all-to-all; every avialable rank participates in all-to-all, even a rank that doesn't load any data
    self.logger.log_step_time("SETUP_CHUNKS")
    self.setup_hkl_chunks()
    self.logger.log_step_time("SETUP_CHUNKS", True)

    # for the ranks, which have loaded the data, distribute the reflections over the hkl chunks
    self.logger.log_step_time("DISTRIBUTE_OVER_CHUNKS")
    self.distribute_reflections_over_hkl_chunks(reflections=reflections)
    self.logger.log_step_time("DISTRIBUTE_OVER_CHUNKS", True)

    # run all-to-all
    if self.params.parallel.a2a == 1: # 1 means: the number of slices in each chunk is 1, i.e. alltoall is done on the whole chunks
      alltoall_reflections = self.get_reflections_from_alltoall()
    else: # do alltoall on chunk slices - useful if the run-time memory is not sufficient to do alltoall on the whole chunks
      alltoall_reflections = self.get_reflections_from_alltoall_sliced(number_of_slices=self.params.parallel.a2a)

    self.logger.log_step_time("SORT")
    self.logger.log("Sorting consolidated reflection table...")
    alltoall_reflections.sort('miller_index_asymmetric')
    self.logger.log_step_time("SORT", True)

    self.logger.log_step_time("GROUP", True)

    return experiments, alltoall_reflections

  def setup_hkl_chunks(self):
    '''Set up a list of reflection tables, or chunks, for distributing reflections'''
    # split the full miller set into chunks; the number of chunks is equal to the number of ranks
    import numpy as np
    self.hkl_split_set = np.array_split(self.params.scaling.miller_set.indices(), self.mpi_helper.size)

    # initialize a list of hkl chunks - reflection tables to store distributed reflections
    self.hkl_chunks = []
    for i in range(len(self.hkl_split_set)):
      self.hkl_chunks.append(self.distribute_reflection_table())

  def distribute_reflections_over_hkl_chunks(self, reflections):
    '''Distribute reflections, according to their HKLs, over pre-set HKL chunks'''
    total_reflection_count = reflections.size()
    total_distributed_reflection_count = 0

    if total_reflection_count > 0:

      # set up two lists to be passed to the C++ extension: HKLs and chunk ids. It's basically a hash table to look up chunk ids by HKLs
      hkl_list = flex.miller_index()
      chunk_id_list = flex.int()

      for i in range(len(self.hkl_split_set)):
        for j in range(len(self.hkl_split_set[i])):
          hkl = (self.hkl_split_set[i][j][0], self.hkl_split_set[i][j][1], self.hkl_split_set[i][j][2])
          hkl_list.append(hkl)
          chunk_id_list.append(i)

      # distribute reflections over hkl chunks
      from xfel.merging import get_hkl_chunks_cpp

      get_hkl_chunks_cpp(reflections, hkl_list, chunk_id_list, self.hkl_chunks)

      for i in range(len(self.hkl_chunks)):
        total_distributed_reflection_count += len(self.hkl_chunks[i])

    self.logger.log("Distributed %d out of %d reflections"%(total_distributed_reflection_count, total_reflection_count))
    self.logger.log(get_memory_usage())

    reflections.clear()

  def get_reflections_from_alltoall(self):
    '''Use MPI alltoall method to gather all reflections with the same asu hkl from all ranks at a single rank'''
    result_reflections = self.distribute_reflection_table()

    self.logger.log_step_time("ALL-TO-ALL")
    self.logger.log("Executing MPI all-to-all...")

    received_hkl_chunks = self.mpi_helper.comm.alltoall(self.hkl_chunks)

    self.logger.log("Received %d hkl chunks after all-to-all"%len(received_hkl_chunks))

    self.logger.log_step_time("ALL-TO-ALL", True)

    self.logger.log_step_time("CONSOLIDATE")
    self.logger.log("Consolidating reflection tables...")

    for chunk in received_hkl_chunks:
      result_reflections.extend(chunk)

    self.logger.log_step_time("CONSOLIDATE", True)

    return result_reflections

  def get_reflections_from_alltoall_sliced(self, number_of_slices):
    '''Split each hkl chunk into N slices. This is needed to address the MPI alltoall memory problem'''

    result_reflections = self.distribute_reflection_table() # all reflections that the rank will receive from alltoall

    list_of_sliced_hkl_chunks = [] # if self.hkl_chunks is [A,B,C...], this list will be [[A1,A2,...,An], [B1,B2,...,Bn], [C1,C2,...,Cn], ...], where n is the number of chunk slices
    for i in range(len(self.hkl_chunks)):
      hkl_chunk_slices = []
      for chunk_slice in self.get_next_reflection_table_slice(self.hkl_chunks[i], number_of_slices):
        hkl_chunk_slices.append(chunk_slice)
      list_of_sliced_hkl_chunks.append(hkl_chunk_slices)

    self.logger.log("Ready for all-to-all...")
    self.logger.log(get_memory_usage())

    for j in range(number_of_slices):
      hkl_chunks_for_alltoall = list()
      for i in range(len(self.hkl_chunks)):
        hkl_chunks_for_alltoall.append(list_of_sliced_hkl_chunks[i][j]) # [Aj,Bj,Cj...]

      self.logger.log_step_time("ALL-TO-ALL")
      self.logger.log("Executing MPI all-to-all...")
      self.logger.log(get_memory_usage())

      received_hkl_chunks = comm.alltoall(hkl_chunks_for_alltoall)

      self.logger.log("After all-to-all received %d hkl chunks" %len(received_hkl_chunks))
      self.logger.log_step_time("ALL-TO-ALL", True)

      self.logger.log_step_time("CONSOLIDATE")
      self.logger.log("Consolidating reflection tables...")

      for chunk in received_hkl_chunks:
        result_reflections.extend(chunk)

      self.logger.log_step_time("CONSOLIDATE", True)

    return result_reflections

  @staticmethod
  def get_next_reflection_table_slice(self, reflections, n_slices):
    '''Generate an exact number of slices from a reflection table. Make slices as even as possible. If not enough reflections, generate empty tables'''
    assert n_slices >= 0

    if n_slices == 1:
      yield reflections
    else:
      import math

      generated_slices = 0
      count = len(reflections)

      if count > 0:
        # how many non-empty slices should we generate and with what stride?
        nonempty_slices = min(count, n_slices)
        stride = int(math.ceil(count / nonempty_slices))

        # generate all non-empty slices
        for i in range(0, count, stride):
          generated_slices += 1
          i2 = i + stride
          if generated_slices == nonempty_slices:
            i2 = count
          yield reflections[i:i2]

      # generate some empty slices if necessary
      empty_slices = max(0, n_slices - generated_slices)
      for i in range(empty_slices):
        yield self.distribute_reflection_table()

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(hkl_group)

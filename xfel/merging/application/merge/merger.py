from __future__ import print_function, division
from dials.array_family import flex
from xfel.merging.application.worker import worker
from six.moves import range

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

class merger(worker):
  """
  Merges multiple measurements of symmetry-reduced hkl's.
  """

  def merging_reflection_table(self):
    '''Create a reflection table for storing merged reflections'''
    table = flex.reflection_table()
    table['miller_index'] = flex.miller_index()
    table['intensity']    = flex.double()
    table['esd']          = flex.double()
    table['rmsd']         = flex.double()
    table['multiplicity'] = flex.int()
    return table

  def distribute_reflection_table(self):
    '''Create a reflection table for storing reflections distributed over hkl chunks'''
    table = flex.reflection_table()
    table['miller_index_asymmetric']  = flex.miller_index()
    table['intensity.sum.value']      = flex.double()
    table['intensity.sum.variance']   = flex.double()
    return table

  def calc_reflection_intensity_stats(self, reflections):
    '''Do intensity statistics on reflection table'''

    multiplicity = len(reflections)
    assert multiplicity != 0

    stats = flex.mean_and_variance(reflections['intensity.sum.value'])
    propagated_esd = (flex.sum(reflections['intensity.sum.variance']) ** 0.5)/ multiplicity

    rmsd = 0.0
    if multiplicity > 1:
      rmsd = stats.unweighted_sample_standard_deviation()

    return {'intensity'     : stats.mean(),
            'esd'           : propagated_esd,
            'rmsd'          : rmsd,
            'multiplicity'  : multiplicity}

  def distribute_reflections_over_hkl_chunks(self, reflections):

    total_reflection_count = reflections.size()
    total_distributed_reflection_count = 0

    if total_reflection_count > 0:

      # set up two lists to be passed to the C++ extension: hkl's and chunk ids. It's basically a hash table to look up chunk ids by hkl's
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
    '''For alltoall, split each hkl chunk into N slices. This is needed to address the MPI alltoall memory problem.'''

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

  def generate_all_miller_indices(self):
    '''Given a unit cell, space group info, and resolution, generate all symmetry-reduced miller indices'''

    from cctbx import miller
    from cctbx.crystal import symmetry

    target_unit_cell = self.params.scaling.unit_cell
    target_space_group_info = self.params.scaling.space_group
    symm = symmetry(unit_cell = target_unit_cell, space_group_info = target_space_group_info)
    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=1000.0, d_min=self.params.merging.d_min)
    self.logger.log("Generated %d miller indices"%miller_set.indices().size())

    return miller_set

  def setup_hkl_chunks(self):
    '''Set up a list of reflection tables for distributing loaded reflections'''
    # generate all expected symmetry-reduced hkls
    full_miller_set = self.generate_all_miller_indices()

    # split the hkl set into chunks; the number of chunks is equal to the number of ranks
    import numpy as np
    self.hkl_split_set = np.array_split(full_miller_set.indices(), self.mpi_helper.size)

    # initialize a list of hkl chunks - reflection tables to store distributed reflections
    self.hkl_chunks = []
    for i in range(len(self.hkl_split_set)):
      self.hkl_chunks.append(self.distribute_reflection_table())

  def get_next_hkl_reflection_table(self, reflections):
    '''Slice reflection table into sub-tables by symmetry-reduced hkl's'''
    i_begin = 0
    hkl_ref = reflections[0].get('miller_index_asymmetric')

    for i in range(len(reflections)):
      hkl = reflections[i].get('miller_index_asymmetric')
      if hkl == hkl_ref:
        continue
      else:
        yield reflections[i_begin:i]
        i_begin = i
        hkl_ref = hkl

    yield reflections[i_begin:i+1]

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

  def output_merged_reflections(self, reflections):
    merged_reflections_file_path = self.params.output.output_dir + '/merge.out'
    merged_file = open(merged_reflections_file_path, 'w')

    for ref in reflections:
      merged_file.write("%s %f %f %f %d\n"%(ref.get('miller_index'), ref.get('intensity'), ref.get('esd'), ref.get('rmsd'), ref.get('multiplicity')))

    merged_file.close()

  def run(self, experiments, reflections):

    # set up hkl chunks to be used for all-to-all; every avialable rank participates in all-to-all, even if that rank doesn't load any data
    self.logger.log_step_time("SETUP_HKL_CHUNKS")
    self.setup_hkl_chunks()
    self.logger.log_step_time("SETUP_HKL_CHUNKS", True)

    # for the ranks, which have loaded the data, distribute the reflections over the hkl chunks
    self.logger.log_step_time("DISTRIBUTE")
    self.distribute_reflections_over_hkl_chunks(reflections=reflections)
    self.logger.log_step_time("DISTRIBUTE", True)

    # run all-to-all
    if self.params.parallel.a2a == 1:
      alltoall_reflections = self.get_reflections_from_alltoall()
    else: # if encountered alltoall memory problem, do alltoall on chunk slices
      alltoall_reflections = self.get_reflections_from_alltoall_sliced(number_of_slices=self.params.parallel.a2a)

    # sort reflections - necessary for the next step, merging of intensities
    self.logger.log_step_time("SORT")
    self.logger.log("Sorting consolidated reflection table...")
    if len(alltoall_reflections) > 0:
      alltoall_reflections.sort('miller_index_asymmetric')
    self.logger.log_step_time("SORT", True)

    # merge reflection intensities: calculate average, std, etc.
    self.logger.log_step_time("AVERAGE")
    self.logger.log("Averaging intensities...")
    all_rank_merged_reflections = self.merging_reflection_table()

    if len(alltoall_reflections) > 0:
      for hkl_reflection_table in self.get_next_hkl_reflection_table(reflections=alltoall_reflections):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=hkl_reflection_table)
        intensity_stats['miller_index'] = hkl_reflection_table[0].get('miller_index_asymmetric')
        all_rank_merged_reflections.append(intensity_stats)

    self.logger.log("Merged intensities for %d HKLs"%(all_rank_merged_reflections.size()))
    self.logger.log_step_time("AVERAGE", True)

    # gather all merged intensities at rank 0
    self.logger.log_step_time("GATHER")
    if self.mpi_helper.rank != 0:
      self.logger.log("Executing MPI gathering of all reflection tables at rank 0...")
    all_merged_reflection_tables = self.mpi_helper.comm.gather(all_rank_merged_reflections, root = 0)
    self.logger.log_step_time("GATHER", True)

    # rank 0: concatenate all merged intensities into the final table
    if self.mpi_helper.rank == 0:
      self.logger.log_step_time("MERGE")
      final_merged_reflection_table = self.merging_reflection_table()

      self.logger.log("Performing final merging of reflection tables received from all ranks...")
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      self.logger.log("Total merged HKLs: {}".format(final_merged_reflection_table.size()))
      self.logger.log_step_time("MERGE", True)

      # write the final merged reflection table out to an ASCII file
      self.logger.log_step_time("WRITE")
      self.output_merged_reflections(final_merged_reflection_table)
      self.logger.log_step_time("WRITE", True)

    return None, None

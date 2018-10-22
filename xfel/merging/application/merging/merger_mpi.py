from __future__ import division
from six.moves import range
import sys

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

"""
Functions for merging reflections
"""

from xfel.merging.application.input.file_loader_mpi import Script as Script_Base
from dials.array_family import flex

# create a reflection table for storing reflections distributed over hkl chunks
def distribute_reflection_table():
  table = flex.reflection_table()
  table['miller_index_asymmetric']  = flex.miller_index()
  table['intensity.sum.value']      = flex.double()
  table['intensity.sum.variance']   = flex.double()
  return table

# create a reflection table for storing merged reflectionsf
def merging_reflection_table():
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index()
  table['intensity']    = flex.double()
  table['esd']          = flex.double()
  table['rmsd']         = flex.double()
  table['multiplicity'] = flex.int()
  return table

class Script(Script_Base):
  '''A class for running the script.'''

  # add a "symmetry-reduced hkl" column to the reflection table
  def add_asu_miller_indices_column(self, experiments, reflections):

    from cctbx import miller
    from cctbx.crystal import symmetry
    from dials.array_family import flex

    new_reflections = flex.reflection_table()
    for experiment_id, experiment in enumerate(experiments):
      refls = reflections.select(reflections['id'] == experiment_id)
      SYM = symmetry(unit_cell = experiment.crystal.get_unit_cell(),
                                space_group = str(self.params.filter.unit_cell.value.target_space_group))
      mset = SYM.miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), indices=refls['miller_index'])

      refls['miller_index_original'] = refls['miller_index']
      del refls['miller_index']

      refls['miller_index_asymmetric'] = mset.map_to_asu().indices()
      new_reflections.extend(refls)

    return new_reflections

  # given a unit cell, space group info, and resolution, generate all symmetry-reduced miller indices
  def generate_all_miller_indices(self):

    from cctbx import miller
    from cctbx.crystal import symmetry

    self.log_step_time("GENERATE_MILLER_SET")
    unit_cell = self.params.filter.unit_cell.value.target_unit_cell
    space_group_info = self.params.filter.unit_cell.value.target_space_group
    symm = symmetry(unit_cell = unit_cell, space_group_info = space_group_info)
    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=1000.0, d_min=self.params.filter.resolution.d_min)
    self.mpi_log_write("\nRank 0 generated %d miller indices"%miller_set.indices().size())
    self.log_step_time("GENERATE_MILLER_SET", True)

    return miller_set

  # set up a list of reflection tables for distributing loaded reflections
  def setup_hkl_chunks(self):

    # generate all expected symmetry-reduced hkls
    full_miller_set = self.generate_all_miller_indices()

    # split the hkl set into chunks; the number of chunks is equal to the number of ranks
    self.log_step_time("SPLIT_MILLER_SET")
    import numpy as np
    self.hkl_split_set = np.array_split(full_miller_set.indices(), self.rank_count)
    self.log_step_time("SPLIT_MILLER_SET", True)

    # initialize a list of hkl chunks - reflection tables to store distributed reflections
    self.log_step_time("INIT_CHUNKS")
    self.hkl_chunks = []
    for i in range(len(self.hkl_split_set)):
      self.hkl_chunks.append(distribute_reflection_table())
    self.log_step_time("INIT_CHUNKS", True)

  # slice reflection table into sub-tables by symmetry-reduced hkl
  def get_next_hkl_reflection_table(self, reflections):

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

  # Generate an exact number of slices from a reflection table. Make slices as even as possible. If don't have enough reflections - generate empty tables.
  def get_next_reflection_table_slice(self, reflections, n_slices):
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
        yield distribute_reflection_table()

  # remove reflection table columns which are not relevant for merging
  def prune_reflection_columns(self, reflections):

    self.log_step_time("PRUNE")
    if reflections.size() > 0:
      all_keys = list()
      for key in reflections[0]:
        all_keys.append(key)

      for key in all_keys:
        if key != 'intensity.sum.value' and key != 'intensity.sum.variance' and key != 'miller_index_asymmetric':
          del reflections[key]
    self.log_step_time("PRUNE", True)

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

    self.mpi_log_write("\nRank %d managed to distribute %d out of %d reflections"%(self.rank, total_distributed_reflection_count, total_reflection_count))
    self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

    reflections.clear()

  def get_reflections_alltoall(self, comm):

    result_reflections = distribute_reflection_table()

    self.log_step_time("ALL-TO-ALL")
    self.mpi_log_write("\nRank %d executing MPI all-to-all..."%self.rank)

    received_hkl_chunks = comm.alltoall(self.hkl_chunks)

    self.mpi_log_write("\nRank %d after all-to-all received %d hkl chunks" % (self.rank, len(received_hkl_chunks)) )

    comm.barrier()
    self.log_step_time("ALL-TO-ALL", True)

    self.log_step_time("CONSOLIDATE")
    self.mpi_log_write("\nRank %d consolidating reflection tables..."%self.rank)

    for chunk in received_hkl_chunks:
      result_reflections.extend(chunk)

    comm.barrier()
    self.log_step_time("CONSOLIDATE", True)

    return result_reflections

  # For alltoall, split each hkl chunk into N slices. This is needed to address an MPI alltoall memory problem.
  def get_reflections_alltoall_sliced(self, comm, number_of_slices):

    result_reflections = distribute_reflection_table() # all reflections that the rank will receive from alltoall

    list_of_sliced_hkl_chunks = [] # if self.hkl_chunks is [A,B,C...], this list will be [[A1,A2..,An], [B1,B2,Bn], [C1,C2,Cn], ...], where n is the number of chunk slices
    for i in range(len(self.hkl_chunks)):
      hkl_chunk_slices = []
      for chunk_slice in self.get_next_reflection_table_slice(self.hkl_chunks[i], number_of_slices):
        hkl_chunk_slices.append(chunk_slice)
      list_of_sliced_hkl_chunks.append(hkl_chunk_slices)

    self.mpi_log_write("\nRANK %d is ready for all-to-all"%self.rank  )
    self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

    for j in range(number_of_slices):
      hkl_chunks_for_alltoall = list()
      for i in range(len(self.hkl_chunks)):
        hkl_chunks_for_alltoall.append(list_of_sliced_hkl_chunks[i][j]) # [Aj,Bj,Cj...]

      self.log_step_time("ALL-TO-ALL")
      self.mpi_log_write("\nRank %d executing MPI all-to-all..."%self.rank)
      self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

      received_hkl_chunks = comm.alltoall(hkl_chunks_for_alltoall)

      self.mpi_log_write("\nRank %d after all-to-all received %d hkl chunks" % (self.rank, len(received_hkl_chunks)) )

      comm.barrier()
      self.log_step_time("ALL-TO-ALL", True)

      self.log_step_time("CONSOLIDATE")
      self.mpi_log_write("\nRank %d consolidating reflection tables..."%self.rank)

      for chunk in received_hkl_chunks:
        result_reflections.extend(chunk)

      comm.barrier()
      self.log_step_time("CONSOLIDATE", True)

    return result_reflections

  # do intensity statistics on reflection table
  def calc_reflection_intensity_stats(self, reflections):

    multiplicity = len(reflections)
    assert multiplicity != 0

    stats = flex.mean_and_variance(reflections['intensity.sum.value'])
    propagated_esd = (flex.sum(reflections['intensity.sum.variance']) ** 0.5)/ multiplicity

    rmsd = 0.0
    if multiplicity > 1:
      rmsd = stats.unweighted_sample_standard_deviation()

    return {'average'       : stats.mean(),
            'esd'           : propagated_esd,
            'rmsd'          : rmsd,
            'multiplicity'  : multiplicity}

  # log elapsed time for a step in the merging process
  def log_step_time(self, step, finished=False):

    if not self.do_timing:
      return

    if not finished: # a step has started - just cache its start time and return
      if not step in self.timing_table:
        self.timing_table[step] = dict({'single_step':dict({'start':MPI.Wtime(), 'elapsed':0.0}),
                                        'cumulative': 0.0
                                       }
                                      )
      else: # if a step is executed repeatedly - re-use the existent step key
        self.timing_table[step]['single_step']['start'] = MPI.Wtime()
      return

    # a step has finished - calculate its elapsed and cumulative time (the latter is needed when the step is executed repeatedly)
    if not step in self.timing_table:
      assert False, "A step has finished, but doesn't have its starting time entry"
      return

    self.timing_table[step]['single_step']['elapsed'] = MPI.Wtime() - self.timing_table[step]['single_step']['start']
    self.timing_table[step]['cumulative'] += self.timing_table[step]['single_step']['elapsed']

    # output the elapsed and cumulative time
    if len(self.mpi_timing_file_path) == 0: # the input hasn't been received yet by this rank, so we don't know the file path yet
      return
    else:
        log_file = open(self.mpi_timing_file_path,'a')
        log_file.write("RANK %d %s: %f s %f s\n"%(self.rank, step, self.timing_table[step]['single_step']['elapsed'], self.timing_table[step]['cumulative']))
        log_file.close()

  def mpi_log_write(self, string):
    mpi_log_file_handle = open(self.mpi_log_file_path, 'a')
    mpi_log_file_handle.write(string)
    mpi_log_file_handle.close()

  def output_merged_reflections(self, reflections):
    merged_file = open(self.merged_reflections_file_path, 'w')

    for ref in reflections:
      merged_file.write("%s %f %f %f %d\n"%(ref.get('miller_index'), ref.get('intensity'), ref.get('esd'), ref.get('rmsd'), ref.get('multiplicity')))

    merged_file.close()

  ###################
  # SCRIPT RUN METHOD
  def run(self, comm, timing=True):

    # create timing table
    self.do_timing = timing
    if self.do_timing:
      self.timing_table = dict()

    # get rank and rank count
    self.rank = comm.Get_rank()
    self.rank_count = comm.Get_size()

    self.log_step_time("TOTAL")

    transmitted = None
    if self.rank == 0:
      # read and parse phil
      self.initialize()
      self.validate()

      # initialize output paths
      if self.do_timing:
        self.mpi_timing_file_path = self.params.output.output_dir + '/timing_%06d_%06d.out'%(self.rank_count, self.rank)
      self.mpi_log_file_path = self.params.output.output_dir + '/rank_%06d_%06d.out'%(self.rank_count, self.rank)
      self.merged_reflections_file_path = self.params.output.output_dir + '/merge.out'

      ###########################
      # GET JSON/PICKLE FILE LIST
      self.log_step_time("LIST_FILES")
      file_list = self.get_list()
      self.mpi_log_write("\nRank 0 generated a list of %d file items"%len(file_list))
      self.log_step_time("LIST_FILES", True)

      # prepare for transmitting the job to all ranks
      self.mpi_log_write("\nRank 0 transmitting file list of length %d"%(len(file_list)))
      transmitted = dict(params = self.params, options = self.options, file_list = file_list)

    #########################################
    # BROADCAST WORK FROM RANK 0 TO ALL RANKS
    comm.barrier()
    self.log_step_time("BROADCAST")
    transmitted = comm.bcast(transmitted, root = 0)

    if self.rank == 0:
      del file_list[:]

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['file_list'][self.rank::self.rank_count]

    if self.do_timing:
      self.mpi_timing_file_path = self.params.output.output_dir + '/timing_%06d_%06d.out'%(self.rank_count, self.rank)
    self.mpi_log_file_path = self.params.output.output_dir + '/rank_%06d_%06d.out'%(self.rank_count, self.rank)

    self.mpi_log_write ("\nRank %d received a file list of %d json-pickle file pairs" % (self.rank, len(new_file_list)))
    comm.barrier()
    self.log_step_time("BROADCAST", True)

    ################################################################################################
    # EACH RANK: SET UP HKL CHUNKS TO BE USED FOR DISTRIBUTING LOADED REFLECTIONS AND FOR ALL-TO-ALL
    self.setup_hkl_chunks() # even if a rank got no files to load, it still has to participate in all-to-all

    # initialize data loading statistics
    rank_experiment_count = 0
    rank_image_count = 0
    rank_reflection_count = 0

    if len(new_file_list) > 0:
      self.mpi_log_write("\nRank %d: first file to load is: %s"%(self.rank, str(new_file_list[0])))

      ######################
      # EACH RANK: LOAD DATA
      self.log_step_time("LOAD")
      experiments, reflections = self.load_data(new_file_list)

      self.mpi_log_write ('\nRank %d has read %d experiments consisting of %d reflections'%(self.rank, len(experiments), len(reflections)))
      self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

      comm.barrier()
      self.log_step_time("LOAD", True)

      ###############################################################################
      # EACH RANK: GET THE TOTAL NUMBER OF LOADED EXPERIMENTS, IMAGES AND REFLECTIONS
      rank_experiment_count = len(experiments)

      # count number of images
      all_imgs = []
      for iset in experiments.imagesets():
        all_imgs.extend(iset.paths())
      rank_image_count = len(set(all_imgs))

      rank_reflection_count = len(reflections)

      #####################################################
      # EACH RANK: ADD A COLUMN WITH SYMMETRY-REDUCED HKL's
      reflections = self.add_asu_miller_indices_column(experiments, reflections)

      #########################################################################
      # EACH RANK: PRUNE REFLECTION COLUMNS, WHICH ARE NOT RELEVANT FOR MERGING
      self.prune_reflection_columns(reflections=reflections)

      #######################################################
      # EACH RANK: DISTRIBUTE REFLECTIONS OVER ALL HKL CHUNKS
      self.log_step_time("DISTRIBUTE")
      self.distribute_reflections_over_hkl_chunks(reflections=reflections)
      comm.barrier()
      self.log_step_time("DISTRIBUTE", True)
    else:
      self.mpi_log_write ("\nRank %d received no data" % self.rank)
      comm.barrier()
      comm.barrier()

    ########################################
    # MPI-REDUCE ALL DATA LOADING STATISTICS
    total_experiment_count  = comm.reduce(rank_experiment_count, MPI.SUM, 0)
    total_image_count       = comm.reduce(rank_image_count, MPI.SUM, 0)
    total_reflection_count  = comm.reduce(rank_reflection_count, MPI.SUM, 0)
    max_reflection_count    = comm.reduce(rank_reflection_count, MPI.MAX, 0)
    min_reflection_count    = comm.reduce(rank_reflection_count, MPI.MIN, 0)

    ########################################################################################
    # RANK 0: LOG ALL DATA LOADING STATISTICS, INITIALIZE THE TOTAL MERGING REFLECTION TABLE
    if self.rank == 0:
      self.mpi_log_write('\nAll ranks have read %d experiments'%total_experiment_count)
      self.mpi_log_write('\nAll ranks have read %d images'%total_image_count)
      self.mpi_log_write('\nAll ranks have read %d reflections'%total_reflection_count)
      self.mpi_log_write('\nThe maximum number of reflections loaded per rank is: %d reflections'%max_reflection_count)
      self.mpi_log_write('\nThe minimum number of reflections loaded per rank is: %d reflections'%min_reflection_count)

      final_merged_reflection_table = merging_reflection_table()

    #############################################################
    # EACH RANK: GET A REFLECTION TABLE FOR MERGING BY ALL-TO-ALL
    if self.params.parallel.a2a == 1:
      alltoall_reflections = self.get_reflections_alltoall(comm=comm)
    else: # if encountered alltoall memory problem, do alltoall on chunk slices
      alltoall_reflections = self.get_reflections_alltoall_sliced(comm=comm, number_of_slices=self.params.parallel.a2a)

    #####################################################
    # EACH RANK: SORT THE REFLECTION TABLE BEFORE MERGING
    self.log_step_time("SORT")
    self.mpi_log_write("\nRank %d sorting consolidated reflection table..."%self.rank)

    if len(alltoall_reflections) > 0:
      alltoall_reflections.sort('miller_index_asymmetric')

    comm.barrier()
    self.log_step_time("SORT", True)

    ####################################################
    # EACH RANK: DO STATISTICS ON REFLECTION INTENSITIES
    self.log_step_time("AVERAGE")
    self.mpi_log_write("\nRank %d doing intensity statistics..."%self.rank)
    all_rank_merged_reflections = merging_reflection_table()

    if len(alltoall_reflections) > 0:
      for hkl_reflection_table in self.get_next_hkl_reflection_table(reflections=alltoall_reflections):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=hkl_reflection_table)
        all_rank_merged_reflections.append({'miller_index': hkl_reflection_table[0].get('miller_index_asymmetric'),
                                            'intensity': intensity_stats['average'],
                                            'esd' : intensity_stats['esd'],
                                            'rmsd' : intensity_stats['rmsd'],
                                            'multiplicity': intensity_stats['multiplicity']})

    self.mpi_log_write ("\nRank %d merged intensities for %d HKLs"%(self.rank, all_rank_merged_reflections.size()))
    comm.barrier()
    self.log_step_time("AVERAGE", True)

    ####################################################
    # EACH RANK: SEND MERGED REFLECTION TABLES TO RANK 0
    self.log_step_time("GATHER")
    if self.rank != 0:
      self.mpi_log_write("\nRank %d executing MPI gathering of all reflection tables at rank 0..."%self.rank)
    all_merged_reflection_tables = comm.gather(all_rank_merged_reflections, root = 0)
    comm.barrier()
    self.log_step_time("GATHER", True)

    ####################################
    # RANK 0: DO FINAL MERGING OF TABLES
    if self.rank == 0:
      self.log_step_time("MERGE")
      self.mpi_log_write ("\nRank 0 doing final merging of reflection tables received from all ranks...")
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      self.mpi_log_write("\nRank 0 total merged HKLs: {}".format(final_merged_reflection_table.size()))
      self.log_step_time("MERGE", True)

      # write the final merged reflection table out to an ASCII file
      self.log_step_time("WRITE")
      self.output_merged_reflections(final_merged_reflection_table)
      self.log_step_time("WRITE", True)

    comm.barrier()
    self.log_step_time("TOTAL", True)

    MPI.Finalize()

    return

if __name__ == '__main__':

  from mpi4py import MPI

  comm = MPI.COMM_WORLD

  script = Script()

  result = script.run(comm=comm)
  if result is None:
    sys.exit(1)

  print ("OK")

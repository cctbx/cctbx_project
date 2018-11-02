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

from xfel.merging.application.phil.phil import Script as Script_Base
from xfel.merging.application.input.file_loader import file_loader
from xfel.merging.application.input.calculate_file_load import file_load_calculator

from dials.array_family import flex

def distribute_reflection_table():
  '''Create a reflection table for storing reflections distributed over hkl chunks'''
  table = flex.reflection_table()
  table['miller_index_asymmetric']  = flex.miller_index()
  table['intensity.sum.value']      = flex.double()
  table['intensity.sum.variance']   = flex.double()
  return table

def merging_reflection_table():
  '''Create a reflection table for storing merged reflections'''
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index()
  table['intensity']    = flex.double()
  table['esd']          = flex.double()
  table['rmsd']         = flex.double()
  table['multiplicity'] = flex.int()
  return table

class Script(Script_Base):
  '''A class for running the script.'''

  def load_data(self, file_list):
    from dxtbx.model.experiment_list import ExperimentList
    from dials.array_family import flex
    all_experiments = ExperimentList()
    all_reflections = flex.reflection_table()

    loader = file_loader(self.params)

    for experiments_filename, reflections_filename in file_list:

      experiments, reflections = loader.load_data(experiments_filename, reflections_filename)

      for experiment_id, experiment in enumerate(experiments):
        all_experiments.append(experiment)
        refls = reflections.select(reflections['id'] == experiment_id)

        if len(refls) > 0:
          refls['id'] = flex.int(len(refls), len(all_experiments)-1)
          all_reflections.extend(refls)

    return all_experiments, all_reflections

  def add_asu_miller_indices_column(self, experiments, reflections):
    '''Add a "symmetry-reduced hkl" column to the reflection table'''
    from cctbx import miller
    from cctbx.crystal import symmetry

    # assert that the space group used to integrate the data is consistent with the input space group for merging
    target_symmetry = symmetry(unit_cell = self.params.filter.unit_cell.value.target_unit_cell, space_group = str(self.params.filter.unit_cell.value.target_space_group))
    target_space_group = target_symmetry.space_group()
    target_patterson_group_info = target_space_group.build_derived_patterson_group().info().symbol_and_number()

    for experiment in experiments:
      experiment_patterson_group_info = experiment.crystal.get_space_group().build_derived_patterson_group().info().symbol_and_number()
      if target_patterson_group_info != experiment_patterson_group_info:
        assert False, "Target patterson group %s is different from an experiment patterson group %s"%(target_patterson_group_info, experiment_patterson_group_info)

    target_miller = target_symmetry.miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), indices=reflections['miller_index'])

    # change the name of the original hkl column
    reflections['miller_index_original'] = reflections['miller_index']
    del reflections['miller_index']

    # add a new hkl column
    reflections['miller_index_asymmetric'] = target_miller.map_to_asu().indices()

  def generate_all_miller_indices(self):
    '''Given a unit cell, space group info, and resolution, generate all symmetry-reduced miller indices'''

    from cctbx import miller
    from cctbx.crystal import symmetry

    unit_cell = self.params.filter.unit_cell.value.target_unit_cell
    space_group_info = self.params.filter.unit_cell.value.target_space_group
    symm = symmetry(unit_cell = unit_cell, space_group_info = space_group_info)
    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=1000.0, d_min=self.params.filter.resolution.d_min)
    self.mpi_log_write("\nRank 0 generated %d miller indices"%miller_set.indices().size())

    return miller_set

  def setup_hkl_chunks(self):
    '''Set up a list of reflection tables for distributing loaded reflections'''
    # generate all expected symmetry-reduced hkls
    full_miller_set = self.generate_all_miller_indices()

    # split the hkl set into chunks; the number of chunks is equal to the number of ranks
    import numpy as np
    self.hkl_split_set = np.array_split(full_miller_set.indices(), self.rank_count)

    # initialize a list of hkl chunks - reflection tables to store distributed reflections
    self.hkl_chunks = []
    for i in range(len(self.hkl_split_set)):
      self.hkl_chunks.append(distribute_reflection_table())

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
        yield distribute_reflection_table()

  def prune_reflection_columns(self, reflections):
    '''Remove reflection table columns which are not relevant for merging'''
    if reflections.size() > 0:
      all_keys = list()
      for key in reflections[0]:
        all_keys.append(key)

      for key in all_keys:
        if key != 'intensity.sum.value' and key != 'intensity.sum.variance' and key != 'miller_index_asymmetric' and key != 'id':
          del reflections[key]

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

  def get_reflections_from_alltoall(self, comm):

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

  def get_reflections_from_alltoall_sliced(self, comm, number_of_slices):
    '''For alltoall, split each hkl chunk into N slices. This is needed to address the MPI alltoall memory problem.'''

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

  def calc_reflection_intensity_stats(self, reflections):
    '''Do intensity statistics on reflection table'''

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

  def log_step_time(self, step, finished=False):
    '''Log elapsed time for a step in the merging process'''

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
    self.rank       = comm.Get_rank()
    self.rank_count = comm.Get_size()

    self.log_step_time("TOTAL")

    # prepare a list of input params and data files to be transmitted to all worker ranks
    transmitted = {}
    if self.rank == 0:
      self.log_step_time("PARSE_INPUT_PARAMS")

      # read and parse phil
      self.initialize()
      self.validate()

      # initialize output paths; doing it here for rank 0 only - to be able to log timing, etc., before every rank initializes the paths
      if self.do_timing:
        self.mpi_timing_file_path = self.params.output.output_dir + '/timing_%06d_%06d.out'%(self.rank_count, self.rank)
      self.mpi_log_file_path = self.params.output.output_dir + '/rank_%06d_%06d.out'%(self.rank_count, self.rank)
      self.merged_reflections_file_path = self.params.output.output_dir + '/merge.out'

      self.log_step_time("PARSE_INPUT_PARAMS", True)

      ###########################
      # GET JSON/PICKLE FILE LIST
      self.log_step_time("LIST_FILES")
      loader = file_loader(self.params)
      file_list = loader.filename_lister()

      self.log_step_time("CALCULATE_FILE_LOAD")
      load_calculator = file_load_calculator(self.params, file_list)
      rank_files = load_calculator.calculate_file_load(self.rank_count)

      total_file_pairs = 0
      for key, value in rank_files.items():
        total_file_pairs += len(value)
      self.mpi_log_write("\nRank 0 generated a list of %d file items for %d ranks"%(total_file_pairs,len(rank_files)))
      self.log_step_time("CALCULATE_FILE_LOAD", True)

      self.log_step_time("LIST_FILES", True)

      # prepare for transmitting the job to all ranks
      self.mpi_log_write("\nRank 0 transmitting file list of length %d"%(len(file_list)))
      transmitted = dict(params = self.params, options = self.options, rank_files = rank_files)

    #########################################
    # BROADCAST WORK FROM RANK 0 TO ALL RANKS
    comm.barrier()
    self.log_step_time("BROADCAST")

    transmitted = comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['rank_files'][self.rank]

    if self.do_timing:
      self.mpi_timing_file_path = self.params.output.output_dir + '/timing_%06d_%06d.out'%(self.rank_count, self.rank)
    self.mpi_log_file_path = self.params.output.output_dir + '/rank_%06d_%06d.out'%(self.rank_count, self.rank)

    self.mpi_log_write ("\nRank %d received a file list of %d json-pickle file pairs" % (self.rank, len(new_file_list)))
    comm.barrier()
    self.log_step_time("BROADCAST", True)

    ################################################################################################
    # EACH RANK: SET UP HKL CHUNKS TO BE USED FOR DISTRIBUTING LOADED REFLECTIONS AND FOR ALL-TO-ALL
    self.log_step_time("SETUP_HKL_CHUNKS")
    self.setup_hkl_chunks() # even if a rank got no files to load, it still has to participate in all-to-all
    comm.barrier()
    self.log_step_time("SETUP_HKL_CHUNKS", True)

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
      self.log_step_time("LOAD_STATISTICS")
      rank_experiment_count = len(experiments)

      # count number of images
      all_imgs = []
      for iset in experiments.imagesets():
        all_imgs.extend(iset.paths())
      rank_image_count = len(set(all_imgs))

      rank_reflection_count = len(reflections)
      comm.barrier()
      self.log_step_time("LOAD_STATISTICS", True)

      #####################################################
      # EACH RANK: ADD A COLUMN WITH SYMMETRY-REDUCED HKL's
      self.log_step_time("ADD_ASU_HKL_COLUMN")
      #reflections = self.add_asu_miller_indices_column(experiments, reflections)
      self.add_asu_miller_indices_column(experiments, reflections)
      comm.barrier()
      self.log_step_time("ADD_ASU_HKL_COLUMN", True)

      #########################################################################
      # EACH RANK: PRUNE REFLECTION COLUMNS, WHICH ARE NOT RELEVANT FOR MERGING
      self.log_step_time("PRUNE")
      self.prune_reflection_columns(reflections=reflections)
      comm.barrier()
      self.log_step_time("PRUNE", True)
    else:
      reflections = distribute_reflection_table()
      self.mpi_log_write ("\nRank %d received no data" % self.rank)
      comm.barrier()
      comm.barrier()
      comm.barrier()
      comm.barrier()

    ########################################
    # MPI-REDUCE ALL DATA LOADING STATISTICS
    self.log_step_time("REDUCE_LOAD_STATS")
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

    comm.barrier()
    self.log_step_time("REDUCE_LOAD_STATS", True)

    #######################################################
    # EACH RANK: DISTRIBUTE REFLECTIONS OVER ALL HKL CHUNKS
    self.log_step_time("DISTRIBUTE")
    self.distribute_reflections_over_hkl_chunks(reflections=reflections)
    comm.barrier()
    self.log_step_time("DISTRIBUTE", True)

    #############################################################
    # EACH RANK: GET A REFLECTION TABLE FOR MERGING BY ALL-TO-ALL
    if self.params.parallel.a2a == 1:
      alltoall_reflections = self.get_reflections_from_alltoall(comm=comm)
    else: # if encountered alltoall memory problem, do alltoall on chunk slices
      alltoall_reflections = self.get_reflections_from_alltoall_sliced(comm=comm, number_of_slices=self.params.parallel.a2a)

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
      final_merged_reflection_table = merging_reflection_table()

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

  from libtbx.mpi4py import MPI

  comm = MPI.COMM_WORLD

  script = Script()

  result = script.run(comm=comm)
  if result is None:
    sys.exit(1)

  print ("OK")

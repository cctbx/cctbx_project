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

def merging_reflection_table():
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index()
  table['intensity']    = flex.double()
  table['esd']          = flex.double()
  table['rmsd']         = flex.double()
  table['multiplicity'] = flex.int()
  return table

from itertools import islice

def split_dictionary(data, stride):
    it = iter(data)
    for i in range(0, len(data), stride):
        yield {k:data[k] for k in islice(it, stride)}

class Script(Script_Base):
  '''A class for running the script.'''

  def sort_reflections_by_asu_miller_index(self, experiments, reflections):
    from cctbx import miller
    from cctbx.crystal import symmetry
    from dials.array_family import flex

    sorted_reflections = flex.reflection_table()

    for experiment_id, experiment in enumerate(experiments):

      refls = reflections.select(reflections['id'] == experiment_id)

      SYM = symmetry(unit_cell = experiment.crystal.get_unit_cell(), space_group = str(self.params.filter.unit_cell.value.target_space_group))
      mset = SYM.miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), indices=refls['miller_index'])

      refls['miller_index_original'] = refls['miller_index']
      del refls['miller_index']

      refls['miller_index_asymmetric'] = mset.map_to_asu().indices()
      sorted_reflections.extend(refls)

    sorted_reflections.sort('miller_index_asymmetric')

    return sorted_reflections

  # given a unit cell, space group info, and resolution, generate all symmetry-reduced miller indices
  def generate_all_miller_indices(self):

    from cctbx import miller
    from cctbx.crystal import symmetry

    unit_cell = self.params.filter.unit_cell.value.target_unit_cell
    space_group_info = self.params.filter.unit_cell.value.target_space_group
    symm = symmetry(unit_cell = unit_cell, space_group_info = space_group_info)
    miller_set = symm.build_miller_set(anomalous_flag=(not self.params.merging.merge_anomalous), d_max=1000.0, d_min=self.params.filter.resolution.d_min)

    return miller_set

  # return a reflection table which stores measurements of a given symmetry-reduced hkl
  def get_hkl_chunk(self, hkl, chunks):

    for chunk in chunks:
      if hkl in chunk:
        return chunk[hkl]

    return None

  # Knowing that the input reflection table is sorted on hkl's, identify blocks of symmetry-reduced hkl's. Append each block to the corresponding reflection table in the input list of hkl chunks.
  def distribute_reflections_over_hkl_chunks(self, reflections, chunks):
    total_reflection_count = reflections.size()
    total_distributed_reflection_count = 0

    if total_reflection_count > 0:

      # remove columns which are not relevant for merging
      all_keys = list()
      for key in reflections[0]:
        all_keys += key
      for key in all_keys:
        if( key != 'intensity.sum.value' and key != 'intensity.sum.variance' and key != 'intensity.sum.value' ):
          del reflections[key]

      # loop through the table to identify blocks of symmetry-reduced hkls
      index_min = 0 # the first index in the reflection list, where a new "current" HKL starts
      hkl_cur = reflections[0].get('miller_index_asymmetric') # current HKL
      hkl = None # HKL corresponding to the looping index

      for i in range(0, total_reflection_count+1): # +1 is needed to handle the last HKL block in the list
        if i < total_reflection_count:
          hkl = reflections[i].get('miller_index_asymmetric')
        else:
          hkl = None

        if hkl != hkl_cur: # a new HKL block has started

          # get a chunk reflection table to append the previous HKL block
          chunk_hkl_reflection_table = self.get_hkl_chunk(hkl=hkl_cur, chunks=chunks)

          # append the previous HKL block to its corresponding chunk
          if chunk_hkl_reflection_table != None:
            chunk_hkl_reflection_table.extend(reflections[index_min:i])

            total_distributed_reflection_count += (i - index_min)

          # begin tracking new HKL block
          hkl_cur = hkl
          index_min = i

    self.mpi_log_write("\nRank %d managed to distribute %d out of %d reflections"%(self.rank, total_distributed_reflection_count, len(reflections)))

    self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

    reflections.clear()

  # Given a reflection table, do the intensity statistics. Use case: multiple measurements of the same symmetry-reduced hkl.
  def calc_reflection_intensity_stats(self, reflections):

    multiplicity = len(reflections)
    if( multiplicity == 0 ):
      return dict()

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

    if not finished: # a step has started - just cache its start time and leave
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
    if( len(self.mpi_timing_file_path) == 0 ): # the input hasn't been received yet by this rank, so we don't know the file path yet
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

    # ger rank and rank count
    self.rank = comm.Get_rank()
    self.rank_count = comm.Get_size()

    self.log_step_time("TOTAL")

    if( self.rank == 0 ):
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
      self.log_step_time("SETUP")
      file_list = self.get_list()
      self.mpi_log_write("\nRank 0 generated a list of %d file items"%len(file_list))
      #print("\nRank 0 generated a list of %d file items"%len(file_list))

      #############################
      # GENERATE ALL MILLER INDICES
      full_miller_set = self.generate_all_miller_indices()
      self.mpi_log_write("\nRank 0 generated %d miller indices"%full_miller_set.indices().size())
      #print("\nRank 0 generated %d miller indices"%full_miller_set.indices().size())

      ###################################################################
      # SPLIT MILLER INDICES INTO CHUNKS TO BE DISTRIBUTED OVER THE RANKS
      import numpy as np
      split_set = np.array_split(full_miller_set.indices(), self.rank_count)

      # within chunks, create an empty reflection table for every hkl
      list_of_chunks_of_hkl_dictionaries = []
      for chunk in split_set:
        chunk_of_hkl_dictionaries = dict()
        for hkl_array in chunk:
          hkl = (hkl_array[0], hkl_array[1], hkl_array[2])
          chunk_of_hkl_dictionaries[hkl] = flex.reflection_table()
        list_of_chunks_of_hkl_dictionaries.append(chunk_of_hkl_dictionaries)

      self.mpi_log_write("\nRank 0 generated %d chunks of hkl containers"%len(list_of_chunks_of_hkl_dictionaries))
      self.log_step_time("SETUP", True)

      # prepare for transmitting the job to all ranks
      self.mpi_log_write("\nRank 0 transmitting file list of length %d"%(len(file_list)))
      transmitted = dict(params = self.params, options = self.options, file_list = file_list, list_of_chunks_of_hkl_dictionaries = list_of_chunks_of_hkl_dictionaries)
    else:
      transmitted = None

    #########################################
    # BROADCAST WORK FROM RANK 0 TO ALL RANKS
    comm.barrier()
    self.log_step_time("BROADCAST")
    transmitted = comm.bcast(transmitted, root = 0)

    if self.rank == 0:
      del file_list[:]
      del list_of_chunks_of_hkl_dictionaries[:]
      #full_miller_set.clear()
      del split_set[:]

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['file_list'][self.rank::self.rank_count]
    chunks = transmitted['list_of_chunks_of_hkl_dictionaries']

    if self.do_timing:
      self.mpi_timing_file_path = self.params.output.output_dir + '/timing_%06d_%06d.out'%(self.rank_count, self.rank)
    self.mpi_log_file_path = self.params.output.output_dir + '/rank_%06d_%06d.out'%(self.rank_count, self.rank)

    self.mpi_log_write ("\nRank %d received a file list of %d json-pickle file pairs and %d chunks of hkl containers" % (self.rank, len(new_file_list), len(chunks)))
    comm.barrier()
    self.log_step_time("BROADCAST", True)

    rank_experiment_count = 0
    rank_image_count = 0
    rank_reflection_count = 0

    if( len(new_file_list) > 0 ):
      self.mpi_log_write("\nRank %d: first file to load is: %s" % (self.rank, str(new_file_list[0])))

      ######################
      # EACH RANK: LOAD DATA
      comm.barrier()
      self.log_step_time("LOAD")
      experiments, reflections = self.load_data(new_file_list)
      self.mpi_log_write ('\nRank %d has read %d experiments consisting of %d reflections'%(self.rank, len(experiments), len(reflections)))

      self.mpi_log_write("\nRANK %d: %s"%(self.rank, get_memory_usage()))

      # count number of images
      all_imgs = []
      for iset in experiments.imagesets():
        all_imgs.extend(iset.paths())

      self.log_step_time("LOAD", True)

      ######################
      # EACH RANK: SORT DATA
      comm.barrier()
      self.log_step_time("SORT")
      if( len(reflections) > 0 ):
        reflections = self.sort_reflections_by_asu_miller_index(experiments, reflections)
      self.log_step_time("SORT", True)

      comm.barrier()

      ####################################################################
      # GET THE TOTAL NUMBER OF LOADED EXPERIMENTS, IMAGES AND REFLECTIONS
      rank_experiment_count = len(experiments)
      rank_image_count = len(set(all_imgs))
      rank_reflection_count = len(reflections)

      #######################################################
      # EACH RANK: DISTRIBUTE REFLECTIONS OVER ALL HKL CHUNKS
      self.log_step_time("DISTRIBUTE")

      self.distribute_reflections_over_hkl_chunks(reflections=reflections, chunks=chunks)

      comm.barrier()
      self.log_step_time("DISTRIBUTE", True)
    else:
      self.mpi_log_write ("\nRank %d received no data" % self.rank)
      comm.barrier()
      comm.barrier()
      comm.barrier()
      comm.barrier()

    ###################################################
    # MPI-REDUCE THE TOTAL NUMBER OF LOADED EXPERIMENTS
    total_experiment_count = comm.reduce(rank_experiment_count, MPI.SUM, 0)
    if( self.rank == 0 ):
      self.mpi_log_write('\nAll ranks have read %d experiments'%total_experiment_count)

    ###################################################
    # MPI-REDUCE THE TOTAL NUMBER OF LOADED IMAGES
    total_image_count = comm.reduce(rank_image_count, MPI.SUM, 0)
    if( self.rank == 0 ):
      self.mpi_log_write('\nAll ranks have read %d images'%total_image_count)

    ###################################################
    # MPI-REDUCE THE TOTAL NUMBER OF LOADED REFLECTIONS
    total_reflection_count = comm.reduce(rank_reflection_count, MPI.SUM, 0)
    if( self.rank == 0 ):
      self.mpi_log_write('\nAll ranks have read %d reflections'%total_reflection_count)

    ###################################################
    # MPI-REDUCE THE MAX NUMBER OF LOADED REFLECTIONS
    max_reflection_count = comm.reduce(rank_reflection_count, MPI.MAX, 0)
    if( self.rank == 0 ):
      self.mpi_log_write('\nThe maximum number of reflections loaded per rank is: %d reflections'%max_reflection_count)

    ###################################################
    # MPI-REDUCE THE MIN NUMBER OF LOADED REFLECTIONS
    min_reflection_count = comm.reduce(rank_reflection_count, MPI.MIN, 0)
    if( self.rank == 0 ):
      self.mpi_log_write('\nThe minimum number of reflections loaded per rank is: %d reflections'%min_reflection_count)

    #################################################################################################################################
    # BREAK THE HKL CHUNKS INTO N PORTIONS AND EXECUTE THE (ALLTOALL-CONSOLIDATE-AVERAGE-MERGE) SEQUENCE N times - 1 FOR EACH PORTION

    # prepare the final reflection table for rank 0 to collect all merged HKLs
    if self.rank == 0:
      final_merged_reflection_table = merging_reflection_table()

    # To address memory problem during MPI alltoall execution, split each hkl chunk into N portions
    number_of_chunk_portions = self.params.parallel.a2a
    list_of_chunks_portions = list()
    for j in range(0, number_of_chunk_portions):
      list_of_chunks_portions.append(list())

    for chunk in chunks:
      i = 0
      for chunk_piece in split_dictionary(chunk, int(len(chunk)/number_of_chunk_portions)+1):
        list_of_chunks_portions[i].append(chunk_piece)
        i += 1

    for portion in list_of_chunks_portions:
      if len(portion) == 0:
        continue

      #######################################################################
      # EACH RANK: RECEIVE ALL HKL CHUNKS (WITH THE SAME HKLS) FROM ALL RANKS
      self.log_step_time("ALL-TO-ALL")
      self.mpi_log_write("\nRank %d executing MPI all-to-all..."%self.rank)

      received_chunks = comm.alltoall(portion)

      received_chunks_count = len(received_chunks)
      self.mpi_log_write ("\nRank %d after all-to-all received %d chunks of hkl containers" % (self.rank, received_chunks_count) )

      comm.barrier()
      self.log_step_time("ALL-TO-ALL", True)

      ######################################################################
      # EACH RANK: CONSOLIDATE ALL REFLECTION LISTS FROM ALL RECEIVED CHUNKS
      self.log_step_time("CONSOLIDATE")
      self.mpi_log_write("\nRank %d consolidating reflection lists..."%self.rank)
      for hkl in received_chunks[0]:
        for chunk_index in range(1,received_chunks_count):
          received_chunks[0][hkl].extend(received_chunks[chunk_index][hkl])
      comm.barrier()
      self.log_step_time("CONSOLIDATE", True)

      ####################################################
      # EACH RANK: DO STATISTICS ON REFLECTION INTENSITIES
      self.log_step_time("AVERAGE")
      self.mpi_log_write("\nRank %d doing intensity statistics..."%self.rank)
      all_rank_merged_reflections = merging_reflection_table()

      for hkl in received_chunks[0]:
        if( len(received_chunks[0][hkl]) > 0 ):
          intensity_stats = self.calc_reflection_intensity_stats(reflections=received_chunks[0][hkl])
          all_rank_merged_reflections.append({'miller_index': hkl,
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
        #final_merged_reflection_table = merging_reflection_table()
        for table in all_merged_reflection_tables:
          final_merged_reflection_table.extend(table)
        self.mpi_log_write("\nRank 0 total merged HKLs: {}".format(final_merged_reflection_table.size()))
        self.log_step_time("MERGE", True)

    # write the final merged reflection table out to an ASCII file
    if self.rank == 0:
      self.mpi_log_write("\nRank 0 final total merged HKLs: {}".format(final_merged_reflection_table.size()))
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

from __future__ import division
import sys

"""
Functions for merging reflections
"""

from xfel.merging.application.input.file_loader_mpi import Script as Script_Base
from dials.array_family import flex

def merging_reflection_table():
  table = flex.reflection_table()
  table['miller_index'] = flex.miller_index()
  table['intensity']    = flex.double()
  table['multiplicity'] = flex.int()
  table['esd']          = flex.double()
  table['rmsd']         = flex.double()
  return table

class Script(Script_Base):
  '''A class for running the script.'''

  # given unit cell, space group info, and resolution, generate all miller indices
  def generate_all_miller_indices(self):

    from cctbx import miller
    from cctbx.crystal import symmetry

    unit_cell = self.params.filter.unit_cell.value.target_unit_cell
    space_group_info = self.params.filter.unit_cell.value.target_space_group
    symm = symmetry(unit_cell = unit_cell, space_group_info = space_group_info)
    miller_set = symm.build_miller_set(anomalous_flag=False, d_max=1000.0, d_min=self.params.filter.resolution.d_min)

    return miller_set

  # given an hkl and all chunks of hkl's, return the chunk containing the given hkl
  def get_hkl_chunk(self, hkl, chunks):

    for chunk in chunks:
      if( hkl in chunk ):
        return chunk

    #print ("\nHKL: %s is not present in any of the chunks"%str(hkl))

    return None

  # For each reflection find an hkl chunk, containing the reflection's hkl; append that reflection to the list of reflections for that hkl
  # Here we are assuming the input reflection list is sorted on miller indices
  def distribute_reflections_over_hkl_chunks(self, reflections, chunks):

    total_reflection_count = reflections.size()
    total_distributed_reflection_count = 0

    if total_reflection_count > 0:

      index_min = 0 # the first index in the reflection list, where a new "current" HKL starts, assuming the list is HKL-sorted
      hkl_cur = reflections[0].get('miller_index_asymmetric') # current HKL
      hkl = None # HKL corresponding to the looping index

      for i in range (0, total_reflection_count+1): # +1 is needed to handle the last "current" HKL in the list
        if i < total_reflection_count:
          hkl = reflections[i].get('miller_index_asymmetric')
        else:
          hkl = None

        if hkl != hkl_cur: # a new current HKL has started

          # get a chunk for the new current HKL
          chunk = self.get_hkl_chunk(hkl=hkl_cur, chunks=chunks)

          # append current HKL to its corresponding chunk
          if chunk != None:
            for j in range(index_min, i):
              chunk[hkl_cur].append(reflections[j])
            total_distributed_reflection_count += (i - index_min)

          # begin tracking new current HKL
          hkl_cur = hkl
          index_min = i

    print("Rank %d managed to distribute %d of %d reflections"%(self.rank, total_distributed_reflection_count, len(reflections)))

    '''
    number_of_allocated_reflections = dict()
    for chunk_index in range(0, len(chunks)):
      number_of_allocated_reflections[chunk_index] = 0
      for hkl in chunks[chunk_index]:
        if len(chunks[chunk_index][hkl]) > 0:
          number_of_allocated_reflections[chunk_index] += len(chunks[chunk_index][hkl])
      print("\nRank %d chunk %d allocated reflections: %d"%(self.rank, chunk_index, number_of_allocated_reflections[chunk_index]))
    '''


  # Given a list of reflections, calculate their intensity statistics. Use case: a list of symmetry-equivalent reflections
  def calc_reflection_intensity_stats(self, reflections):

    if( len(reflections) == 0 ):
      return dict()


    # TODO: use library methods for the statistics


    count = 0
    average_intensity = 0.0
    esd = 0.0
    rmsd = 0.0

    for ref in reflections:
      #hkl = ref.get('miller_index_asymmetric')
      average_intensity   += ref.get('intensity.sum.value')
      esd                 += ref.get('intensity.sum.variance')
      count               += 1

    average_intensity /= count
    esd /= count

    for ref in reflections:
      rmsd += (ref.get('intensity.sum.value') - average_intensity) ** 2

    rmsd /= count
    rmsd = rmsd ** 0.5

    return {'average' : average_intensity,
            'count' : count,
            'esd' : esd,
            'rmsd': rmsd}

  # create a table to store timings of all steps
  def create_timing_table(self):
    self.timing_table = dict({'TOTAL':       [0, 0.0, 0.0],
                                                 'SETUP':       [1, 0.0, 0.0],
                                                 'BROADCAST':   [2, 0.0, 0.0],
                                                 'LOAD':        [3, 0.0, 0.0],
                                                 'DISTRIBUTE':  [4, 0.0, 0.0],
                                                 'ALL-TO-ALL':  [5, 0.0, 0.0],
                                                 'CONSOLIDATE': [6, 0.0, 0.0],
                                                 'AVERAGE':     [7, 0.0, 0.0],
                                                 'GATHER':      [8, 0.0, 0.0],
                                                 'MERGE':       [9, 0.0, 0.0]})

  # log time for a particular step
  def log_step_time(self, comm, step, finished=False):

    if not self.do_timing:
      return

    if not finished: # started - cache the start time
      self.timing_table[step][1] = MPI.Wtime()
    else: # finished - calculate the spent time
      self.timing_table[step][2] = MPI.Wtime() - self.timing_table[step][1]

  def save_log(self):
    if not self.do_timing: # for now the log is used for timimg only
      return

    filename = "timing_" + str(self.rank_count) + "_" + str(self.rank) + ".out"
    log_file = open(filename,'w')

    for step, value in sorted(self.timing_table.iteritems(), key=lambda (k,v): (v,k)):
      log_file.write("RANK %d %s: %f s\n"%(self.rank, step, value[2]))

    log_file.close()


  def run(self, comm, timing=True):

    self.do_timing = timing
    if self.do_timing:
      self.create_timing_table()

    self.rank = comm.Get_rank()
    self.rank_count = comm.Get_size()

    self.log_step_time(comm, "TOTAL")

    if( self.rank == 0 ):
      self.initialize()
      self.validate()

      ################
      # GET FILE LIST
      self.log_step_time(comm, "SETUP")
      file_list = self.get_list()
      print("rank 0 got a list of %d items"%len(file_list))

      #############################
      # GENERATE ALL MILLER INDICES
      full_miller_set = self.generate_all_miller_indices()
      print("\nGenerated %d miller indices"%full_miller_set.indices().size())

      ###################################################################
      # SPLIT MILLER INDICES INTO CHUNKS TO BE DISTRIBUTED OVER THE RANKS
      import numpy as np
      split_set = np.array_split(full_miller_set.indices(), self.rank_count)

      # within chunks, create an empty container for every hkl
      list_of_chunks_of_hkl_dictionaries = []
      for chunk in split_set:
        chunk_of_hkl_dictionaries = dict()
        for hkl_array in chunk:
          hkl = (hkl_array[0], hkl_array[1], hkl_array[2])
          chunk_of_hkl_dictionaries[hkl] = list()
        list_of_chunks_of_hkl_dictionaries.append(chunk_of_hkl_dictionaries)

      print("\nGenerated %d chunks of hkl containers"%len(list_of_chunks_of_hkl_dictionaries))

      #self.params.input.path = None # the input is already parsed
      self.log_step_time(comm, "SETUP", True)

      #########################################
      # BROADCAST WORK FROM RANK 0 TO ALL RANKS
      print ('\nTransmitting file list of length %d'%(len(file_list)))

      transmitted = dict(params = self.params, options = self.options, file_list = file_list, list_of_chunks_of_hkl_dictionaries = list_of_chunks_of_hkl_dictionaries)

    else:
      transmitted = None

    comm.barrier()
    self.log_step_time(comm, "BROADCAST")
    transmitted = comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['file_list'][self.rank::self.rank_count]
    chunks = transmitted['list_of_chunks_of_hkl_dictionaries']
    print ("\nRank %d received a file list of %d json-pickle file pairs and %d chunks of hkl containers" % (self.rank, len(new_file_list), len(chunks)))
    comm.barrier()
    self.log_step_time(comm, "BROADCAST", True)

    if( len(new_file_list) > 0 ):
      print ("\nRank %d: first file to load is: %s" % (self.rank, str(new_file_list[0])))

      ######################
      # EACH RANK: LOAD DATA
      comm.barrier()
      self.log_step_time(comm, "LOAD")
      experiments, reflections = self.load_data(new_file_list)
      print ('\nRank %d has read %d experiments consisting of %d reflections'%(self.rank, len(experiments), len(reflections)))
      comm.barrier()
      self.log_step_time(comm, "LOAD", True)

      #######################################################
      # EACH RANK: DISTRIBUTE REFLECTIONS OVER ALL HKL CHUNKS
      self.log_step_time(comm, "DISTRIBUTE")
      for chunk_index in range(0, len(chunks)):
        min_hkl = str(min(chunks[chunk_index]))
        max_hkl = str(max(chunks[chunk_index]))
        print("\nRank %d received chunk: %d: HKL count: %d; Min HKL: %s; Max HKL: %s" % (self.rank, chunk_index, len(chunks[chunk_index]), min_hkl, max_hkl))
      self.distribute_reflections_over_hkl_chunks(reflections=reflections, chunks=chunks)
      comm.barrier()
      self.log_step_time(comm, "DISTRIBUTE", True)
    else:
      print ("\nRank %d received no data" % self.rank)
      comm.barrier()
      comm.barrier()
      comm.barrier()

    #######################################################################
    # EACH RANK: RECEIVE ALL HKL CHUNKS (WITH THE SAME HKLS) FROM ALL RANKS
    self.log_step_time(comm, "ALL-TO-ALL")
    print("\nRank %d executing MPI all-to-all...\n"%self.rank)
    received_chunks = comm.alltoall(chunks)
    received_chunks_count = len(received_chunks)
    print ("\nAfter all-to-all rank %d received %d chunks of hkl containers" % (self.rank, received_chunks_count) )
    for chunk_index in range(0,received_chunks_count):
      min_hkl = str(min(received_chunks[chunk_index]))
      max_hkl = str(max(received_chunks[chunk_index]))
      print("\nAfter all-to-all rank %d received chunk %d: HKL count: %d; Min HKL: %s; Max HKL: %s" % (self.rank, chunk_index, len(received_chunks[chunk_index]), min_hkl, max_hkl))
    comm.barrier()
    self.log_step_time(comm, "ALL-TO-ALL", True)

    ######################################################################
    # EACH RANK: CONSOLIDATE ALL REFLECTION LISTS FROM ALL RECEIVED CHUNKS
    self.log_step_time(comm, "CONSOLIDATE")
    print("\nRank %d consolidating reflection lists...\n"%self.rank)
    for hkl in received_chunks[0]:
      for chunk_index in range(1,received_chunks_count):
        received_chunks[0][hkl] += received_chunks[chunk_index][hkl]
    comm.barrier()
    self.log_step_time(comm, "CONSOLIDATE", True)

    ####################################################
    # EACH RANK: DO STATISTICS ON REFLECTION INTENSITIES
    self.log_step_time(comm, "AVERAGE")
    print("\nRank %d doing intensity statistics...\n"%self.rank)
    all_rank_merged_reflections = merging_reflection_table()

    for hkl in received_chunks[0]:
      if( len(received_chunks[0][hkl]) > 0 ):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=received_chunks[0][hkl])
        all_rank_merged_reflections.append({'miller_index': hkl,
                                            'intensity': intensity_stats['average'],
                                            'multiplicity': intensity_stats['count'],
                                            'esd' : intensity_stats['esd'],
                                            'rmsd' : intensity_stats['rmsd']})

    print ("\nr\Rank %d merged %d reflection intensities"%(self.rank, all_rank_merged_reflections.size()))
    comm.barrier()
    self.log_step_time(comm, "AVERAGE", True)

    #############################################
    # EACH RANK: SEND REFLECTION TABLES TO RANK 0
    self.log_step_time(comm, "GATHER")
    if self.rank != 0:
      print("\nRank %d executing MPI gathering of all reflection tables at rank 0..."%self.rank)
    all_merged_reflection_tables = comm.gather(all_rank_merged_reflections, root = 0)
    comm.barrier()
    self.log_step_time(comm, "GATHER", True)

    ####################################
    # RANK 0: DO FINAL MERGING OF TABLES
    if self.rank == 0:
      self.log_step_time(comm, "MERGE")
      print ("\nRank 0 doing final merging of reflection tables received from all ranks...")
      final_merged_reflection_table = merging_reflection_table()
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      print("Rank 0 total merged reflections: {}".format(final_merged_reflection_table.size()) )
      self.log_step_time(comm, "MERGE", True)

    comm.barrier()
    self.log_step_time(comm, "TOTAL", True)

    self.save_log()

    return

if __name__ == '__main__':

  from mpi4py import MPI

  comm = MPI.COMM_WORLD

  script = Script()

  result = script.run(comm=comm)
  if result is None:
    sys.exit(1)

  print ("OK")

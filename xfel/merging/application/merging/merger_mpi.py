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

    return None

  # For each reflection find an hkl chunk, containing the reflection's hkl; append that reflection to the list of reflections for that hkl
  # Here we are assuming the input reflection list is sorted on miller indices
  def distribute_reflections_over_hkl_chunks(self, rank, reflections, chunks):

    # initializations
    hkl_cur = None # current hkl
    chunk_cur = None # chunk containing current hkl
    index_min = -1 # the first index in the sorted reflection list, where we encounter an hkl
    index_max = -1 # the (last + 1) index where we encounter the hkl

    #from IPython import embed; embed()

    total_distributed_reflections = 0
    for i, ref in enumerate(reflections):
      hkl = ref.get('miller_index_asymmetric')
      if( hkl_cur != hkl ): # encountered a new hkl
        if( hkl_cur != None ): # if there exists a previous hkl and a chunk for it, add all of that hkl's reflections to the corresponding chunk
          assert(index_min >= 0)
          if( chunk_cur != None ):
            index_max = i
            #TODO: to avoid reflection list copying, use extend() instead of append()
            #chunk_cur[hkl].extend(reflections[index_min:index_max])
            for j in range(index_min, index_max):
              chunk_cur[hkl_cur].append(reflections[j])
            total_distributed_reflections += (index_max - index_min)
          else:
            #print("\nHKL: %s is not found in any of the chunks"%str(hkl))
            pass

        # update current hkl, its starting index, and its chunk id
        hkl_cur = hkl
        index_min = i
        chunk_cur = self.get_hkl_chunk(hkl=hkl, chunks=chunks)

    print("Rank %d managed to distribute %d of %d reflections"%(rank, total_distributed_reflections, len(reflections)))

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

  def run(self, comm):

    rank = comm.Get_rank()
    rank_count = comm.Get_size()

    if( rank == 0 ):
      ################
      # GET FILE LIST
      self.initialize()
      self.validate()
      file_list = self.get_list()
      print("rank 0 got a list of %d items"%len(file_list))

      #############################
      # GENERATE ALL MILLER INDICES
      full_miller_set = self.generate_all_miller_indices()
      print("\nGenerated %d miller indices"%full_miller_set.indices().size())

      ###################################################################
      # SPLIT MILLER INDICES INTO CHUNKS TO BE DISTRIBUTED OVER THE RANKS
      import numpy as np
      split_set = np.array_split(full_miller_set.indices(), rank_count)

      # within chunks, create an empty container for every hkl
      list_of_chunks_of_hkl_dictionaries = []
      for chunk in split_set:
        chunk_of_hkl_dictionaries = dict()
        for hkl_array in chunk:
          hkl = (hkl_array[0], hkl_array[1], hkl_array[2])
          chunk_of_hkl_dictionaries[hkl] = list()
        list_of_chunks_of_hkl_dictionaries.append(chunk_of_hkl_dictionaries)

      print("\nGenerated %d chunks of hkl containers"%len(list_of_chunks_of_hkl_dictionaries))

      #from IPython import embed; embed()

      #self.params.input.path = None # the input is already parsed

      #########################################
      # BROADCAST WORK FROM RANK 0 TO ALL RANKS
      print ('\nTransmitting file list of length %d'%(len(file_list)))

      transmitted = dict(params = self.params, options = self.options, file_list = file_list, list_of_chunks_of_hkl_dictionaries = list_of_chunks_of_hkl_dictionaries)

    else:
      transmitted = None

    transmitted = comm.bcast(transmitted, root = 0)

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['file_list'][rank::rank_count]
    chunks = transmitted['list_of_chunks_of_hkl_dictionaries']
    #from IPython import embed; embed()
    print ("\nRank %d received a file list of %d json-pickle file pairs and %d chunks of hkl containers" % (rank, len(new_file_list), len(chunks)))

    if( len(new_file_list) > 0 ):
      print ("\nRank %d: first file to load is: %s" % (rank, str(new_file_list[0])))

      ######################
      # EACH RANK: LOAD DATA
      experiments, reflections = self.load_data(new_file_list)
      print ('\nRank %d has read %d experiments consisting of %d reflections'%(rank, len(experiments), len(reflections)))

      #######################################################
      # EACH RANK: DISTRIBUTE REFLECTIONS OVER ALL HKL CHUNKS
      for chunk_index in range(0, len(chunks)):
        min_hkl = str(min(chunks[chunk_index]))
        max_hkl = str(max(chunks[chunk_index]))
        print("\nRank %d received chunk: %d: HKL count: %d; Min HKL: %s; Max HKL: %s" % (rank, chunk_index, len(chunks[chunk_index]), min_hkl, max_hkl))

      self.distribute_reflections_over_hkl_chunks(rank=rank, reflections=reflections, chunks=chunks)
    else:
      print ("\nRank %d received no data" % rank)

    #######################################################################
    # EACH RANK: RECEIVE ALL HKL CHUNKS (WITH THE SAME HKLS) FROM ALL RANKS
    print("\nRank %d executing MPI all-to-all...\n"%rank)
    received_chunks = comm.alltoall(chunks)
    received_chunks_count = len(received_chunks)
    print ("\nAfter all-to-all rank %d received %d chunks of hkl containers" % (rank, received_chunks_count) )
    for chunk_index in range(0,received_chunks_count):
      min_hkl = str(min(received_chunks[chunk_index]))
      max_hkl = str(max(received_chunks[chunk_index]))
      print("\nAfter all-to-all rank %d received chunk %d: HKL count: %d; Min HKL: %s; Max HKL: %s" % (rank, chunk_index, len(received_chunks[chunk_index]), min_hkl, max_hkl))

    ######################################################################
    # EACH RANK: CONSOLIDATE ALL REFLECTION LISTS FROM ALL RECEIVED CHUNKS
    print("\nRank %d consolidating reflection lists...\n"%rank)
    for hkl in received_chunks[0]:
      for chunk_index in range(1,received_chunks_count):
        received_chunks[0][hkl] += received_chunks[chunk_index][hkl]

    #########################################
    # EACH RANK: MERGE REFLECTION INTENSITIES
    print("\nRank %d doing intensity statistics...\n"%rank)
    all_rank_merged_reflections = merging_reflection_table()

    #from IPython import embed; embed()


    for hkl in received_chunks[0]:
      if( len(received_chunks[0][hkl]) > 0 ):
        intensity_stats = self.calc_reflection_intensity_stats(reflections=received_chunks[0][hkl])
        all_rank_merged_reflections.append({'miller_index': hkl,
                                            'intensity': intensity_stats['average'],
                                            'multiplicity': intensity_stats['count'],
                                            'esd' : intensity_stats['esd'],
                                            'rmsd' : intensity_stats['rmsd']})

    #############################################
    # EACH RANK: SEND REFLECTION TABLES TO RANK 0
    if rank != 0:
      print("\nRank %d executing MPI gathering of all reflection tables at rank 0..."%rank)
    all_merged_reflection_tables = comm.gather(all_rank_merged_reflections, root = 0)
    #from IPython import embed; embed()

    ####################################
    # RANK 0: DO FINAL MERGING OF TABLES
    if rank == 0:
      print ("\nRank 0 doing final merging of reflection tables received from all ranks...")
      final_merged_reflection_table = merging_reflection_table()
      for table in all_merged_reflection_tables:
        final_merged_reflection_table.extend(table)
      print("Rank 0 total merged reflections: {}".format(final_merged_reflection_table.size()) )

    #from IPython import embed; embed()

    return

if __name__ == '__main__':

  from mpi4py import MPI

  comm = MPI.COMM_WORLD

  script = Script()

  result = script.run(comm=comm)
  if result is None:
    sys.exit(1)

  print ("OK")

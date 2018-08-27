from __future__ import division
import sys

"""
Functions for merging reflections
"""

class hkl_intensity_merger(object):
  def __init__(self, params, reflections):
    self.params = params
    self.reflections = reflections
 
  def merge(self, comm):
  
    rank = comm.Get_rank()
  
    count = {}
    average_intensity = {}
    esd = {}
    rmsd = {}
    hkl_cur = (0,0,0) # current hkl

    #from IPython import embed; embed()	

    for i, ref in enumerate(self.reflections):
      hkl = ref.get('miller_index_asymmetric')
      #print('\nhkl: %s'% str(hkl))      

      if( hkl_cur != hkl ): # encountered a new hkl
        
        # verify that the new hkl is not present in the dictionary
        if( hkl in average_intensity ):
          print ("HKL %s was encountered previously" % str(hkl))
          
        if(hkl_cur != {0,0,0} and hkl < hkl_cur):
          print("HKL out of order: hkl_cur=%s; hkl=%s"%(hkl_cur,hkl))
        
        hkl_cur = hkl # update current hkl
        average_intensity[hkl]     = ref.get('intensity.sum.value')
        esd[hkl]                   = ref.get('intensity.sum.variance')
        rmsd[hkl]                  = 0
        count[hkl]                 = 1

      else: # hkl hasn't changed - keep averaging for the current hkl
        average_intensity[hkl]    += ref.get('intensity.sum.value')
        esd[hkl]                  += ref.get('intensity.sum.variance')
        count[hkl]                += 1

    for hkl in average_intensity:
        average_intensity[hkl] /= count[hkl]
        esd[hkl] /= count[hkl]

    for i, ref in enumerate(self.reflections):
        hkl = ref.get('miller_index_asymmetric')
        rmsd[hkl] += (ref.get('intensity.sum.value') - average_intensity[hkl]) ** 2

    for hkl in rmsd:
        rmsd[hkl] /= count[hkl]
        rmsd[hkl] = rmsd[hkl] ** 0.5
        #if(rmsd[hkl]==0):
        #  print('\nZero rmsd: hkl: %s; count: %d'% (str(hkl),count[hkl]) )     
 
    #print ('Symmetry-independent relections intensities - lenghts of output arrays: %d; esd: %d; rmsd: %d'%(len(average_intensity), len(esd), len(rmsd)))
    
    print("\nRank %d has merged a list of %d reflections into a list of %d reflections; Multiplicity (min,max): (%d,%d); Intensity (min,max): (%d,%d); ESD (min,max): (%d,%d); RMSD (min,max): (%d,%d)"
      % (rank, len(self.reflections), len(average_intensity), min(count.values()), max(count.values()),\
      min(average_intensity.values()), max(average_intensity.values()), min(esd.values()), max(esd.values()), min(rmsd.values()), max(rmsd.values())))

    #print ('Minimum and maximum reflection multiplicity: %d, %d'% (min(count.values()), max(count.values())) )
    #print ('Minimum and maximum reflection intensity: %d, %d'% (min(average_intensity.values()), max(average_intensity.values())) )
    #print ('Minimum and maximum reflection esd: %d, %d'% (min(esd.values()), max(esd.values())) )
    #print ('Minimum and maximum reflection rmsd: %d, %d'% (min(rmsd.values()), max(rmsd.values())) )

    return

from xfel.merging.application.input.file_loader_mpi import Script as Script_Base
from xfel.merging.application.input.file_loader import file_loader
from xfel.merging.application.phil.phil import master_phil

class Script(Script_Base):
  '''A class for running the script.'''

  def merge_data(self, comm, reflections):
    merger = hkl_intensity_merger(self.params, reflections)
    merger.merge(comm = comm)
    return
  
  # given unit cell and space group info, generate all miller indexes
  def generate_all_miller_indices(self):
    
    #from IPython import embed; embed()
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
  
  # for each reflection find an hkl chunk, containing the reflection's hkl; append that reflection to the list of reflections for that hkl
  def distribute_reflections_over_hkl_chunks(self, rank, reflections, chunks):
    hkl_cur = (0,0,0) # current hkl
    chunk_cur = None # chunk containing current hkl 

    #from IPython import embed; embed()	
    
    total_distributed_reflections = 0
    for i, ref in enumerate(reflections):
      hkl = ref.get('miller_index_asymmetric')
      #print('\nhkl: %s'% str(hkl))      

      if( hkl_cur != hkl ): # encountered a new hkl
                
        hkl_cur = hkl # update current hkl and its chunk     
        chunk_cur = self.get_hkl_chunk(hkl=hkl, chunks=chunks)
      
      # append a reflection to the list of reflections for the current hkl 
      if( chunk_cur != None ):
        chunk_cur[hkl].append(ref)
        total_distributed_reflections += 1
    
    print("Rank %d managed to distribute %d of %d reflections"%(rank, total_distributed_reflections, len(reflections)))
    
    def average_reflection_intensities(self, reflections):
         count = {}
    average_intensity = {}
    esd = {}
    rmsd = {}
    hkl_cur = (0,0,0) # current hkl

    #from IPython import embed; embed()	

    for i, ref in enumerate(self.reflections):
      hkl = ref.get('miller_index_asymmetric')
      #print('\nhkl: %s'% str(hkl))      

      if( hkl_cur != hkl ): # encountered a new hkl
        
        # verify that the new hkl is not present in the dictionary
        if( hkl in average_intensity ):
          print ("HKL %s was encountered previously" % str(hkl))
          
        if(hkl_cur != {0,0,0} and hkl < hkl_cur):
          print("HKL out of order: hkl_cur=%s; hkl=%s"%(hkl_cur,hkl))
        
        hkl_cur = hkl # update current hkl
        average_intensity[hkl]     = ref.get('intensity.sum.value')
        esd[hkl]                   = ref.get('intensity.sum.variance')
        rmsd[hkl]                  = 0
        count[hkl]                 = 1

      else: # hkl hasn't changed - keep averaging for the current hkl
        average_intensity[hkl]    += ref.get('intensity.sum.value')
        esd[hkl]                  += ref.get('intensity.sum.variance')
        count[hkl]                += 1

    for hkl in average_intensity:
        average_intensity[hkl] /= count[hkl]
        esd[hkl] /= count[hkl]

    for i, ref in enumerate(self.reflections):
        hkl = ref.get('miller_index_asymmetric')
        rmsd[hkl] += (ref.get('intensity.sum.value') - average_intensity[hkl]) ** 2

    for hkl in rmsd:
        rmsd[hkl] /= count[hkl]
        rmsd[hkl] = rmsd[hkl] ** 0.5
        #if(rmsd[hkl]==0):
        #  print('\nZero rmsd: hkl: %s; count: %d'% (str(hkl),count[hkl]) )     

    

    
  def run(self, comm):
     
    rank = comm.Get_rank()
    rank_count = comm.Get_size()

    if( rank == 0 ):
      self.initialize()
      self.validate()
      file_list = self.get_list()
      
      print("rank 0 got a list of %d items"%len(file_list))
      
      ####################################################
      # GENERATE ALL MILLER INDEXES
      full_miller_set = self.generate_all_miller_indices()
      ####################################################
      
      #################################################################
      # SPLIT MILLER INDEXES INTO CHUNKS
      import numpy as np
      split_set = np.array_split(full_miller_set.indices(), rank_count)
      #from IPython import embed; embed()
      #print("\nSplit Miller set into %d sets"%split_set.size())
      #################################################################
     
      list_of_chunks_of_hkl_dictionaries = []
      for chunk in split_set:
        chunk_of_hkl_dictionaries = {}
        for hkl_array in chunk:
          hkl = (hkl_array[0], hkl_array[1], hkl_array[2])
          chunk_of_hkl_dictionaries[hkl] = list()
        list_of_chunks_of_hkl_dictionaries.append(chunk_of_hkl_dictionaries)
      
      print("\nGenerated %d chunks of hkl dictionaries"%len(list_of_chunks_of_hkl_dictionaries))
     
      #from IPython import embed; embed()
     
      #self.params.input.path = None # the input is already parsed
      ################################################
      # BROADCASTING
      print ('Transmitting file list of length %d'%(len(file_list)))

      transmitted = dict(params = self.params, options = self.options, file_list = file_list, list_of_chunks_of_hkl_dictionaries = list_of_chunks_of_hkl_dictionaries)

    else:
      transmitted = None
    
    transmitted = comm.bcast(transmitted, root = 0)
    ################################################

    self.params = transmitted['params']
    self.options = transmitted['options']
    new_file_list = transmitted['file_list'][rank::rank_count]
    chunks = transmitted['list_of_chunks_of_hkl_dictionaries']
    #from IPython import embed; embed()
    print ("\nRank %d received a file list of %d jason-pickle file pairs as well as %d chunks of hkl dictionaries" % (rank, len(new_file_list), len(chunks)))   
	
    if( len(new_file_list) > 0 ):
      print ("\nRank %d: first file to load is: %s" % (rank, str(new_file_list[0])))
      
      #########################################################
      # EACH RANK LOADING DATA
      experiments, reflections = self.load_data(new_file_list)
      print ('\nRank %d has read %d experiments consisting of %d reflections'%(rank, len(experiments), len(reflections)))
      #########################################################
     
      #############################################################################
      # EACH RANK DISTRIBUTING ITS REFLECTIONS, OVER ALL HKL CHUNKS
      #############################################################################  
      
      for chunk_index in range(0, len(chunks)):
        min_hkl = str(min(chunks[chunk_index]))
        max_hkl = str(max(chunks[chunk_index]))
        print("Rank: %d; Chunk: %d; Size: %d; Min HKL: %s; Max HKL: %s" % (rank, chunk_index, len(chunks[chunk_index]), min_hkl, max_hkl))

      self.distribute_reflections_over_hkl_chunks(rank=rank, reflections=reflections, chunks=chunks)
      
      #######################################
      # ALL TO ALL     
      received_chunks = comm.alltoall(chunks)
      received_chunks_count = len(received_chunks)
      print ("\nAfter ALL-TO-ALL rank %d received %d chunks of hkl-intensity dictionaries" % (rank, received_chunks_count) )      
      for chunk_index in range(0,received_chunks_count):
        min_hkl = str(min(received_chunks[chunk_index]))
        max_hkl = str(max(received_chunks[chunk_index]))
        print("Rank: %d; Chunk: %d; Size: %d; Min HKL: %s; Max HKL: %s" % (rank, chunk_index, len(received_chunks[chunk_index]), min_hkl, max_hkl))
      ########################################
      
      ###################################################################
      # EACH RANK CONSOLIDATING ALL REFLECTION LISTS FROM RECEIVED CHUNKS
      for hkl in received_chunks[0]:
        for chunk_index in range(1,received_chunks_count):
          received_chunks[0][hkl] += received_chunks[chunk_index][hkl]
      ###################################################################
      
      ########################################################
      # EACH RANK MERGING
      for hkl in received_chunks[0]:
        merge_intensities(reflections=received_chunks[0][hkl])
      ########################################################
    else:
       print ("\nRank %d received no data" % rank)   

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

from __future__ import division
from scitbx.array_family import flex


"""The mapper_unmapper class solves the problem of selecting a subset of the data
   that is "visited" in the simulation.  In other words, selecting
   Miller indices that are actually measured, and frames that actually have observations.
   Otherwise, there would be I-parameters and G-parameters with undefined gradients.
"""

def mapper_factory(base_class):
  class mapper_unmapper(base_class):
    def __init__(self,Ibase,Gbase,I_visited,G_visited,FSIM,**kwargs):
      g_counter=0; forward_map_G=flex.size_t(len(G_visited)); backward_map_G=flex.size_t()
      for s in xrange(len(G_visited)):
        #print s, G_visited[s], c[len_I + s], c[len_I + len(Gbase) + s]
        if G_visited[s]:
          forward_map_G[s] = g_counter
          backward_map_G.append(s)
          g_counter+=1

      subsetGbase = Gbase.select(backward_map_G)
      remapped_frame = forward_map_G.select(FSIM.frame)

      i_counter=0; forward_map_I=flex.size_t(len(I_visited)); backward_map_I=flex.size_t()
      for s in xrange(len(I_visited)):
        #print s,I_visited[s], c[s]
        if I_visited[s]:
          forward_map_I[s] = i_counter
          backward_map_I.append(s)
          i_counter+=1
      subsetIbase = Ibase.select(backward_map_I)
      remapped_miller = forward_map_I.select(FSIM.miller)

      from cctbx.examples.merging import intensity_data
      remapped_FSIM = intensity_data()
      remapped_FSIM.raw_obs = FSIM.raw_obs
      remapped_FSIM.exp_var = FSIM.exp_var
      remapped_FSIM.stol_sq = FSIM.stol_sq
      remapped_FSIM.frame   = remapped_frame
      remapped_FSIM.miller  = remapped_miller

      base_class.__init__(self,subsetIbase,subsetGbase,remapped_FSIM,**kwargs)
      fitted_I,fitted_G,fitted_B = self.unpack()

      self.expanded_G = flex.double(len(Gbase))
      self.expanded_B = flex.double(len(Gbase))
      for s in xrange(len(G_visited)):
        if G_visited[s]:
          self.expanded_G[s]=fitted_G[ forward_map_G[s] ]
          self.expanded_B[s]=fitted_B[ forward_map_G[s] ]

      self.expanded_I = flex.double(len(Ibase))
      for s in xrange(len(I_visited)):
        if I_visited[s]:
          self.expanded_I[s]=fitted_I[ forward_map_I[s] ]

      print "DONE UNMAPPING HERE"

    def e_unpack(self): return (self.expanded_I, self.expanded_G,self.expanded_B)

  return mapper_unmapper

def mapper_factory_with_explicit_B(base_class):
  class mapper_unmapper(base_class):
    def __init__(self,Ibase,Gbase,Bbase,I_visited,G_visited,FSIM,**kwargs):
      g_counter=0; forward_map_G=flex.size_t(len(G_visited)); backward_map_G=flex.size_t()
      for s in xrange(len(G_visited)):
        #print s, G_visited[s], c[len_I + s], c[len_I + len(Gbase) + s]
        if G_visited[s]:
          forward_map_G[s] = g_counter
          backward_map_G.append(s)
          g_counter+=1

      subsetGbase = Gbase.select(backward_map_G)
      subsetBbase = Bbase.select(backward_map_G)
      remapped_frame = forward_map_G.select(FSIM.frame)

      i_counter=0; forward_map_I=flex.size_t(len(I_visited)); backward_map_I=flex.size_t()
      for s in xrange(len(I_visited)):
        #print s,I_visited[s], c[s]
        if I_visited[s]:
          forward_map_I[s] = i_counter
          backward_map_I.append(s)
          i_counter+=1
      subsetIbase = Ibase.select(backward_map_I)
      remapped_miller = forward_map_I.select(FSIM.miller)

      from cctbx.examples.merging import intensity_data
      remapped_FSIM = intensity_data()
      remapped_FSIM.raw_obs = FSIM.raw_obs
      remapped_FSIM.exp_var = FSIM.exp_var
      remapped_FSIM.stol_sq = FSIM.stol_sq
      remapped_FSIM.frame   = remapped_frame
      remapped_FSIM.miller  = remapped_miller

      base_class.__init__(self,subsetIbase,subsetGbase,subsetBbase,remapped_FSIM,**kwargs)
      fitted_I,fitted_G,fitted_B = self.unpack()

      self.expanded_G = flex.double(len(Gbase))
      self.expanded_B = flex.double(len(Gbase))
      for s in xrange(len(G_visited)):
        if G_visited[s]:
          self.expanded_G[s]=fitted_G[ forward_map_G[s] ]
          self.expanded_B[s]=fitted_B[ forward_map_G[s] ]

      self.expanded_I = flex.double(len(Ibase))
      for s in xrange(len(I_visited)):
        if I_visited[s]:
          self.expanded_I[s]=fitted_I[ forward_map_I[s] ]

      print "DONE UNMAPPING HERE"

    def e_unpack(self): return (self.expanded_I, self.expanded_G,self.expanded_B)

  return mapper_unmapper

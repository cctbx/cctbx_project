from __future__ import division
from __future__ import print_function

class small_cell_orientation:
 """ Class for determining an orientation matrix given a set of reflections """

 def __init__(self,miller_indices, u_vectors, sym):
  """
  @param miller_indices indexed miller indices
  @param u_vectors reciprocal space vectors for the given miller indices
  @sym cctbx symmetry object
  """
  from libtbx import adopt_init_args
  adopt_init_args(self,locals())
  # assumes miller_indices and u-vectors are flex.vec3_doubles
  # sym is a cctbx.crystal.symmetry

 def unrestrained_setting(self):
  """ Calculate the basis vectors from N spots using numpy. This is equation 5 in
  Brewster et. al 2015.
  @return a cctbx crystal_orientation object
  """
  from dials.array_family import flex
  from scitbx.matrix import sqr
  import numpy as np

  if len(self.miller_indices) < 3:
    return None

  N = self.miller_indices.size()

  hkl = self.miller_indices.as_double()
  hkl.reshape(flex.grid((N,3)))
  hkl = hkl.as_numpy_array()

  xyz = self.u_vectors.as_double()
  xyz.reshape(flex.grid((N,3)))
  xyz = xyz.as_numpy_array()

  try:
    result = np.linalg.lstsq(hkl,xyz)
  except Exception as e:
    print("Exception while calculating basis vectors: %s"%e.message)
    return None

  solution,self.residuals,rank,singular = result[0],result[1],result[2],result[3]
  if len(self.residuals) == 0:
    self.residuals = [0,0,0] # happens when only 3 spots in the max clique

  print("Summed squared residuals of x,y,z for %d spots in 1/angstroms: %.7f, %.7f, %.7f"%(N,self.residuals[0],self.residuals[1],self.residuals[2]))
  Amatrix = sqr(solution.flatten()).transpose()

  from cctbx import crystal_orientation
  ori = crystal_orientation.crystal_orientation(Amatrix, crystal_orientation.basis_type.reciprocal)
  return ori

if __name__=="__main__":
  import libtbx.easy_pickle as ep
  data = ep.load("r0013_shot-s00-20130311223605649.example")
  miller_indices,u_vectors,symmetry = data[0],data[1],data[2]
  S = small_cell_orientation(miller_indices, u_vectors, symmetry)
  print(S.unrestrained_setting())

from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from rstbx.array_family import flex
from rstbx.dps_core import Direction
from annlib_ext import AnnAdaptor
from rstbx.dps_core import sampling
from rstbx.dps_core import Direction, directional_show
#from libtbx.development.timers import Timer

# various sampling algorithms to cover the directional hemisphere

class SimpleSampler(sampling.SimpleSamplerTool):
  # replaces the setSampling_internal() algorithm from dptbx.cpp, with
  # no confinement directions

  def __init__(self,maxgrid,charactergrid):
    self.max_grid = maxgrid  # the maximum allowable grid, corresponding
                             # to the grid sampling used in the Rossman DPS
                             # paper; approx 7900 directions; 0.4 seconds

    self.characteristic_grid = charactergrid # the natural grid sampling
                             # reflective of the problem at hand =
                             # obs. resolution limit / most conservative cell

    sampling.SimpleSamplerTool.__init__(self,
      min( self.max_grid, self.characteristic_grid ))
      # initial directional sampling in radians

  def get_top_solutions(self,ai,input_directions,size,cutoff_divisor):
    # size was 30 in the Rossmann DPS paper
    kval_cutoff = ai.getXyzSize()/cutoff_divisor;

    hemisphere_solutions = flex.Direction();
    hemisphere_solutions.reserve(size);

    for i in range(len(input_directions)):
      sampled_direction = ai.fft_result(input_directions[i])
      if sampled_direction.kval > kval_cutoff:
        hemisphere_solutions.append(sampled_direction)

    if (hemisphere_solutions.size()<3):
      return hemisphere_solutions
    kvals = flex.double([
  hemisphere_solutions[x].kval for x in range(len(hemisphere_solutions))])

    perm = flex.sort_permutation(kvals,True)

    #  conventional algorithm; just take the top scoring hits.
    #  use this for quick_grid, when it is known ahead of time
    #  that many (hundreds) of directions will be retained

    hemisphere_solutions_sort = flex.Direction(
        [hemisphere_solutions[p] for p in perm[0:min(size,len(perm))]])

    return hemisphere_solutions_sort;

  def hemisphere(self,ai,size,cutoff_divisor):
    #this replaces the hemisphere() function in dptbx.cpp
    return self.get_top_solutions(ai,self.angles,size,cutoff_divisor,
                                  grid=self.incr)

class HemisphereSampler(SimpleSampler):
  def __init__(self,max_grid,characteristic_grid,quick_grid):
    self.quick_grid = quick_grid
    SimpleSampler.__init__(self,max_grid,characteristic_grid)

    if self.quick_grid >= self.incr:
      #some quick benchmarks:
      estimate_number_of_target_grids = 2.*math.pi / (self.incr*self.incr)

      estimate_number_of_quick_grids = 2.*math.pi / (self.quick_grid*self.quick_grid)

      self.second_round_sampling = (
        math.pi*self.quick_grid*self.quick_grid)/(self.incr*self.incr)

      self.target_size_for_quick_grid_answer_vector = int(
        estimate_number_of_quick_grids /2./self.second_round_sampling)

  def get_top_solutions(self,ai,input_directions,size,cutoff_divisor,grid):
    # size was 30 in the Rossmann DPS paper

    kval_cutoff = ai.getXyzSize()/cutoff_divisor;

    hemisphere_solutions = flex.Direction();
    hemisphere_solutions.reserve(size);
    for i in range(len(input_directions)):
      D = sampled_direction = ai.fft_result(input_directions[i])

      if D.real < self.max_cell_input and sampled_direction.kval > kval_cutoff:
        #directional_show(D, "ddd %5d"%i)
        hemisphere_solutions.append(sampled_direction)

    if (hemisphere_solutions.size()<3):
      return hemisphere_solutions
    kvals = flex.double([
  hemisphere_solutions[x].kval for x in range(len(hemisphere_solutions))])

    perm = flex.sort_permutation(kvals,True)

    #  need to be more clever than just taking the top 30.
    #  one huge cluster around a strong basis direction could dominate the
    #  whole hemisphere map, preventing the discovery of three basis vectors

    perm_idx = 0
    unique_clusters = 0
    hemisphere_solutions_sort = flex.Direction()
    while perm_idx < len(perm) and \
          unique_clusters < size:
      test_item = hemisphere_solutions[perm[perm_idx]]
      direction_ok = True
      for list_item in hemisphere_solutions_sort:
        distance = math.sqrt(math.pow(list_item.dvec[0]-test_item.dvec[0],2) +
                     math.pow(list_item.dvec[1]-test_item.dvec[1],2) +
                     math.pow(list_item.dvec[2]-test_item.dvec[2],2)
                     )
        if distance < 0.087: #i.e., 5 degrees
           direction_ok=False
           break
      if direction_ok:
        unique_clusters+=1
        hemisphere_solutions_sort.append(test_item)
      perm_idx+=1

    return hemisphere_solutions_sort;

  def hemisphere(self,ai,max_cell,size = 30,cutoff_divisor = 4.):
    self.max_cell_input = max_cell
    if self.quick_grid >= self.incr: #trigger adaptive sampling
      quick_directions = SimpleSampler.construct_hemisphere_grid(
        self,self.quick_grid)
      quick_sampling = SimpleSampler.get_top_solutions(
        self,ai,quick_directions,
        self.target_size_for_quick_grid_answer_vector,
        cutoff_divisor = 10.)
      self.angles = SimpleSampler.construct_hemisphere_grid(
        self,self.incr)
      self.restrict_target_angles_to_quick_winners(quick_sampling)
      # now the hemisphere search will be restricted to the adaptive
      # sample rather than the entire hemisphere
    else:
      self.angles = SimpleSampler.construct_hemisphere_grid(
        self,self.incr)

    unrefined_basis_vectors = self.get_top_solutions(ai,self.angles,size,
      cutoff_divisor,grid=self.incr)

    #now refine them by the first method
    grid_refined_vectors = self.refine_top_solutions_by_grid_search(
      old_solutions = unrefined_basis_vectors, ai = ai)

    ai.setSolutions(grid_refined_vectors)

  def restrict_target_angles_to_quick_winners(self,quick_sampling):

    # find which target_grid directions are near neighbors to
    # the quick_sampling winners from the first round of testing.
    # in the second round, restrict the calculation of FFTs to those directions

    self.restricted = flex.Direction()
    target = flex.double()
    for iz in self.angles:
      # fairly inefficient in Python; expect big improvement in C++
      v = iz.dvec
      target.append(v[0]);target.append(v[1]);target.append(v[2]);
    kn = int(2 *self.second_round_sampling)

    #construct k-d tree for the reference set
    A = AnnAdaptor(data = target, dim = 3, k = kn)

    query = flex.double()
    for j in quick_sampling:
      v = j.dvec
      query.append(v[0]);query.append(v[1]); query.append(v[2]);
      if abs((math.pi/2.) - j.psi) < 0.0001: #take care of equatorial boundary
        query.append(-v[0]);query.append(-v[1]); query.append(-v[2]);

    A.query(query) #find nearest neighbors of query points
    neighbors = flex.sqrt(A.distances)
    neighborid = A.nn

    accept_flag = flex.bool(len(self.angles))
    for idx in range(len(neighbors)):
      # use small angle approximation to test if target is within desired radius
      if neighbors[idx] < self.quick_grid:
        accept_flag[neighborid[idx]] = True

    #go through all of the original target angles
    for iz in range(len(self.angles)):
      if accept_flag[iz]:
        self.restricted.append(self.angles[iz])

    self.angles = self.restricted

  def refine_top_solutions_by_grid_search(self, old_solutions, ai):
    new_solutions = flex.Direction()
    final_target_grid = self.characteristic_grid / 100.
    for item in old_solutions:
      new_solutions.append(ai.refine_direction(
        candidate = item,
        current_grid = self.incr,
        target_grid = final_target_grid))
    return new_solutions

def hemisphere_shortcut(ai,characteristic_sampling,max_cell):
    H = HemisphereSampler(
      max_grid = 0.029,
      characteristic_grid = characteristic_sampling,
      quick_grid = 0.016) # all grid parameters in radians
    H.hemisphere(ai,max_cell,size=30,cutoff_divisor=4.)
def hemisphere_refine_shortcut(ai,characteristic_sampling,unrefined):
    H = HemisphereSampler(
      max_grid = 0.029,
      characteristic_grid = characteristic_sampling,
      quick_grid = 0.016) # all grid parameters in radians
    cutoff_divisor=4.
    H.kval_cutoff = ai.getXyzSize()/cutoff_divisor;
    return H.refine_top_solutions_by_grid_search(unrefined,ai)


"""New adaptive algorithms to provide directional sampling of FFT's; July 2006.

1) The legacy algorithm is moved from the C++ function
   AutoIndexEngine::hemisphere() to the Python class SimpleSampler.  Moving
   to Python makes it simpler to rapidly prototype new methods.  Alogrithm
   provides uniform sampling over the hemisphere with a default spacing of
   0.029 radians as described in the Rossmann DPS paper.  An override
   "characteristic grid" is allowed in case the cell lengths are very large.

2) The Python implementation is somewhat slower than C++ (~20% factor) probably
   owing to the many Boost.Python interface calls.  It is the eventual
   intention to port everything back to C++ once the redesign is stable.

3) The new algorithms are embodied in the derived class HemisphereSampler,
   with the intention of balancing competing design requirements of fine
   sampling for large unit cells, and sparse enough sampling for ideal
   performance.  Diffraction cases are now divided into three categories,
   and are treated differently depending on the likely characteristic
   scale determined in the tnear2.py module.

   a) Problems with a characteristic scale > 0.029 radians (parameter max_grid)
      are still treated exactly according to Rossmann DPS, giving about 7900
      sampling directions.

   b) For problems with a characteristic scale < max_grid, the Simple sampling
      method is still used, but with a sampling fineness set to the
      characteristic scale; provided that the characteristic scale is still
      > 0.016 radians (parameter quick_grid).  At this limit there are about
      25000 sampling directions, using ~2 seconds of CPU when encoded in C++
      and when run on the fastest processors. Large-cell problems tend to take
      somewhat longer as the transforms are larger and have more data. This is
      the upper acceptable time elapsed for the sampling process.

   c) For still larger problems, a two-step process is employed. i) The
      limiting resolution of the dataset is degraded until the characteristic
      scale is quick_grid.  The Rossmann DPS sampling is then run on this
      data subset.  Promising direction vectors are selected.  ii) Now
      using the full dataset, the hemisphere is resampled at the original
      characteristic scale (for the full-dataset); but ONLY in areas immediately
      surrounding the promising directions from step i).  This two-step process
      is very efficient, and executes in <10 seconds even for the
      largest-cell problems examined.

4) Workaround for coarse gridding.  In the two-step sampling for large-cell
   problems, there is a danger that the quick_grid sampling of step i) will
   be too coarse. Coarse sampling may completely overlook a sharp peak
   coresponding to a true unit cell basis vector.  Two workarounds are employed.
   First, the cutoff for defining a direction vector of interest is relaxed
   in step i).  Instead of using kval > xyzdata.size()/4 (as in the legacy
   algorithm) we take kval > xyzdata.size()/10.  Secondly, we allow a very large
   number of directions of interest.  Instead of the 30 peaks used by Rossmann
   DPS, we allow up to several thousand.  The determining parameter is that
   when the interesting directions are finely sampled in step ii), we still want
   to keep the total number of FFTs under about 12500 (half of the quick_grid
   burden).  Therefore the quick_grid parameter ultimately limits the accessible
   problem size, but none of the crystallographic structures to date are limited
   by this implementation.

5) Workaround for unequal basis vector scores.  Even when sampled at a known
   basis vector direction, the kval score of the FFT is never 100% of the
   number of input data.  Sometimes the kval score for one or two of the
   basis vectors is significantly lower.  Moreover the FWHM can differ
   dramatically among the three basis vectors.  In extreme cases, one huge
   cluster around a strong basis direction can dominate the entire hemisphere
   map, precluding the discovery of the other two basis vectors.  In the
   Rossmann DPS paper it is stated simply that the top 30 sampled directions
   are taken, but this clearly needs to be modified in view of these results.
   Instead, we modify the algorithm to accept a maximum of 30 (parameter "size")
   unique clusters, with an unlimited number of directions within each cluster.
   By definition a cluster encompasses any top-scoring direction vector that
   is within 5-degrees of the first-identified (top-score) vector within that
   cluster.  This procedure gives the needed flexibility to find the difficult basis
   vectors.  After finding the top-scoring sampled directions, the kval scores
   are refined by a grid-search.  A legacy C++ function (AutoIndexEngine::
   setSolutions) then resorts the vectors, and picks the top vector in each
   cluster.  This is completely general since Niggli cell basis vectors are
   never within 5-degrees of each other.

6) Improvement in vector refinement.  In the legacy C++ procedure,
   AutoIndexEngine::hemisphere_refine_by_grid_search(), the directional vector
   is refined with nested grid searches until the granularity reaches 0.0001
   radians.  In the new Python procedure refine_top_solutions_by_grid_search(),
   the final granularity is customized to be (1/50) of the characteristic
   scale.  From the experience of publications/largecell1/fwhm.py, this
   granularity should be enough to correctly sample the kval peak, and in routine
   cases the value will be >0.0001 and thus more efficient.
"""

import math
from rstbx.array_family import flex
from rstbx.dps_core import Direction

# sampling algorithm to cover the directional hemisphere

class SimpleSamplerTool:

  def __init__(self,characteristic_grid):
    # the maximum allowable characteristic grid should be about 0.029 radians, corresponding
    # to the grid sampling used in the Rossman DPS paper; approx 7900 directions; 0.4 seconds

    # the characteristic grid sampling should be reflective of the problem at hand =
    # approximately the observed resolution limit / most conservative (largest) cell estimate

    self.incr = characteristic_grid
    # initial directional sampling in radians

    self.angles = self.construct_hemisphere_grid(self.incr)

  def construct_hemisphere_grid(self,sampling):
    #psi is the equivalent of latitude, measured as an angle from the North pole
    psi_index_range = int (0.5 + math.pi/2./sampling) #equvlnt of python round
    adjusted_psi_incr = math.pi/2./psi_index_range #adjusted for integral no.
    angles = flex.Direction()
    angles.reserve(4*psi_index_range*psi_index_range)

    for x in xrange(psi_index_range+1):
      psi = x * adjusted_psi_incr
      if psi > math.pi: eps = 1E-4; psi=math.pi-eps

      #phi is the equivalent of longitude
      if (psi==0):
        phi=0.;
        angles.append(Direction(psi=psi,phi=phi));
      else:
        phi_index_range = int (0.5 + 2.*math.pi*math.sin(psi)/sampling);
        adjusted_phi_incr = 2.*math.pi/phi_index_range;
        for y in xrange(phi_index_range):
          phi = y * adjusted_phi_incr;
          angles.append(Direction(psi=psi,phi=phi));

    return angles

class HemisphereSamplerBase(SimpleSamplerTool):
  def __init__(self,characteristic_grid,max_cell):
    SimpleSamplerTool.__init__(self,characteristic_grid)
    self.max_cell=max_cell
  def get_top_solutions(self,ai,input_directions,size,cutoff_divisor,grid):

    kval_cutoff = ai.getXyzSize()/cutoff_divisor;

    hemisphere_solutions = flex.Direction();
    hemisphere_solutions.reserve(size);
    print "# of input directions", len(input_directions)
    for i in xrange(len(input_directions)):
      D = sampled_direction = ai.fft_result(input_directions[i])

      #hardcoded parameter in for silicon example.  Not sure at this point how
      #robust the algorithm is when this value is relaxed.

      if D.real < self.max_cell and sampled_direction.kval > kval_cutoff:
        print i ,"%.4f %.2f %.2f kmax=%d kval=%f"%(D.real,
    180*D.psi/math.pi, 180.*D.phi/math.pi, D.kmax, D.kval)
        hemisphere_solutions.append(sampled_direction)
    if (hemisphere_solutions.size()<3):
      return hemisphere_solutions
    kvals = flex.double([
  hemisphere_solutions[x].kval for x in xrange(len(hemisphere_solutions))])

    perm = flex.sort_permutation(kvals,1)

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
        if distance < 0.087: #i.e., 5 degrees radius for clustering analysis
           direction_ok=False
           break
      if direction_ok:
        unique_clusters+=1
      hemisphere_solutions_sort.append(test_item)
      perm_idx+=1

    return hemisphere_solutions_sort;

  def hemisphere(self,ai,size = 30,cutoff_divisor = 4.,verbose=True):

    unrefined_basis_vectors = self.get_top_solutions(ai,self.angles,size,
      cutoff_divisor,grid=self.incr)

    if verbose:
      for i in xrange(len(unrefined_basis_vectors)):
        D = unrefined_basis_vectors[i];
        print i ,"%7.4f %8.2f %8.2f kmax=%2d kval=%5.1f"%(D.real,
            180*D.psi/math.pi, 180.*D.phi/math.pi, D.kmax, D.kval)

    ai.setSolutions(unrefined_basis_vectors)

def hemisphere_shortcut(ai,characteristic_sampling,max_cell):
    H = HemisphereSamplerBase(characteristic_grid = characteristic_sampling,max_cell=max_cell)
    H.hemisphere(ai,size=480,cutoff_divisor=1.5)
    #1.0/cutoff divisor is the height requirement (as a fraction of total Bragg spots)
    #required to consider a direction as a candidate basis vector.  4.0 good for
    #macromolecular work, 1.5 ok for small molecules; probably 1.1 would do.

from __future__ import absolute_import, division, print_function
from six.moves import range
import math
import os
from libtbx.test_utils import approx_equal
from rstbx.array_family import flex
from rstbx.dps_core import Direction, directional_show, SimpleSamplerTool

# sampling algorithm to cover the directional hemisphere

diagnostic=False

def test_simple_sampler():
  # SimpleSamplerTool migrated from Python to C++.
  SST = SimpleSamplerTool(0.014)
  SST.construct_hemisphere_grid(SST.incr)
  result = approx_equal(SST.incr, 0.014)
  result &= (len(SST.angles) == 32227)
  reference_psi = [0.0, 0.014024967203525862, 0.028049934407051724, 0.042074901610577586, 0.056099868814103448, 0.070124836017629311, 0.084149803221155173, 0.098174770424681035, 0.1121997376282069, 0.12622470483173276, 0.14024967203525862, 0.15427463923878448, 0.16829960644231035, 0.18232457364583621, 0.19634954084936207, 0.21037450805288793, 0.22439947525641379, 0.23842444245993966, 0.25244940966346552, 0.26647437686699138, 0.28049934407051724, 0.2945243112740431, 0.30854927847756897, 0.32257424568109483, 0.33659921288462069, 0.35062418008814655, 0.36464914729167242, 0.37867411449519828, 0.39269908169872414, 0.40672404890225, 0.42074901610577586, 0.43477398330930173, 0.44879895051282759, 0.46282391771635345, 0.47684888491987931, 0.49087385212340517, 0.50489881932693104, 0.5189237865304569, 0.53294875373398276, 0.54697372093750862, 0.56099868814103448, 0.57502365534456035, 0.58904862254808621, 0.60307358975161207, 0.61709855695513793, 0.6311235241586638, 0.64514849136218966, 0.65917345856571552, 0.67319842576924138, 0.68722339297276724, 0.70124836017629311, 0.71527332737981897, 0.72929829458334483, 0.74332326178687069, 0.75734822899039655, 0.77137319619392242, 0.78539816339744828, 0.79942313060097414, 0.8134480978045, 0.82747306500802587, 0.84149803221155173, 0.85552299941507759, 0.86954796661860345, 0.88357293382212931, 0.89759790102565518, 0.91162286822918104, 0.9256478354327069, 0.93967280263623276, 0.95369776983975862, 0.96772273704328449, 0.98174770424681035, 0.99577267145033621, 1.0097976386538621, 1.0238226058573878, 1.0378475730609138, 1.0518725402644398, 1.0658975074679655, 1.0799224746714913, 1.0939474418750172, 1.1079724090785432, 1.121997376282069, 1.1360223434855947, 1.1500473106891207, 1.1640722778926467, 1.1780972450961724, 1.1921222122996982, 1.2061471795032241, 1.2201721467067501, 1.2341971139102759, 1.2482220811138016, 1.2622470483173276, 1.2762720155208536, 1.2902969827243793, 1.3043219499279051, 1.318346917131431, 1.332371884334957, 1.3463968515384828, 1.3604218187420085, 1.3744467859455345, 1.3884717531490605, 1.4024967203525862, 1.416521687556112, 1.4305466547596379, 1.4445716219631639, 1.4585965891666897, 1.4726215563702154, 1.4866465235737414, 1.5006714907772674, 1.5146964579807931, 1.5287214251843189, 1.5427463923878448, 1.5567713595913708, 1.5707963267948966]
  all_psi = []
  for i in SST.angles:
    if i.psi not in all_psi:  all_psi.append(i.psi)
  result &= approx_equal(all_psi, reference_psi)
  with open(os.devnull, 'w') as devnull:
    one_slice_phi = [i.phi for i in SST.angles if approx_equal(i.psi, 0.014024967203525862, out=devnull)]
  reference_phi = [0.0, 1.0471975511965976, 2.0943951023931953, 3.1415926535897931, 4.1887902047863905, 5.2359877559829879]
  result &= approx_equal(one_slice_phi, reference_phi)
  return result

class HemisphereSamplerBase(SimpleSamplerTool):
  def __init__(self,characteristic_grid,max_cell):
    SimpleSamplerTool.__init__(self,characteristic_grid)
    self.construct_hemisphere_grid(self.incr)
    self.max_cell=max_cell
  def get_top_solutions(self,ai,input_directions,size,cutoff_divisor,grid):

    kval_cutoff = ai.getXyzSize()/cutoff_divisor;

    hemisphere_solutions = flex.Direction();
    hemisphere_solutions.reserve(size);
    #print "# of input directions", len(input_directions)
    for i in range(len(input_directions)):
      D = sampled_direction = ai.fft_result(input_directions[i])

      #hardcoded parameter in for silicon example.  Not sure at this point how
      #robust the algorithm is when this value is relaxed.

      if D.real < self.max_cell and sampled_direction.kval > kval_cutoff:
        if diagnostic:  directional_show(D, "%5d"%i)
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
      for i in range(len(unrefined_basis_vectors)):
        D = unrefined_basis_vectors[i];
        if diagnostic:  directional_show(D, "SE%5d"%i)

    ai.setSolutions(unrefined_basis_vectors)

def hemisphere_shortcut(ai,characteristic_sampling,max_cell):
    H = HemisphereSamplerBase(characteristic_grid = characteristic_sampling,max_cell=max_cell)
    H.hemisphere(ai,size=480,cutoff_divisor=1.5)
    #1.0/cutoff divisor is the height requirement (as a fraction of total Bragg spots)
    #required to consider a direction as a candidate basis vector.  4.0 good for
    #macromolecular work, 1.5 ok for small molecules; probably 1.1 would do.

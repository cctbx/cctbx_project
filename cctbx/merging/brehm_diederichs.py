from __future__ import absolute_import, division, print_function
from cctbx import sgtbx

from cctbx.array_family import flex
from libtbx.development.timers import Profiler
from cctbx import miller as miller_ext
import math
from cctbx.merging import update_wij_rij
from six.moves import range
"""
This module implements algorithm 2 from "Breaking the indexing ambiguity in
serial crystallography", Wolfgang Brehm & Kay Diederichs, Acta Cryst. D70 (2014)
"""
class algorithm2:

  def __init__(self, data, lattice_id, resort=False, verbose=True):
    self.verbose = verbose
    ######  INPUTS #######
    #       data = miller array: ASU miller index & intensity
    #       lattice_id = flex double: assignment of each miller index to a lattice number
    ######################

    if resort:
      order = flex.sort_permutation(lattice_id.data())
      sorted_lattice_id = flex.select(lattice_id.data(), order)
      sorted_data = data.data().select( order)
      sorted_indices = data.indices().select( order)
      self.lattice_id = sorted_lattice_id
      self.data = data.customized_copy(indices = sorted_indices, data = sorted_data)
    else:
      self.lattice_id = lattice_id.data() # type flex int
      self.data = data # type miller array with flex double data
    assert type(self.data.indices()) == type(flex.miller_index())
    assert type(self.data.data()) == type(flex.double())

    # construct a lookup for the separate lattices
    last_id = -1; self.lattices = flex.int()
    for n in range(len(self.lattice_id)):
      if self.lattice_id[n] != last_id:
        last_id = self.lattice_id[n]
        self.lattices.append(n)

  def generate_twin_operators(self):
    from mmtbx.scaling.twin_analyses import twin_laws
    TL = twin_laws(miller_array=self.data)
    if len(TL.operators) == 0:
      from libtbx.utils import Sorry
      raise Sorry("No twin laws are possible for this crystal lattice.")
    for twin_law in TL.operators:
      if self.verbose: print(twin_law.operator.r().as_hkl())
    return TL.operators

  def generate_reindex_sets(self):
    ops = self.generate_twin_operators()
    alternates = {}
    for op in ops:
      hkl = self.data.indices()
      cb_op = sgtbx.change_of_basis_op(op.operator.r().as_hkl())
      hklrev = cb_op.apply(hkl)
      alternates[op.operator.r().as_hkl()] = self.data.customized_copy(indices = hklrev).map_to_asu()
    return alternates

  def generate_comparison_selections(self,nproc=1):
    nproc = nproc-1
    if nproc<4:
      return [self.lattices>-1]
    result = []
    for n in range(nproc):
      array = flex.bool()
      for LL in range(len(self.lattices)):
        array.append((LL-n)%nproc < 3)
      result.append(array)
    return result
    # return value is a list of boolean arrays.  Each of length len(self.lattices),
    # corresponds to the number of lattices.

  def one_lattice_slice(self, indices, lattice_id):
    lower_index = self.lattices[lattice_id]
    if lattice_id < len(self.lattices)-1:
      upper_index = self.lattices[lattice_id+1]
      return indices[lower_index:upper_index]
    else:
      assert lattice_id == len(self.lattices)-1
      return indices[lower_index:]

  def run_core_algorithm(self, group, alternates, use_weights, asymmetric=1,
                         show_plot=True, save_plot=False, plot_name='xy.png'):
    # asymmetric 0=do nothing; 1=up/up; 2=down/up; 3=up/up + down/up
    #T = Profiler("coset")

    # construct rij matrix
    NN = group.count(True)

    if self.verbose: print("IN RUN CORE",group, alternates, use_weights, asymmetric,"Group of %d"%NN)
    index_selected = group.iselection()
    rij = []
    wij = []
    slices = {}
    for coset in ['h,k,l']+list(alternates.keys()):
      slices[coset] = {}
      twin_data = alternates.get(coset, self.data) # i.e., for 'h,k,l' the twin data is self.data
      for itr in range(NN):
        slices[coset][itr] = self.one_lattice_slice(indices = twin_data.indices(), lattice_id=index_selected[itr])

    for twin_law in alternates.keys():
      twin_data = alternates[twin_law]
      rij_ = flex.double(flex.grid(NN,NN),0.)
      wij_ = flex.double(flex.grid(NN,NN),0.)
      for i in range(NN):
        wij_[(i,i)] = 0.0
        indices_i     = slices["h,k,l"][i]
        indices_i_rev = slices[twin_law][i]
        i_start = self.lattices[index_selected[i]]

        # this would be a good compromise point for a detail call to C++ to speed things up XXX
        for j in range(i+1,NN):
          indices_j = slices["h,k,l"][j]
          j_start = self.lattices[index_selected[j]]

          if asymmetric%2 == 1: # up - up
            update_wij_rij(i,j,indices_i,indices_j,self.data.data(), self.data.data(),
                           i_start, j_start, wij_, rij_, 1., use_weights)
            #matches = miller_ext.match_indices(indices_i, indices_j)
            #intensities_i = flex.double()
            #intensities_j = flex.double()
            #for pair in matches.pairs():
               #print indices_i[pair[0]], indices_j[pair[1]],i_start,j_start,self.data.data()[i_start+pair[0]],self.data.data()[j_start+pair[1]]

            #   intensities_i.append( self.data.data()[i_start+pair[0]] )
            #   intensities_j.append( self.data.data()[j_start+pair[1]] )
            #corr = flex.linear_correlation(intensities_i, intensities_j)
            #if corr.is_well_defined():
            #  print i,j,corr.coefficient(),corr.n()
            #  if use_weights:
            #    wij_[(i,j)] = corr.n()
            #    wij_[(j,i)] = corr.n()
            #  rij_[(i,j)] = corr.coefficient()
            #  rij_[(j,i)] = corr.coefficient()

          if asymmetric >= 2: # down - up
            update_wij_rij(i,j,indices_i_rev,indices_j,twin_data.data(), self.data.data(),
                           i_start, j_start, wij_, rij_, -1., use_weights)
            #matches_rev = miller_ext.match_indices(indices_i_rev, indices_j)
            #intensities_i = flex.double()
            #intensities_j = flex.double()
            #for pair in matches_rev.pairs():
               #print indices_i_rev[pair[0]], indices_j[pair[1]],i_start,j_start,twin_data.data()[i_start+pair[0]],self.data.data()[j_start+pair[1]]

            #   intensities_i.append( twin_data.data()[i_start+pair[0]] )
            #   intensities_j.append( self.data.data()[j_start+pair[1]] )
            #corr = flex.linear_correlation(intensities_i, intensities_j)
            #if corr.is_well_defined():
            #  if use_weights:
            #    wij_[(i,j)] += corr.n()
            #    wij_[(j,i)] += corr.n()
            #  rij_[(i,j)] -= corr.coefficient()
            #  rij_[(j,i)] -= corr.coefficient()

      focus = wij_.focus()
      flat_wij = wij_.as_1d()
      selection = (flat_wij>1)
      print("w_ij is a %dx%d matrix with %d/%d >1 elements with average value %4.1f"%(
        focus[0],focus[1],selection.count(True),len(wij_), flex.mean(flat_wij.select(selection))))
      rij.append(rij_)
      wij.append(wij_)
    if self.verbose: print("CONSTRUCTED RIJ")
    xcoord = flex.random_double (NN)
    ycoord = flex.random_double (NN)
    M = minimize(xcoord,ycoord,rij[0],wij[0],self.verbose)
    coord_x = M.x[0:NN]
    coord_y = M.x[NN:2*NN]
    P = minimize_divide(coord_x, coord_y)
    selection = P.plus_minus()
    if show_plot or save_plot:
      import matplotlib
      if not show_plot:
        # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
        matplotlib.use('Agg') # use a non-interactive backend

      from matplotlib import pyplot as plt
      plt.plot(coord_x.select(selection),coord_y.select(selection),"r.", markersize=2.)
      plt.plot(coord_x.select(~selection),coord_y.select(~selection),"k.", markersize=3.)
      plt.axes().set_aspect("equal")
      if save_plot:
        plt.savefig(plot_name,
                    size_inches=(10,10),
                    dpi=300,
                    bbox_inches='tight')
      if show_plot:
        plt.show()
      plt.clf()
    grouped_lattice_ids = self.lattices.select(group)
    assert len(grouped_lattice_ids) == len(selection)

    from libtbx import group_args

    return group_args(reindexing_sets={
      "h,k,l":set(grouped_lattice_ids.select(selection)),
      # FIXME: order of .keys() changes depending on py2/3 , might break if len(key) > 1
      list(alternates.keys())[0]:set(grouped_lattice_ids.select(~selection))},
                      rij=rij,
                      wij=wij,
                      coord_x=coord_x,
                      coord_y=coord_y)

    # FIXME: verify len( keys) > 1 or else keys()[0] might break py2/3 compat
    return {"h,k,l":set(grouped_lattice_ids.select(selection)),
            list(alternates.keys())[0]:set(grouped_lattice_ids.select(~selection))}
    # Values in the return sets are indexes into the self.data data array (or the self.lattice_id array)
    # corresponding to the start of data for a particular lattice.  Need a different function call
    # to actually get the lattice ids.

  def __call__(self, input_queue, *args, **kwargs):
    # Scale frames sequentially within the current process.  The
    # return value is picked up by the callback.  See also
    # self.scale_all_serial()
    from Queue import Empty

    try :
      while True:
        try:
          igroup,group = input_queue.get_nowait()
        except Empty:
          return self

        result = self.run_core_algorithm(group, *args, **kwargs)

        input_queue.task_done()
        return igroup,result

    except Exception as e :
      print("Exception within __call__",e)
      return None

  def report(self,millerlookups):
     # TODO: verify what millerlookups is - is it a dict ?
     result = millerlookups.fromkeys(list(millerlookups.keys()))
     for key in result:
       result[key]=[]
       for item in millerlookups[key]:
         result[key].append(self.lattice_id[item])
     for key in result:
       result[key].sort()
     if self.verbose:
       print("REPORT")
       print(result)
     return result

from scitbx.lbfgs.tst_curvatures import lbfgs_with_curvatures_mix_in
class minimize(lbfgs_with_curvatures_mix_in):
  def __init__(self,xcoord,ycoord,rij_matrix,wij_matrix,verbose):
    self.verbose = verbose
    self.x = xcoord.concatenate(ycoord)
    self.NN = len(xcoord)
    self.rij_matrix = rij_matrix
    self.wij_matrix = wij_matrix
    if self.verbose:print(len(self.x))
    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      traditional_convergence_test_eps=1.0,
      use_curvatures=True)

  def curvatures(self):

    x_sq = self.x * self.x
    x_slice = x_sq[0:self.NN]
    y_slice = x_sq[self.NN:2*self.NN]
    return 2.*(self.wij_matrix.matrix_multiply(x_slice).concatenate(
             self.wij_matrix.matrix_multiply(y_slice)))

  def compute_functional_and_gradients(self):
    coord_x = self.x[0:self.NN]
    coord_y = self.x[self.NN:2*self.NN]

    inner = self.rij_matrix - coord_x.matrix_outer_product(coord_x) - coord_y.matrix_outer_product(coord_y)
    elements = self.wij_matrix*inner*inner
    f = 0.5 * flex.sum(elements)

    # quick gradients
    wrij_matrix = self.wij_matrix * self.rij_matrix
    term_1 = wrij_matrix.matrix_multiply(coord_x).concatenate(wrij_matrix.matrix_multiply(coord_y))
    temp_2 = self.wij_matrix * (coord_x.matrix_outer_product(coord_x))
    term_2x = (temp_2).matrix_multiply(coord_x)
    term_2y = (temp_2).matrix_multiply(coord_y)
    temp_3 = self.wij_matrix * (coord_y.matrix_outer_product(coord_y))
    term_3x = (temp_3).matrix_multiply(coord_x)
    term_3y = (temp_3).matrix_multiply(coord_y)
    term_2 = term_2x.concatenate(term_2y)
    term_3 = term_3x.concatenate(term_3y)
    grad = -2.* ( term_1 - term_2 - term_3 )

    if self.verbose: print("Functional",f)
    #from matplotlib import pyplot as plt
    #plt.plot(coord_x,coord_y,"r.")
    #plt.axes().set_aspect("equal")
    #plt.show()
    return f,grad

class minimize_divide(lbfgs_with_curvatures_mix_in):
  def __init__(self,xcoord,ycoord):
    self.xcoord = xcoord
    self.ycoord = ycoord
    self.x = flex.double([0])
    from scitbx.matrix import col
    self.center_of_mass = col((flex.mean(self.xcoord), flex.mean(self.ycoord)))

    grid = [ self.functional_only(t*math.pi/180.) for t in range(0,180,5) ]

    minvalue = min(grid)
    self.x = flex.double([5*grid.index(minvalue)*math.pi/180.])
    #print "grid_minimum",list(self.x)
    # XXX Looks like LBFGS is not performing minimization; investigate later
    # XXX can probably return here; 5-degree granularity enough; no need for LBFGS

    lbfgs_with_curvatures_mix_in.__init__(self,
      min_iterations=0,
      max_iterations=1000,
      use_curvatures=False)

  def functional_only(self,theta):
    from scitbx.matrix import col
    Bvec = col((math.cos(theta), math.sin(theta)))
    func = 0.
    for p in range(len(self.xcoord)):
      Avec = col((self.xcoord[p], self.ycoord[p]))-self.center_of_mass
      AdotB = Avec.dot(Bvec)
      Pvec = Avec - AdotB*Bvec # distance vector to dividing line
      func -= (Pvec.dot(Pvec))
    return func

  def compute_functional_and_gradients(self):
    from scitbx.matrix import col

    Bvec = col((math.cos(self.x[0]), math.sin(self.x[0])))
    func = 0.
    grad = 0.
    for p in range(len(self.xcoord)):
      Avec = col((self.xcoord[p], self.ycoord[p]))-self.center_of_mass
      AdotB = Avec.dot(Bvec)
      Pvec = Avec - AdotB*Bvec # distance vector to dividing line
      func -= (Pvec.dot(Pvec))
      cos_th = math.cos(self.x[0])
      sin_th = math.sin(self.x[0])
      grad += cos_th*(-Avec[0]*sin_th + Avec[1]*cos_th) - sin_th*(Avec[0]*cos_th + Avec[1]*sin_th)
    #print "Divide Functional",func,"with degrees=%7.2f"%(180.*self.x[0]/math.pi), self.x[0]
    return func,flex.double([2.* Avec.length() * grad])

  def plus_minus(self):
    from scitbx.matrix import col
    result = flex.bool()
    Bvec = col((math.cos(self.x[0]), math.sin(self.x[0]))) # dividing plane
    Nvec = Bvec.rotate_2d(90.,deg=True)
    for p in range(len(self.xcoord)):
      Avec = (col((self.xcoord[p], self.ycoord[p]))-self.center_of_mass).normalize()
      theta = math.acos(Nvec.dot(Avec))
      result.append(theta>=math.pi/2.)
    return result

def reassemble(patchwork,verbose=False):
  resorted = [patchwork[0]]
  # TODO what is resorted[0] ? Is it a dictionary, if not fix the code
  refkeys = resorted[0].keys()
  for ip in range(1,len(patchwork)):
    reference = resorted[-1]
    refkeys = list(reference.keys())
    test = patchwork[ip]
    newdict = {}
    for testkey in test.keys():
      scores = [len(test[testkey].intersection(reference[refkey])) for refkey in refkeys]
      max_score = max(scores)
      max_idx = scores.index(max_score)
      newdict[refkeys[max_idx]]=test[testkey]
    resorted.append(newdict)

  # Now all the patches have been lined up.
  voted = {}
  lookup = {refkeys[0]:0, refkeys[1]:1}
  for group in resorted:
    for key in group.keys():
      for item in group[key]:
        if voted.get(item, None) is None: voted[item]=[0]*len(refkeys)
        voted[item][lookup[key]] += 1
  if verbose:
    for key in voted.keys():
      print("vote",key, voted[key])
  # The votes are tallied, now report them
  result = {}
  for key in refkeys: result[key]=[]
  for key in voted.keys():
    max_score = max(voted[key])
    max_idx = voted[key].index(max_score)
    result[refkeys[max_idx]].append(key)
  return result

def mp_reassemble(mp_patchwork):
  clean_list = [0]*len(mp_patchwork)
  for ituple in mp_patchwork:
    clean_list[ituple[0]]=ituple[1]
  return reassemble(clean_list)

def run(L, asymmetric=3, nproc=1,verbose=True, show_plot=True, save_plot=False):
  algo2 = algorithm2(data = L[0],lattice_id = L[1], resort=True, verbose = verbose)
  alternates = algo2.generate_reindex_sets()
  c_selections = algo2.generate_comparison_selections(nproc=nproc)
  result_sets = []
  for i, group in enumerate(c_selections):
    result_sets.append(
      algo2.run_core_algorithm(
        group, alternates, use_weights=True, asymmetric=asymmetric,
        show_plot=show_plot, save_plot=save_plot,
        plot_name='xy_%i.png' %i).reindexing_sets
    )
  result = reassemble(result_sets)
  return algo2.report(result)

def run_multiprocess(L, asymmetric=3, nproc=20, verbose=False, show_plot=True,
                     save_plot=False):
  try :
      import multiprocessing
  except ImportError as e :
      print("multiprocessing module not available (requires Python >= 2.6)"); exit()
  algo2 = algorithm2(data = L[0],lattice_id = L[1], resort=True, verbose = verbose)
  alternates = algo2.generate_reindex_sets()
  #import libtbx.introspection
  #nproc = libtbx.introspection.number_of_processors(return_value_if_unknown=1)
  c_selections = algo2.generate_comparison_selections(nproc=nproc)
  result_sets = []
  def _mpcallback(result_set):
    import libtbx
    if verbose: print("IN MPCALLBACK with",result_set)
    if type(result_set)==type(()) and type(result_set[1])==libtbx.group_args:
      #result_sets.append(result_set.reindexing_sets)
      result_sets.append((result_set[0], result_set[1].reindexing_sets))
  input_queue = multiprocessing.Manager().JoinableQueue()
  for ic,item in enumerate(c_selections):
    input_queue.put((ic,item))

  pool = multiprocessing.Pool(processes=nproc-1)
  # Each process accumulates its own statistics in serial, and the
  # grand total is eventually collected by the main process'
  #  _mpcallback
  for i in range(nproc-1):
    pool.apply_async(
      func=algo2,
      args=[input_queue, alternates, True],
      kwds={'show_plot': show_plot,
            'save_plot': save_plot,
            'plot_name': 'xy_%i.png' %i,
            'asymmetric': asymmetric},
      callback=_mpcallback)
  pool.close()
  pool.join()

  # Block until the input queue has been emptied.
  input_queue.join()

  result = mp_reassemble(result_sets)
  return algo2.report(result)

def test_reassembly():
  testdata = [{'h,k,l': set([0, 231041, 871318, 324871, 229780, 371093, 830358, 650903, 82465, 694052, 40613, 695591, 831658, 512945, 436, 38837, 180028, 178613, 373056, 467266, 419397, 281288, 651855, 131036, 511458, 607588, 872805, 323944, 464747, 783720, 421748, 232821]), 'h,-h-k,-l': set([81027, 282758, 607499, 131985, 279826, 371459, 420632, 652581, 782887, 510123, 466097, 1075, 83895, 607928, 129721, 830267, 557634, 738378, 559949, 693072, 736859, 784605, 739942, 871528, 176745, 558701, 322682, 39806])}, {'h,k,l': set([231041, 514309, 181511, 832779, 1806, 830358, 85024, 82465, 653730, 694052, 40613, 695591, 284072, 831658, 324871, 512945, 436, 178613, 608822, 180028, 373056, 467266, 281288, 651855, 607588, 131036, 783720, 511458, 374371, 41316, 872805, 323944, 325864, 421748, 232821, 874097]), 'h,-h-k,-l': set([371459, 282758, 131985, 420632, 741284, 652581, 786344, 466097, 1075, 83895, 607928, 738378, 559949, 132944, 234327, 467928, 784605, 422750, 561120, 739942, 696856, 871528, 558701, 39806])}, {'h,k,l': set([282758, 87383, 131985, 696856, 562073, 609818, 741284, 652581, 833841, 786344, 327338, 182705, 1075, 83895, 607928, 742467, 559949, 787150, 132944, 234327, 467928, 784605, 422750, 561120, 739942, 655593, 515438, 42104]), 'h,-h-k,-l': set([235777, 514309, 181511, 832779, 1806, 134171, 85024, 653730, 40613, 284454, 695591, 284072, 831658, 324871, 512945, 608822, 180028, 373056, 467266, 424260, 698181, 376017, 3033, 374371, 41316, 872805, 325864, 468592, 874097, 421748, 232821, 875516])}, {'h,k,l': set([285959, 234327, 696856, 562073, 609818, 741284, 833841, 786344, 656554, 182705, 469941, 425534, 742467, 787150, 132944, 42104, 88662, 87383, 467928, 516442, 563292, 422750, 561120, 655593, 788717, 515438, 611192, 327338]), 'h,-h-k,-l': set([235777, 514309, 181511, 42379, 184076, 743693, 1806, 328984, 134171, 85024, 653730, 284454, 284072, 4010, 237106, 608822, 699703, 832779, 424260, 698181, 377414, 876109, 376017, 835029, 3033, 374371, 41316, 325864, 468592, 874097, 875516, 135466])}, {'h,k,l': set([378625, 657922, 285959, 185363, 562073, 563993, 609818, 836385, 833841, 517852, 327338, 182705, 43829, 286777, 425534, 469941, 742467, 787150, 611192, 426708, 88662, 87383, 789721, 516442, 563292, 89570, 655593, 470891, 788717, 515438, 4592, 42104, 611834, 656554]), 'h,-h-k,-l': set([235777, 876937, 42379, 184076, 743693, 3033, 328984, 134171, 284454, 4010, 237106, 699703, 700992, 424260, 698181, 377414, 876109, 376017, 835029, 745177, 238049, 136036, 468592, 329588, 875516, 135466])}, {'h,k,l': set([378625, 657922, 517852, 285959, 185363, 563993, 790942, 836385, 136873, 656554, 186670, 330676, 43829, 286777, 425534, 469941, 837963, 746693, 378699, 88662, 789721, 516442, 563292, 519263, 89570, 470891, 788717, 238958, 4592, 5876, 471925, 611192, 611834]), 'h,-h-k,-l': set([876937, 42379, 184076, 743693, 328984, 287769, 612636, 658211, 565148, 4010, 237106, 699703, 44729, 700992, 426821, 377414, 876109, 701522, 426708, 835029, 90584, 745177, 238049, 136036, 878191, 329588, 135466])}, {'h,k,l': set([45442, 876937, 188170, 289040, 838292, 472983, 287769, 748059, 612636, 91551, 331555, 565148, 879280, 44729, 700992, 792388, 426821, 658637, 701522, 658211, 90584, 702553, 238049, 136036, 613479, 745177, 520939, 878191, 329588]), 'h,-h-k,-l': set([378625, 657922, 379909, 138253, 185363, 7534, 563993, 790942, 836385, 136873, 186670, 240692, 43829, 286777, 566207, 837963, 746693, 378699, 426708, 789721, 517852, 519263, 89570, 470891, 238958, 4592, 5876, 471925, 611834, 330676, 427773])}, {'h,k,l': set([839169, 332677, 138253, 7534, 381336, 790942, 379909, 136873, 186670, 240692, 330676, 139423, 880574, 566207, 837963, 746693, 378699, 428366, 519263, 660197, 792807, 238958, 5876, 471925, 427773]), 'h,-h-k,-l': set([474753, 45442, 46596, 242569, 188170, 289040, 189528, 838292, 472983, 8856, 287769, 748059, 565148, 91551, 879280, 331555, 612636, 615088, 521778, 792388, 44729, 567106, 749016, 426821, 92748, 658637, 701522, 658211, 90584, 702553, 290526, 613479, 520939, 703045, 878191])}, {'h,k,l': set([474753, 45442, 881131, 46596, 242569, 188170, 289040, 749016, 838292, 750681, 472983, 8856, 748059, 91551, 879280, 331555, 661033, 615088, 521778, 568758, 704373, 567106, 792388, 140613, 92748, 658637, 793681, 703045, 189528, 702553, 290526, 47203, 613479, 520939, 292077, 476014, 9973]), 'h,-h-k,-l': set([839169, 332677, 138253, 244496, 191379, 381336, 379909, 616230, 240692, 382264, 139423, 880574, 566207, 428366, 429325, 334392, 522587, 660197, 792807, 93803, 7534, 840562, 427773])}, {'h,k,l': set([474753, 46596, 242569, 749016, 335637, 8856, 794909, 703045, 383521, 246309, 661033, 615088, 521778, 568758, 704373, 567106, 140613, 92748, 141901, 793681, 193238, 189528, 750681, 290526, 292962, 47203, 881131, 292077, 476014, 94963, 9973, 751785, 841343]), 'h,-h-k,-l': set([839169, 93803, 332677, 429325, 244496, 191379, 381336, 139423, 430756, 48293, 616230, 662062, 382264, 570940, 880574, 882111, 616651, 428366, 334392, 477112, 522587, 705637, 11618, 660197, 792807, 524011, 840562])}, {'h,k,l': set([796681, 617998, 568758, 335637, 96151, 794909, 49825, 431396, 246309, 616230, 661033, 336307, 142388, 294454, 704373, 140613, 383521, 141901, 793681, 193238, 750681, 525530, 292962, 47203, 881131, 292077, 476014, 94963, 9973, 751785, 841343]), 'h,-h-k,-l': set([93803, 429325, 384656, 191379, 706465, 194722, 430756, 48293, 13222, 662062, 753073, 382264, 883515, 570940, 882111, 842565, 616651, 663202, 477112, 572630, 522587, 705637, 478304, 244496, 11618, 247397, 334392, 524011, 840562])}, {'h,k,l': set([572630, 384656, 884885, 663877, 195361, 194722, 430756, 48293, 13222, 662062, 753073, 479158, 477112, 883515, 570940, 882111, 574021, 706465, 616651, 663202, 337489, 526678, 798169, 705637, 478304, 11618, 247397, 618087, 524011, 843503, 842565, 14591]), 'h,-h-k,-l': set([754305, 144135, 796681, 617998, 433170, 335637, 96151, 295704, 794909, 97438, 383521, 431396, 246309, 248487, 751785, 336307, 142388, 294454, 49825, 141901, 193238, 525530, 51166, 292962, 385636, 707695, 94963, 841343])}, {'h,k,l': set([480640, 572630, 435082, 384656, 296979, 884885, 663877, 195361, 194722, 13222, 144555, 753073, 479158, 52282, 883515, 574782, 574021, 706465, 845253, 663202, 337489, 799701, 526678, 798169, 478304, 619364, 247397, 618087, 843503, 842565, 14591]), 'h,-h-k,-l': set([754305, 886147, 144135, 796681, 15426, 617998, 433170, 96151, 295704, 707695, 97438, 49825, 98979, 431396, 248487, 336307, 142388, 294454, 709184, 338882, 249922, 665418, 756301, 527758, 386980, 525530, 51166, 385636, 196719])}, {'h,k,l': set([754305, 886147, 145924, 144135, 15426, 527758, 433170, 576278, 295704, 707695, 100124, 97438, 98979, 386980, 248487, 436231, 757174, 801087, 709184, 619969, 338882, 249922, 665418, 756301, 710866, 51166, 16735, 385636, 196719, 340090, 481661, 846334]), 'h,-h-k,-l': set([480640, 666838, 52873, 435082, 296979, 884885, 574021, 195361, 144555, 251181, 479158, 52282, 574782, 663877, 845253, 528717, 337489, 799701, 526678, 388184, 798169, 197471, 619364, 618087, 297582, 843503, 887167, 14591])}, {'h,k,l': set([886147, 145924, 436231, 15426, 527758, 621456, 576278, 100124, 388639, 711584, 98979, 386980, 757174, 198968, 801087, 709184, 619969, 338882, 53571, 249922, 665418, 252491, 756301, 710866, 16735, 18402, 196719, 340090, 481661, 846334]), 'h,-h-k,-l': set([480640, 52873, 101386, 296979, 147735, 668696, 437539, 144555, 251181, 802485, 52282, 435082, 574782, 577355, 845253, 298315, 528717, 888274, 799701, 666838, 388184, 482521, 197471, 619364, 297582, 847347, 529611, 758649, 341372, 887167])}, {'h,k,l': set([102401, 712706, 145924, 669782, 436231, 621456, 530833, 576278, 100124, 388639, 711584, 804018, 757174, 198968, 484665, 801087, 619969, 53571, 252491, 298962, 54486, 578265, 16735, 18402, 710866, 759921, 439032, 340090, 253051, 481661, 846334]), 'h,-h-k,-l': set([52873, 101386, 147735, 668696, 437539, 251181, 848433, 200629, 889275, 802485, 577355, 18629, 298315, 528717, 390350, 888274, 666838, 148695, 388184, 482521, 342748, 197471, 622819, 297582, 847347, 529611, 758649, 341372, 887167])}, {'h,k,l': set([101386, 147735, 668696, 671262, 149664, 437539, 532263, 253737, 343850, 55601, 200629, 202424, 485433, 19771, 849854, 802485, 577355, 18629, 298315, 390350, 888274, 848433, 148695, 482521, 623963, 342748, 622819, 847347, 529611, 758649, 341372, 889275, 299391]), 'h,-h-k,-l': set([102401, 712706, 54486, 621456, 530833, 760982, 388639, 711584, 103332, 804018, 198968, 484665, 579521, 53571, 252491, 805199, 391632, 298962, 439764, 669782, 890456, 578265, 18402, 713702, 759921, 439032, 253051])}, {'h,k,l': set([102401, 712706, 54486, 762632, 851081, 625419, 530833, 760982, 254751, 103332, 486695, 804018, 484665, 579521, 104258, 56647, 805199, 391632, 298962, 439764, 669782, 393047, 890456, 578265, 20703, 713702, 580328, 759921, 439032, 253051, 203646]), 'h,-h-k,-l': set([300545, 715140, 671262, 149664, 672037, 532263, 253737, 343850, 55601, 200629, 202424, 485433, 19771, 849854, 18629, 806089, 390350, 533458, 848433, 148695, 623963, 342748, 622819, 440945, 344569, 151164, 889275, 299391])}, {'h,k,l': set([300545, 715140, 346258, 487831, 671262, 149664, 716065, 762916, 672037, 532263, 253737, 343850, 55601, 442165, 202424, 485433, 19771, 849854, 806089, 105272, 533458, 623963, 204384, 440945, 344569, 151164, 299391]), 'h,-h-k,-l': set([851842, 762632, 851081, 255882, 625419, 626574, 760982, 152859, 254751, 103332, 486695, 581676, 807215, 21471, 394378, 579521, 104258, 56647, 534601, 805199, 391632, 439764, 57686, 393047, 890456, 20703, 713702, 301159, 580328, 673262, 203646])}, {'h,k,l': set([851842, 762632, 851081, 255882, 625419, 626574, 535831, 152859, 254751, 486695, 347818, 763947, 581676, 807215, 257209, 21471, 394378, 104258, 56647, 534601, 57686, 393047, 20703, 301159, 580328, 673262, 153972, 443382, 808825, 203646]), 'h,-h-k,-l': set([300545, 715140, 626950, 346258, 852116, 487831, 106909, 716065, 762916, 672037, 442165, 806089, 105272, 22205, 206277, 489030, 583241, 59082, 674125, 533458, 717049, 302554, 204384, 395883, 440945, 344569, 151164])}, {'h,k,l': set([626950, 346258, 852116, 487831, 628248, 584345, 106909, 716065, 762916, 442165, 105272, 22205, 155202, 206277, 489030, 583241, 59082, 674125, 397135, 302554, 208091, 303452, 258911, 204384, 395883, 717049, 59900]), 'h,-h-k,-l': set([851842, 853675, 255882, 349463, 626574, 810350, 535831, 152859, 675100, 489631, 347818, 763947, 581676, 807215, 257209, 22586, 394378, 766019, 718536, 534601, 57686, 108126, 21471, 537568, 445026, 301159, 673262, 153972, 443382, 808825])}, {'h,k,l': set([626950, 304527, 852116, 628248, 584345, 106909, 766129, 446780, 22205, 155202, 206277, 489030, 583241, 59082, 674125, 397135, 585942, 350681, 302554, 208091, 303452, 258911, 109540, 156391, 675690, 395883, 260336, 629750, 717049, 59900]), 'h,-h-k,-l': set([853675, 349463, 535831, 209374, 675100, 489631, 347818, 763947, 718536, 719326, 257209, 22586, 766019, 61128, 397899, 538960, 490969, 108126, 810463, 537568, 445026, 854371, 810350, 23023, 153972, 443382, 808825])}, {'h,k,l': set([540039, 766866, 349463, 677532, 489631, 587298, 675100, 853675, 718536, 630962, 209374, 398905, 22586, 766019, 61128, 397899, 538960, 719326, 490969, 352221, 108126, 810463, 537568, 445026, 854371, 210919, 157674, 810350, 23023, 811762]), 'h,-h-k,-l': set([305286, 261512, 447170, 304527, 628248, 584345, 23850, 766129, 446780, 155202, 855750, 62414, 397135, 585942, 350681, 208091, 303452, 258911, 109540, 156391, 492265, 675690, 720999, 111599, 260336, 629750, 59900])}, {'h,k,l': set([540039, 24456, 212233, 766866, 677532, 587298, 630962, 719326, 262665, 398905, 61128, 397899, 588578, 538960, 447828, 493526, 490969, 352221, 209374, 810463, 399969, 854371, 767588, 305894, 210919, 157674, 23023, 811762]), 'h,-h-k,-l': set([305286, 261512, 112267, 540610, 304527, 158997, 632218, 813085, 678944, 23850, 766129, 856633, 446780, 447170, 353475, 855750, 63048, 721994, 62414, 585942, 350681, 109540, 156391, 492265, 675690, 720999, 111599, 260336, 629750])}, {'h,k,l': set([540039, 24456, 212233, 766866, 633493, 677532, 495007, 587298, 768811, 630962, 262665, 398905, 814279, 588578, 542287, 447828, 493526, 25690, 352221, 399969, 767588, 305894, 210919, 157674, 857968, 811762, 160501, 307450, 64253]), 'h,-h-k,-l': set([305286, 261512, 112267, 540610, 158997, 632218, 590363, 448668, 813085, 678944, 679558, 401702, 354727, 23850, 213935, 856633, 113345, 447170, 353475, 855750, 63048, 721994, 62414, 264166, 720999, 492265, 111599, 723955])}, {'h,k,l': set([24456, 262665, 633493, 634520, 114588, 495007, 588578, 725467, 768811, 212233, 859333, 161734, 814279, 449871, 215248, 447828, 493526, 25690, 542287, 399969, 767588, 305894, 264294, 25961, 857968, 160501, 307450, 64253]), 'h,-h-k,-l': set([679558, 65415, 112267, 681229, 158997, 632218, 590363, 448668, 813085, 678944, 401702, 354727, 308700, 213935, 403255, 856633, 591932, 113345, 540610, 353475, 63048, 721994, 770256, 496091, 356188, 264166, 543335, 723955, 815350])}, {'h,k,l': set([679558, 65415, 682378, 681229, 265746, 162836, 590363, 448668, 401702, 354727, 308700, 213935, 403255, 591932, 113345, 66506, 404557, 27214, 770256, 593745, 450642, 496091, 356188, 264166, 543335, 723955, 635765, 815350, 116220]), 'h,-h-k,-l': set([633493, 634520, 114588, 727325, 495007, 816475, 768811, 356674, 859333, 161734, 814279, 860492, 496846, 449871, 544208, 25690, 725467, 309468, 215248, 264294, 216680, 25961, 771695, 857968, 542287, 160501, 307450, 64253])}, {'h,k,l': set([28172, 406157, 683537, 634520, 114588, 727325, 67998, 816475, 309468, 310573, 356674, 859333, 161734, 117323, 860492, 496846, 449871, 215248, 594898, 725467, 637404, 357469, 544208, 264294, 216680, 25961, 771695, 164336, 218100]), 'h,-h-k,-l': set([65415, 682378, 267659, 681229, 265746, 162836, 861718, 544669, 496091, 308700, 451756, 403255, 591932, 66506, 404557, 27214, 817103, 770256, 593745, 450642, 498267, 356188, 772964, 543335, 635765, 815350, 727545, 116220])}, {'h,k,l': set([28172, 406157, 683537, 773400, 727325, 67998, 69792, 309468, 310573, 118067, 499767, 683708, 269250, 818757, 860492, 117323, 356674, 496846, 544208, 594898, 165975, 816475, 637404, 357469, 29413, 216680, 311406, 771695, 164336, 218100]), 'h,-h-k,-l': set([358921, 682378, 267659, 596241, 265746, 162836, 861718, 544669, 451756, 407606, 66506, 218443, 638668, 404557, 27214, 817103, 593745, 450642, 545486, 498267, 728289, 772964, 453482, 635765, 727545, 116220, 862335])}, {'h,k,l': set([597760, 28172, 406157, 683537, 820114, 684947, 500886, 773400, 67998, 69792, 164336, 310573, 118067, 774070, 499767, 683708, 269250, 818757, 219463, 639756, 117323, 546383, 594898, 70998, 165975, 863322, 637404, 357469, 29413, 408934, 311406, 167152, 218100, 31102, 729215]), 'h,-h-k,-l': set([358921, 267659, 596241, 119186, 861718, 544669, 454942, 451756, 359858, 407606, 312135, 218443, 638668, 545486, 817103, 498267, 728289, 772964, 453482, 727545, 269694, 862335])}, {'h,k,l': set([597760, 408934, 639756, 820114, 684947, 456213, 500886, 773400, 410777, 821658, 69792, 774832, 118067, 168245, 774070, 502071, 683708, 269250, 818757, 219463, 499767, 546383, 70998, 165975, 863322, 641372, 361053, 29413, 270182, 311406, 167152, 31102, 729215]), 'h,-h-k,-l': set([220804, 358921, 596241, 119186, 686227, 312858, 454942, 730539, 32178, 359858, 547252, 407606, 312135, 218443, 638668, 545486, 72413, 728289, 121189, 453482, 598902, 864764, 269694, 862335])}, {'h,k,l': set([313984, 597760, 865794, 639756, 820114, 684947, 456213, 500886, 410777, 821658, 362399, 169636, 774832, 168245, 774070, 502071, 219463, 548813, 546383, 687480, 70998, 863322, 641372, 361053, 270182, 408934, 167152, 122360, 31102, 729215]), 'h,-h-k,-l': set([457729, 730539, 220804, 822925, 412046, 119186, 686227, 312858, 503195, 454942, 775844, 271529, 73159, 32178, 731357, 359858, 547252, 312135, 33108, 600535, 72413, 121189, 222439, 642924, 598902, 864764, 269694])}, {'h,k,l': set([313984, 865794, 550030, 456213, 410777, 821658, 362399, 170915, 169636, 688550, 641372, 34602, 774832, 314929, 168245, 502071, 601537, 731976, 644169, 548813, 687480, 504796, 361053, 223713, 270182, 122360, 458873, 74236]), 'h,-h-k,-l': set([457729, 220804, 123788, 822925, 412046, 686227, 413718, 312858, 503195, 775844, 271529, 730539, 824110, 731357, 32178, 547252, 273222, 73159, 33108, 363990, 600535, 72413, 121189, 222439, 642924, 867183, 598902, 777210, 864764])}, {'h,k,l': set([457729, 35973, 123788, 822925, 412046, 824110, 413718, 503195, 688669, 778399, 775844, 271529, 274247, 125102, 316209, 645687, 414778, 273222, 73159, 171340, 505551, 33108, 363990, 600535, 731357, 222439, 867817, 642924, 867183, 602353, 777210]), 'h,-h-k,-l': set([313984, 865794, 225027, 550030, 362399, 170915, 169636, 688550, 34602, 314929, 460212, 732981, 825396, 601537, 731976, 644169, 365386, 75211, 548813, 551760, 122360, 504796, 223713, 687480, 458873, 74236])}, {'h,k,l': set([869249, 225027, 550030, 275471, 415521, 170915, 688550, 366887, 34602, 314929, 460212, 732981, 825396, 601537, 731976, 644169, 365386, 75211, 551760, 506745, 504796, 223713, 317416, 461431, 458873, 74236, 780031]), 'h,-h-k,-l': set([603778, 824110, 76932, 35973, 123788, 646933, 413718, 552964, 36634, 688669, 778399, 733466, 273222, 826663, 689964, 125102, 316209, 126259, 645687, 414778, 173126, 274247, 171340, 505551, 363990, 867817, 867183, 226800, 602353, 777210])}, {'h,k,l': set([869249, 174210, 225027, 735113, 368268, 275471, 462992, 415521, 366887, 276776, 460212, 732981, 79799, 825396, 507977, 365386, 75211, 551760, 127318, 604760, 317416, 319724, 461431, 506745, 780031]), 'h,-h-k,-l': set([603778, 76932, 35973, 416909, 646933, 781080, 552964, 36634, 227866, 778399, 733466, 826663, 689964, 125102, 688669, 691376, 316209, 126259, 647990, 645687, 414778, 173126, 274247, 171340, 505551, 554322, 827962, 36958, 867817, 226800, 602353, 869495])}, {'h,k,l': set([603778, 76932, 416909, 870161, 646933, 781080, 552964, 733466, 36634, 369438, 826663, 689964, 691376, 736049, 829362, 126259, 647990, 827962, 173126, 227866, 554322, 36958, 278112, 128482, 226800, 869495, 228217]), 'h,-h-k,-l': set([869249, 174210, 321159, 735113, 368268, 275471, 462992, 649624, 415521, 781817, 366887, 276776, 417322, 79799, 555834, 175676, 606145, 507977, 37460, 127318, 508759, 604760, 464143, 80605, 317416, 319724, 691825, 461431, 506745, 780031])}, {'h,k,l': set([81027, 607499, 416909, 870161, 279826, 781080, 227866, 369438, 782887, 510123, 691376, 829362, 647990, 129721, 827962, 830267, 557634, 693072, 554322, 736859, 36958, 278112, 128482, 176745, 869495, 228217, 322682]), 'h,-h-k,-l': set([0, 174210, 321159, 735113, 368268, 464143, 462992, 229780, 371093, 871318, 650903, 649624, 276776, 417322, 736049, 38837, 79799, 555834, 175676, 606145, 419397, 507977, 37460, 127318, 508759, 604760, 80605, 464747, 319724, 691825, 781817])}, {'h,k,l': set([0, 231041, 871318, 321159, 464143, 229780, 371093, 830358, 650903, 649624, 82465, 694052, 417322, 736049, 436, 38837, 555834, 175676, 178613, 606145, 419397, 281288, 651855, 37460, 508759, 131036, 80605, 783720, 511458, 607588, 323944, 464747, 691825, 781817]), 'h,-h-k,-l': set([81027, 607499, 870161, 279826, 371459, 420632, 369438, 782887, 510123, 466097, 829362, 129721, 830267, 557634, 738378, 693072, 736859, 278112, 128482, 871528, 176745, 558701, 228217, 322682, 39806])}]
  result = reassemble(testdata,verbose=False)
  assert len(result["h,k,l"])==392
  assert len(result["h,-h-k,-l"])==368
  print("OK")

if __name__=="__main__":
  test_reassembly()

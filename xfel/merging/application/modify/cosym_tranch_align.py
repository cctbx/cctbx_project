from __future__ import absolute_import, division, print_function
from itertools import permutations
import copy
from cctbx.array_family import flex
import numpy as np

"""Details of tranch alignment."""

def get_proposal_score(reports_as_lists):
    n_proposal_score = 0
    NN = len(reports_as_lists) # = n_tranches x n_permutations
    rij = flex.double(flex.grid(NN,NN),0.)
    wij = flex.double(flex.grid(NN,NN),1.)

    for ix in range(len(reports_as_lists)):
      base_report = reports_as_lists[ix]
      # compute the unions
      base_all = set(base_report[0])
      for icoset in range(1, len(base_report)):
        base_all = base_all.union(set(base_report[icoset]))

      for iy in range(len(reports_as_lists)):
        test_report = reports_as_lists[iy]
        matches = 0
        # compute the unions
        test_all = set(test_report[0])
        for icoset in range(1, len(test_report)):
          test_all = test_all.union(set(test_report[icoset]))
        # total overlap between base and test irrespective of cosets;
        total_overlay = len(base_all.intersection(test_all))

        for icoset in range(len(test_report)):
          matches += len(set(base_report[icoset]).intersection(set(test_report[icoset])))
        print ('%3d/%3d'%(matches, total_overlay),end=' ')
        rij [(ix,iy)] = matches/total_overlay if total_overlay>0 else 0.
        wij [(ix,iy)] = total_overlay if total_overlay>0 else 1.
        n_proposal_score += matches
      print()
    return rij, wij

def alignment_by_embedding(reports,plot=False):
  #from IPython import embed; embed()
  """reports is a list of tranch results, one list item per composite tranch
     Each item is a list, over cosets, e.g. two elements for a merohedral twinning op
     Each element is itself a list of uuid's assigned to that coset.
  """
  n_tranches = len(reports)
  n_operators = len(reports[0])
  reports = copy.deepcopy(reports)
  # will now amend the reports-list so that it has a new section for each permutation of cosets.

  for itranch in range(len(reports)):
    cache_permutations = list(permutations(reports[itranch]))
    assert len(reports[itranch]) == n_operators
    for jperm in range(1,len(cache_permutations)):
      reports.append(list(cache_permutations[jperm]))
    # Code is only valid for single twin operator, see alignment_by_seqential_trial() below

  rij, wij = get_proposal_score(reports)
  mt = flex.mersenne_twister(seed=0)
  #from IPython import embed; embed()
  NN = len(reports)
  xcoord = mt.random_double (size=NN)
  ycoord = mt.random_double (size=NN)
  if plot:
    from matplotlib import pyplot as plt
    plt.title("Initial spread before refinement")
    plt.plot(xcoord, ycoord, "g.")
    plt.show()

  from cctbx.merging.brehm_diederichs import minimize as mz, minimize_divide
  M = mz(xcoord,ycoord,rij,wij,verbose=True)
  coord_x = M.x[0:NN]
  coord_y = M.x[NN:2*NN]
  if plot:
    plt.title("Raw refinement")
    plt.plot(coord_x,coord_y,"r.", markersize=1.)
    plt.show()
  P = minimize_divide(coord_x, coord_y)
  selection = P.plus_minus()
  if plot:
    plt.title("After refinement")
    plt.plot(coord_x.select(selection),coord_y.select(selection),"r.", markersize=2.)
    plt.plot(coord_x.select(~selection),coord_y.select(~selection),"k.", markersize=3.)
    plt.show()

  print (list(selection))
  reformed_reports = [[] for i in range(n_tranches)] # output should have as many reports as tranches
  n_permutations = NN // n_tranches
  for iflag, select_flag in enumerate(selection):
    if select_flag:
      itranch = iflag % n_tranches
      #print( itranch, iflag, n_permutations)
      reformed_reports[itranch] = reports[iflag]

  assert [] not in reformed_reports # all tranches must have coset assignments
  return reformed_reports

def counts_for_sequential_trial(reports, active_subset):
  """Once the tranches are aligned
         tranche0 tranche1 tranche2
   coset0    A       C        E
   coset1    B       D        F
   The total number of datasets is union (A,....,F)
   This should also be equal to union (A,C,E) + union (B,D,F)
   Not true prior to alignment.
  """
  n_tranches = len(reports)
  n_cosets = len(reports[0])
  coset_totals = [set() for c in range(n_cosets)]
  total = set()
  for itranche in range(n_tranches):
    if not active_subset[itranche]: continue
    for icoset in range(n_cosets):
      cut = set(reports[itranche][icoset])
      total = total.union(cut)
      coset_totals[icoset] = coset_totals[icoset].union(cut)
  test = 0
  for icoset in range(n_cosets):
    test += len(coset_totals[icoset])
  print("Aligned total is",test,"grand total is",len(total), "delta",test-len(total))
  return test-len(total)

def grand_total_for_sequential_trial(reports, active_subset):
  n_tranches = len(reports)
  n_cosets = len(reports[0])
  total = set()
  for itranche in range(n_tranches):
    if not active_subset[itranche]: continue
    for icoset in range(n_cosets):
      cut = set(reports[itranche][icoset])
      total = total.union(cut)
  return len(total)

def alignment_by_sequential_trial(reports, plot=False):
  reports = copy.deepcopy(reports)
  active_subset_base = [False for i in reports]
  active_subset_base[0]=True # always start with tranche 0

  while active_subset_base.count(True)<len(active_subset_base):
    print("BEGIN LOOP with",active_subset_base.count(True),"aligned tranches")

    if False: # shortcut, fails
      selected_index = active_subset_base.index(False)
    else:
      #try to find the max overlapping tranche candidate (smallest grand total)
      base_total = grand_total_for_sequential_trial(reports, active_subset_base)
      candidates = np.array([np.iinfo(np.int64).max]*len(reports)) # all have max values
      for itranche in range(len(reports)):
        if not active_subset_base[itranche]:
          active_subset = copy.deepcopy(active_subset_base)
          active_subset[itranche]=True
          candidates[itranche]=grand_total_for_sequential_trial(reports, active_subset)
          print(itranche,base_total,candidates[itranche])
      selected_index = np.argmin(candidates)
    print("The selected tranche is",selected_index)

    # now having identified the candidate
    active_subset_base[selected_index]=True
    cache_permutations = list(permutations(reports[selected_index]))
    n_perm = len(cache_permutations)
    candidates = np.array([np.iinfo(np.int64).max]*n_perm) # all have max values
    for iperm in range(n_perm):
      reformed = copy.deepcopy(reports)
      reformed[selected_index] = list(cache_permutations[iperm])
      candidates[iperm] = counts_for_sequential_trial(reformed, active_subset_base)
    selected_perm = np.argmin(candidates)
    print ("The selected perm is",selected_perm)
    reports[selected_index] = list(cache_permutations[selected_perm])
    print("END LOOP with",active_subset_base.count(True),"aligned tranches")
  return reports

from __future__ import absolute_import, division, print_function
from itertools import permutations
import copy
from cctbx.array_family import flex

"""Details of tranch alignment."""

def get_proposal_score(reports_as_lists):
    n_proposal_score = 0
    NN = len(reports_as_lists) # = n_tranches x n_operators
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
  reports = copy.deepcopy(reports)
  # will now amend the reports-list so that it has a new section for each permutation of cosets.

  for itranch in range(len(reports)):
    cache_permutations = list(permutations(reports[itranch]))
    reports.append(cache_permutations[1])
    # will have to rewrite this code if there is more that one symmetry operator XXX FIXME

  rij, wij = get_proposal_score(reports)
  mt = flex.mersenne_twister(seed=0)
  #from IPython import embed; embed()
  NN = len(reports)
  xcoord = mt.random_double (size=NN)
  ycoord = mt.random_double (size=NN)
  if plot:
    from matplotlib import pyplot as plt
    plt.plot(xcoord, ycoord, "g.")
    plt.show()

  from cctbx.merging.brehm_diederichs import minimize as mz, minimize_divide
  M = mz(xcoord,ycoord,rij,wij,verbose=True)
  coord_x = M.x[0:NN]
  coord_y = M.x[NN:2*NN]
  P = minimize_divide(coord_x, coord_y)
  selection = P.plus_minus()
  if plot:
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

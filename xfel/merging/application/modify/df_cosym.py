from __future__ import absolute_import, division, print_function
from pandas import DataFrame as df
import pandas as pd
from scitbx.array_family import flex
from collections import OrderedDict
from itertools import permutations

class reconcile_cosym_reports:
  def __init__(self, reports):
    self.reports = reports
    # assume a list of reports, each report has (dataframe, coset_object)
    print("Processing reports from %d ranks"%(len(reports)))
    ireport = 0
    self.anchor_dataframe = reports[0][0]
    self.reference_partitions = reports[0][1].partitions
    self.result_dataframe = None
    self.n_cosets = len(self.reference_partitions)

  def simple_merge(self, voting_method="consensus"):
    print("simple merge")
    pd.set_option('display.max_rows', 300)
    pd.set_option('display.max_columns', 8)
    pd.set_option("display.width", None)
    # create as reports-as-lists data structure.
    # list { one item per rank }
    # item { one experiment-list per coset }
    reports_as_lists = []
    for idx, item in enumerate(self.reports):
      item_df = item[0]
      del(item_df["reindex_op"]) # do not need this column, coset is enough
      one_report = [[] for x in range(self.n_cosets)] # one sublist for every coset
      for lidx in range(len(item_df["experiment"])):
        one_report[item_df["coset"][lidx]].append(item_df["experiment"][lidx])
      reports_as_lists.append(one_report)

    # now modify the reports-as-list structure so ranks line up with respect to their coset assignments
    for idx in range(1, len(reports_as_lists)):
      base_report = reports_as_lists[idx-1]
      matches_vs_perm = []
      cache_permutations = list(permutations(reports_as_lists[idx]))
      for iperm,perm in enumerate(cache_permutations):
        matches_vs_perm.append(0)
        for icoset in range(len(perm)):
          imatches = len(set(base_report[icoset]).intersection(set(perm[icoset])))
          matches_vs_perm[-1] += imatches
          #print("Rank %d perm %d coset %d matches %d"%(idx,iperm,icoset,imatches))
      #print("matches for all permutations",matches_vs_perm, "choose", matches_vs_perm.index(max(matches_vs_perm)))
      # now choose the correct permutation and put it into reports-as-list
      correct_iperm = matches_vs_perm.index(max(matches_vs_perm))
      correct_perm = cache_permutations[correct_iperm]
      reports_as_lists[idx] = correct_perm

    # merge everything together into experiments plus votes
    experiment_lookup = dict()
    for idx in range(1, len(reports_as_lists)):
      for icoset in range(len(reports_as_lists[idx])):
        for uuid_expt in reports_as_lists[idx][icoset]:
          experiment_lookup[uuid_expt] = experiment_lookup.get( uuid_expt, len(experiment_lookup) )
    # a unique integer has now been assigned to every experiment uuid
    coset_assign = flex.size_t(flex.grid((len(experiment_lookup),3)))
    vote_matrix = flex.size_t(flex.grid((len(experiment_lookup),3)))
    vote_population = flex.size_t(len(experiment_lookup), 0)
    for idx in range(len(reports_as_lists)):
      for icoset in range(len(reports_as_lists[idx])):
        for uuid_expt in reports_as_lists[idx][icoset]:
          ekey = experiment_lookup[uuid_expt]
          coset_assign[ ekey, vote_population[ekey] ] = icoset
          vote_matrix[ ekey, vote_population[ekey] ] = idx
          vote_population[ekey]+=1

    #print out the matrices:
    if False:
      for item in experiment_lookup:
        ekey = experiment_lookup[item]
        print ("%4d"%ekey, item, " ", coset_assign[ekey,0], coset_assign[ekey,1], coset_assign[ekey,2],
                                 " ", vote_matrix[ekey,0], vote_matrix[ekey,1], vote_matrix[ekey,2])

    # tally the votes first.
    vote_tallies = flex.size_t(flex.grid((len(experiment_lookup),self.n_cosets)))
    for item in experiment_lookup:
      ekey = experiment_lookup[item]
      for ivote in range(3):
        vote_tallies[ ekey, coset_assign[ekey,ivote] ] += 1
      #print ("tallies %4d"%ekey, item, vote_tallies[ekey,0], vote_tallies[ekey,1]) #printout for two-coset case

    # choose winners, either by consensus or majority
    uuid_consensus = []
    uuid_majority = []
    coset_consensus = []
    coset_majority = []
    for item in experiment_lookup:
      ekey = experiment_lookup[item]
      for icoset in range(self.n_cosets):
        if vote_tallies[ekey,icoset]>=2:
          uuid_majority.append(item)
          coset_majority.append(icoset)
          if vote_tallies[ekey,icoset]==3: #hardcoded as it is assumed we always cast three votes
            uuid_consensus.append(item)
            coset_consensus.append(icoset)
          break
    print("Of %d experiments, majority cosets were chosen for %d and consensus for %d"%(
           len(experiment_lookup),len(coset_majority),len(coset_consensus)
    ))
    if voting_method == "majority":
      keyval = [("experiment", uuid_majority), ("coset", coset_majority)]
    elif voting_method == "consensus":
      keyval = [("experiment", uuid_consensus), ("coset", coset_consensus)]
    raw = OrderedDict(keyval)
    return df(raw)

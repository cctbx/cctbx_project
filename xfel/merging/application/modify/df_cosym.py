from __future__ import absolute_import, division, print_function
from pandas import DataFrame as df
import pandas as pd
import pickle
from scitbx.array_family import flex
import numpy as np
from collections import OrderedDict
from itertools import permutations

def iop():
  with open("../merge_002/myfile","rb") as F:
    refdata = pickle.load(F)
  symops = list(refdata.keys())[1:]

  refdata["max"] = refdata[symops].max(axis=1)
  refdata["min"] = refdata[symops].min(axis=1)
  refdata['diff'] = refdata["max"] - refdata["min"]
  refdata["reindex_op"] = refdata[symops].idxmax(axis=1)

  with open("cosym_myfile","rb") as F:
    cosdata = pickle.load(F)

  #inner join
  merged_inner = pd.merge(left=refdata, right=cosdata, left_on='experiment', right_on='experiment')


  same = []
  diff_same = flex.double()
  diff_diff = flex.double()
  for idx in range(len(merged_inner["reindex_op_x"])):
    X,Y = merged_inner["reindex_op_x"][idx],merged_inner["reindex_op_y"][idx]
    #if (X == "h,k,l" and
    #    Y == "h,k,l") or (
    #    X == "h,-h-k,-l" and
    #    Y == "-h,-l,-k"):
    #  same.append(False)
    #else: same.append(True)
    if (X == Y):
      same.append(True)
      diff_same.append(merged_inner["diff"][idx])
    else:
      same.append(False)
      diff_diff.append(merged_inner["diff"][idx])
  merged_inner['same'] = same


  pd.set_option('display.max_rows', 300)
  # What's the size of the output data?
  print(merged_inner.shape)
  print(merged_inner)
  A = len(merged_inner["same"])
  B = len(merged_inner[merged_inner["same"]==True])
  print("%d of %d comparisons match, %0.2f%%"%(B,A,100.*B/A))
  print("with same answers, diff was ", flex.mean(diff_same))
  print("with opposite answers, diff was ", flex.mean(diff_diff))

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
  def merge_all(self):
    pd.set_option('display.max_rows', 300)
    pd.set_option('display.max_columns', 8)
    pd.set_option("display.width", None)
    for item in self.reports[1:]:
      item_dataframe = item[0]

      reindexed_dataframes = [] # for each report, index it over all possible cosets
      anchor_item_inner = pd.merge(left=self.anchor_dataframe, right=item_dataframe, left_on='experiment', right_on='experiment')
      agreement_scores = np.array([0]*self.n_cosets,dtype=int)
      coset_x = np.array(anchor_item_inner["coset_x"],dtype=int)
      for icos in range(self.n_cosets): # loop over each coset, score agreement with the anchor
        coset_y = (np.array(anchor_item_inner["coset_y"],dtype=int) + icos) % self.n_cosets
        anchor_item_inner["agree"] = (coset_x==coset_y)
        n_agree = anchor_item_inner["agree"].value_counts()[True]
        agreement_scores[icos]=n_agree
        #print("Agree:",n_agree)
        #print(anchor_item_inner)
      consensus_icoset = np.argmax(agreement_scores)
      print("The consensus is ",consensus_icoset)
      #reset the agreement column to consensus coset
      if True:
        coset_y = (np.array(anchor_item_inner["coset_y"],dtype=int) + consensus_icoset) % self.n_cosets
        anchor_item_inner["agree"] = (coset_x==coset_y)
      item_dataframe["coset"] = (np.array(item_dataframe["coset"],dtype=int) + consensus_icoset) % self.n_cosets
      item_dataframe["keep"] = np.array([True]*len(item_dataframe["coset"]),dtype=bool) # fill a dummy value, keep all the new elements
      # update the item_dataframe keep column with agree flags from anchor_item_inner.  Only keep experiments that agree.
      print(item_dataframe)
      print(anchor_item_inner)
      return
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
      print("one report")
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
          print("Rank %d perm %d coset %d matches %d"%(idx,iperm,icoset,imatches))
      print("matches for all permutations",matches_vs_perm, "choose", matches_vs_perm.index(max(matches_vs_perm)))
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
    for item in experiment_lookup:
      ekey = experiment_lookup[item]
      print ("%4d"%ekey, item, " ", coset_assign[ekey,0], coset_assign[ekey,1], coset_assign[ekey,2],
                               " ", vote_matrix[ekey,0], vote_matrix[ekey,1], vote_matrix[ekey,2])

    #now take a vote, either by consensus or majority
    # tally the votes first.
    vote_tallies = flex.size_t(flex.grid((len(experiment_lookup),self.n_cosets)))
    for item in experiment_lookup:
      ekey = experiment_lookup[item]
      for ivote in range(3):
        vote_tallies[ ekey, coset_assign[ekey,ivote] ] += 1
      print ("tallies %4d"%ekey, item, vote_tallies[ekey,0], vote_tallies[ekey,1]) #printout for two-coset case

    # choose winners
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

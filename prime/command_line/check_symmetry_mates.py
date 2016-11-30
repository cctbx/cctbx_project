from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.check_symmetry_mates
"""
Author      : Uervirojnangkoorn, M.
Created     : 10/19/2016
"""
import argparse
import os
import cPickle as pickle
from collections import Counter
from cctbx import miller
from cctbx.array_family import flex
from libtbx.easy_mp import pool_map

def main(pr_pickle_file):
  pres = pickle.load(open(pr_pickle_file, 'rb'))
  if pres is not None:
    #sort according to indices
    perm = pres.observations.sort_permutation(by_value="packed_indices")
    obs_asu = pres.observations.select(perm)
    partiality = pres.partiality.select(perm)
    #correct to full reflections
    obs_asu = obs_asu.customized_copy(data=obs_asu.data()/partiality,
        sigmas=obs_asu.sigmas()/partiality)
    #group by similar indices
    obs_uniq = obs_asu.merge_equivalents().array()
    matches_uniq = miller.match_multi_indices(
        miller_indices_unique=obs_uniq.indices(),
        miller_indices=obs_asu.indices())
    pair_0 = flex.int([pair[0] for pair in matches_uniq.pairs()])
    pair_1 = flex.int([pair[1] for pair in matches_uniq.pairs()])
    group_id_list = flex.int([pair_0[pair_1[i]] for i in range(len(matches_uniq.pairs()))])
    tally = Counter()
    for elem in group_id_list:
      tally[elem] += 1
    #select only tokens with count > 1
    poly_tokens = [k for k,v in tally.iteritems() if v > 1]
    delta_I = flex.double()
    obs_I = flex.double()
    for token in poly_tokens:
      obs_group = obs_asu.select(obs_asu.indices() == obs_uniq.indices()[token])
      I_avg = flex.mean(obs_group.data())
      delta_I.extend(flex.double([abs(I-I_avg) for I in obs_group.data()]))
      obs_I.extend(obs_group.data())
      #if pres.pickle_filename.endswith('int_monarin_1_01380.pickle'):
      #  for ind, I, sigI, dI in zip(obs_group.indices(), obs_group.data(), obs_group.sigmas(), flex.abs(obs_group.data()-I_avg)):
      #    print ind, I, sigI, dI
    if len(obs_I) > 0:
      print pres.pickle_filename, '%6.2f'%(flex.sum(delta_I)*100/flex.sum(flex.abs(obs_I)))


if __name__ == "__main__":
  parser  = argparse.ArgumentParser(
      description = """Read post-refined pickles (.o file)
      then output correlation between symmetry mates on each image."""
  )
  parser.add_argument(
      'data',
      metavar='DATA',
      help='Path to post-refined pickles'
  )
  args = parser.parse_args()
  #build frame object for mproc
  n_proc = 16
  frames = []
  for pr_pickle_file in os.listdir(args.data):
    if pr_pickle_file.endswith('.o'):
      frames.append(os.path.join(args.data, pr_pickle_file))
  cc_results = pool_map(
      iterable=frames,
      func=main,
      processes=n_proc)

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME prime.compare_solution_pickles

from six.moves import cPickle as pickle
import argparse

def main(sol_fname, ind_fname):
  sol_pickle = pickle.load(open(sol_fname, "rb"))
  ind_pickle = pickle.load(open(ind_fname, "rb"))
  cn_match = 0
  for key in sol_pickle:
    if key in ind_pickle:
      if sol_pickle[key] == ind_pickle[key]:
        cn_match += 1
      else:
        print(key, sol_pickle[key], ind_pickle[key])
  print('Found %d images with %d matches'%(len(sol_pickle), cn_match))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(
    description='Compare two indexing ambiguity solution pickles'
  )
  parser.add_argument(
    'pickle1',
    metavar='First Pickle',
    help='Path to the first solution pickle'
  )
  parser.add_argument(
    'pickle2',
    metavar='Second pickle',
    help='Path to the second solution pickle'
  )
  args = parser.parse_args()
  print("Compare two solution pickles.")
  main(args.pickle1, args.pickle2)

from __future__ import absolute_import, division, print_function

from libtbx import easy_pickle
import os
import sys

def run(args):
  assert (len(args) > 0)
  assert os.path.isfile(args[0])
  pickle_file = args[0]
  func = easy_pickle.load(pickle_file)
  return func()

if (__name__ == "__main__"):
  run(sys.argv[1:])

from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.reorder_by_known_angles

import cPickle as pickle
import argparse

def read_pickles(data):
  frame_files = []
  for p in data:
    if os.path.isdir(p) == False:
      #check if list-of-pickle text file is given
      pickle_list_file = open(p,'r')
      pickle_list = pickle_list_file.read().split("\n")
      for pickle_filename in pickle_list:
        if os.path.isfile(pickle_filename):
          frame_files.append(pickle_filename)
    else:
      for pickle_filename in os.listdir(p):
        if pickle_filename.endswith('.pickle'):
          frame_files.append(p+'/'+pickle_filename)
  #check if pickle_dir is given in input file instead of from cmd arguments.
  if len(frame_files)==0:
    print 'No pickle files found.'
    exit()
  return sorted(frame_files)

def main(data, step):
  """
  Input: path pointing to integration pickles arranged in sequencial order
         step angle
  Function: reassign the calculated convention to next sequencial
  """

if __name__ == "__main__":
  parser = argparse.ArgumentParser(
    description='Reorder indexing orientation according to know angles'
  )
  parser.add_argument(
    'data',
    metavar='DATA',
    help='path to integration pickles'
  )
  parser.add_argument(
    'step',
    metavar='STEP',
    type=float,
    help='step angle (degrees)'
  )
  args = parser.parse_args()
  print args.data, args.step

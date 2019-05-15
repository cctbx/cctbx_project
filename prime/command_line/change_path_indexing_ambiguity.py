# LIBTBX_SET_DISPATCHER_NAME prime.change_path_indexing_ambiguity
"""
Author      : Uervirojnangkoorn, M.
Created     : 8/20/2015
Description : read indexing_ambiguity pickle and overwrite oldpath with new path
"""
from __future__ import absolute_import, division, print_function
import sys
from six.moves import cPickle as pickle
from six.moves import range

def read_input(args):
  data = ''
  oldpath = ''
  newpath = ''
  for i in range(len(args)):
    if args[i]=='-h':
      print(txt_help)
      exit()

    pair=args[i].split('=')
    if len(pair) == 2:
      if pair[0]=='data':
        data = pair[1]
      elif pair[0]=='oldpath':
        oldpath = pair[1]
      elif pair[0]=='newpath':
        newpath = pair[1]

  if data == '' or oldpath == '' or newpath == '':
    print("Please all parameters data, oldpath, and newpath.")
    exit()

  return data, oldpath, newpath


if (__name__ == "__main__"):
  #Help message
  txt_help = "Use this command to change cctbx.xfel integration pickle path to a new path\n"
  txt_help += "in case they have been moved to a new location.\n"
  txt_help += "Usage: prime.change_path_indexing_ambiguity data=indexing_ambiguity.pickle oldpath=/old/path/to/pickles newpath=/new/path/to/pickles\n"
  txt_help += "Good luck!\n"

  #Read input parameters and frames (pickle files)
  if len(sys.argv) == 1:
    print(txt_help)
    exit()

  data, oldpath, newpath = read_input(args = sys.argv[1:])

  pickle_data = pickle.load(open(data, "rb"))
  new_pickle_data = {}
  cn_i = 0
  for key in pickle_data.keys():
    basis = pickle_data[key]
    newkey = key.replace(oldpath, newpath)
    new_pickle_data[newkey] = basis
    print(cn_i+1, newkey, basis)
    cn_i +=1

  pickle.dump(new_pickle_data, open("new_indexing_ambiguity.pickle", "wb"))
  print('Found and replace %6d keys'%(cn_i))

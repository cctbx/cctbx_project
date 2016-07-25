from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime._genref_determine_mean_I
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : prime._genref_determine_mean_I is an internal command to suppport
queuing system and mproc.
"""
from libtbx.easy_mp import pool_map
import sys
import cPickle as pickle
from prime.postrefine import postref_handler

def determine_mean_I_mproc(args):
  frame_no, pickle_filename, iparams = args
  prh = postref_handler()
  mean_I, txt_out = prh.calc_mean_intensity(pickle_filename, iparams)
  print frame_no, pickle_filename, mean_I
  if mean_I is not None:
    pickle.dump(mean_I, open(iparams.run_no+'/pickles/'+str(frame_no)+".o","wb"),pickle.HIGHEST_PROTOCOL)

if (__name__ == "__main__"):
  if len(sys.argv)==1:
    print 'Not allowed. This is an internal command for queuing system.'
    exit()
  #load input
  inp_pickle = pickle.load(open(sys.argv[1], "rb"))
  iparams = inp_pickle['iparams']
  frames = inp_pickle['frames']
  #calculate mean of mean I
  mm_I = 0
  determine_mean_I_result = pool_map(
          iterable=frames,
          func=determine_mean_I_mproc,
          processes=iparams.n_processors)

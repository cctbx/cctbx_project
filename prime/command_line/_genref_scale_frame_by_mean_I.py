from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime._genref_scale_frame_by_mean_I
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : prime._genref_scale_frame_by_mean_I is an internal command to
suppport queuing system and mproc.
"""
from libtbx.easy_mp import pool_map
import sys
import cPickle as pickle
from prime.postrefine import postref_handler

def scale_frame_by_mean_I_mproc(args):
  frame_no, pickle_filename, iparams, mm_I = args
  prh = postref_handler()
  pres, txt_out = prh.scale_frame_by_mean_I(frame_no, pickle_filename, iparams, mm_I)
  pickle.dump(pres, open(iparams.run_no+'/pickles/'+str(frame_no)+".o","wb"),pickle.HIGHEST_PROTOCOL)

if (__name__ == "__main__"):
  if len(sys.argv)==1:
    print 'Not allowed. This is an internal command for queuing system.'
    exit()
  #load input
  inp_pickle = pickle.load(open(sys.argv[1], "rb"))
  iparams = inp_pickle['iparams']
  frames = inp_pickle['frames']
  #perform initial scaling
  scale_frame_by_mean_I_result = pool_map(
          iterable=frames,
          func=scale_frame_by_mean_I_mproc,
          processes=iparams.n_processors)

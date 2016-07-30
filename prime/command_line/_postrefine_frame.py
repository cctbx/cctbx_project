from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime._postrefine_frame
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : prime._postrefine_frame is an internal command to
suppport queuing system and mproc.
"""
from libtbx.easy_mp import pool_map
import sys
import cPickle as pickle
from prime.postrefine import postref_handler

def postrefine_by_frame_mproc(args):
  frame_no, pickle_result, iparams, miller_array_ref, avg_mode = args
  prh = postref_handler()
  pres, txt_out = prh.postrefine_by_frame(frame_no, pickle_result, iparams, miller_array_ref, avg_mode)
  #overwrite to the original frame_no file.
  pickle.dump(pres, open(iparams.run_no+'/pickles/'+str(frame_no)+".o","wb"),pickle.HIGHEST_PROTOCOL)

if (__name__ == "__main__"):
  if len(sys.argv)==1:
    print 'Not allowed. This is an internal command for queuing system.'
    exit()
  #load input
  inp_pickle = pickle.load(open(sys.argv[1], "rb"))
  iparams = inp_pickle['iparams']
  frames = inp_pickle['frames']
  #perform post-refinement
  postrefine_by_frame_result = pool_map(
            iterable=frames,
            func=postrefine_by_frame_mproc,
            processes=iparams.n_processors)

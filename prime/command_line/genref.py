from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.genref
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : prime.genref generates scaled observations and writes out pickles.
* suppport queuing system and mproc.
"""
import os, sys
import numpy as np
import cPickle as pickle
from subprocess import call
from prime.postrefine.mod_run import run_handler

def calc_mean_of_mean_I(iparams):
  frames_mean_I = [pickle.load( open( iparams.run_no+'/pickles/'+fname, "rb" ) ) \
    for fname in os.listdir(iparams.run_no+'/pickles/')]
  mean_of_mean_I = np.median(frames_mean_I)
  return mean_of_mean_I

class genref_handler(object):
  """
  A wrapper class for genref handler
  """
  def __init__(self):
    """
    Intialitze parameters
    """

  def run(self, args):
    #read inputs
    from prime.postrefine.mod_input import process_input
    iparams, txt_out_input = process_input(args)
    print txt_out_input
    self.run_by_params(iparams)

  def run_by_params(self, iparams):
    runh = run_handler()
    #read all integration pickles
    from prime.postrefine.mod_input import read_pickles
    frame_files = read_pickles(iparams.data)
    n_frames = len(frame_files)
    if n_frames == 0:
      print "No integration pickle found. Exit program."
      exit()
    frames = [(i, frame_files[i], iparams) for i in range(n_frames)]
    mm_I = 0
    #run command to calculate mean_I
    if iparams.flag_apply_b_by_frame == False:
      inp_pickle = {'iparams':iparams, 'frames':frames}
      pickle.dump(inp_pickle, open(iparams.run_no+'/inputs/0.inp',"wb"))
      call(["prime._genref_determine_mean_I", iparams.run_no+'/inputs/0.inp'])
      runh.check_done(iparams, n_frames)
      mm_I = calc_mean_of_mean_I(iparams)
    #run command for scaling
    if iparams.queue.mode is None:
      #run single node
      frames = [(i, frame_files[i], iparams, mm_I) for i in range(n_frames)]
      inp_pickle = {'iparams':iparams, 'frames':frames}
      pickle.dump(inp_pickle, open(iparams.run_no+'/inputs/0.inp',"wb"))
      call(["prime._genref_scale_frame_by_mean_I", iparams.run_no+'/inputs/0.inp'])
    else:
      #run on n_nodes
      n_imgs_per_node = int(round(n_frames/iparams.queue.n_nodes))
      for i_node in range(iparams.queue.n_nodes):
        start_frame = i_node*n_imgs_per_node
        if i_node < iparams.queue.n_nodes - 1:
          end_frame = start_frame + n_imgs_per_node
        else:
          end_frame = n_frames
        frames = [(i, frame_files[i], iparams, mm_I) for i in range(start_frame, end_frame)]
        inp_pickle = {'iparams':iparams, 'frames':frames}
        pickle.dump(inp_pickle, open(iparams.run_no+'/inputs/'+str(i_node)+'.inp',"wb"))
        call(["bsub","-q",iparams.queue.qname,"-o",iparams.run_no+"/qout/qout_gr.txt","prime._genref_scale_frame_by_mean_I", iparams.run_no+"/inputs/"+str(i_node)+".inp"])
    runh.check_done(iparams, n_frames)
    #write output to logfile
    txt_out = 'Scaling complete. Run prime.merge your_input_phil.phil to merge for a reference set.\n'
    print txt_out
    f = open(iparams.run_no+'/log.txt', 'a')
    f.write(txt_out)
    f.close()


if __name__=="__main__":
  grh = genref_handler()
  grh.run(sys.argv[1:] if len(sys.argv)>1 else None)

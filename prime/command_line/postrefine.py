from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.postrefine
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Description : Commands linked to prime.postrefine libraries.
'''
import os, sys
from prime.postrefine.mod_mx import mx_handler
import cPickle as pickle
from subprocess import call
from prime.postrefine.mod_run import run_handler

class postrefine_handler(object):
  def __init__(self):
    """
    Intialize parameters
    """

  def run(self, args):
    #read inputs
    runh = run_handler()
    from prime.postrefine.mod_input import process_input
    iparams, txt_out_input = process_input(argv=args, flag_check_exist=False)
    iparams.flag_volume_correction = False
    if iparams.partiality_model == "Lognormal":
      iparams.voigt_nu = 0.008 #use voigt_nu as lognpdf zero parameter
    #read all result pickles
    try:
      DIR = iparams.run_no+'/pickles/'
      pickle_results = [pickle.load(open(DIR+fname, "rb")) for fname in os.listdir(DIR)]
      file_no_results = [int(fname.split('.')[0]) for fname in os.listdir(DIR)]
      n_results = len(pickle_results)
    except Exception:
      print "Error reading input pickles."
      print "*VERSION UPGRADE NOTE* use prime.run instead of prime.postrefine to run all processes together."
      exit()
    #get reference file - look for n.mtz with n as maximum number.
    hklrefin = None
    if iparams.hklrefin is None:
      DIR = iparams.run_no+'/mtz/'
      file_no_list = [int(fname.split('.')[0]) for fname in os.listdir(DIR)]
      if len(file_no_list) > 0:
        hklrefin = DIR + str(max(file_no_list)) + '.mtz'
    else:
      hklrefin = iparams.hklrefin
    if hklrefin is None:
      print "No reference set found. Exit program"
    print "Reference set:", hklrefin, " No. of images:", n_results
    mxh = mx_handler()
    flag_hklrefin_found, miller_array_ref = mxh.get_miller_array_from_reflection_file(hklrefin)
    #post-refinement
    avg_mode = 'weighted'
    #run command for post-refinement
    if iparams.queue.mode is None:
      frames = [(file_no_results[i], pickle_results[i], iparams, miller_array_ref, avg_mode) for i in range(n_results)]
      inp_pickle = {'iparams':iparams, 'frames':frames}
      pickle.dump(inp_pickle, open(iparams.run_no+'/inputs/0.inp',"wb"))
      call(["prime._postrefine_frame", iparams.run_no+'/inputs/0.inp'])
    else:
      #run on n_nodes
      n_imgs_per_node = int(round(n_results/iparams.queue.n_nodes))
      for i_node in range(iparams.queue.n_nodes):
        start_frame = i_node*n_imgs_per_node
        if i_node < iparams.queue.n_nodes - 1:
          end_frame = start_frame + n_imgs_per_node
        else:
          end_frame = n_results
        frames = [(i, pickle_results[i], iparams, miller_array_ref, avg_mode) for i in range(start_frame, end_frame)]
        inp_pickle = {'iparams':iparams, 'frames':frames}
        pickle.dump(inp_pickle, open(iparams.run_no+'/inputs/'+str(i_node)+'.inp',"wb"))
        call(["bsub","-q",iparams.queue.qname,"-o",iparams.run_no+"/qout/qout_pr.txt","prime._postrefine_frame", iparams.run_no+"/inputs/"+str(i_node)+".inp"])
    runh.check_done(iparams, n_results)
    print "Post-refinement completed. Run prime.merge for the merged reflection file."

if __name__=="__main__":
  prh = postrefine_handler()
  prh.run(sys.argv[:1])

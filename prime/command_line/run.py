from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.run
"""
Author      : Uervirojnangkoorn, M.
Created     : 7/14/2016
Description : Commands for running prime main class.
"""
from prime.command_line.genref import genref_handler
from prime.command_line.postrefine import postrefine_handler
from prime.command_line.merge import merge_handler
from prime.postrefine.mod_input import process_input
import os, sys, time, logging, shutil

class run_handler(object):
  def __init__(self):
    """
    Intialize parameters
    """

  def run(self, args):
    #generate reference set, prepare scaled pickles
    print "Scaling integration pickles"
    grh = genref_handler()
    grh.run(args)
    #merge for the first reference set
    print "Merging for a reference set."
    mh = merge_handler()
    mh.run(args)
    #start post-refinement loops
    iparams, txt_out_input = process_input(argv=args, flag_check_exist=False)
    prh = postrefine_handler()
    for i_cycle in range(iparams.n_postref_cycle):
      print "Post-refinement cycle ", i_cycle+1
      prh.run(args)
      print "Merging cycle ", i_cycle+1
      if i_cycle < iparams.n_postref_cycle - 1:
        mh.run(args)
      else:
        #final run
        mh.run(args, avg_mode='final')
        #copy the max_no.mtz to postref_final.mtz
        DIR = iparams.run_no+'/mtz/'
        file_no_list = [int(fname.split('.')[0]) for fname in os.listdir(DIR)]
        if len(file_no_list) > 0:
          mtz_final_fname = DIR + str(max(file_no_list)) + '.mtz'
          shutil.copy(mtz_final_fname, iparams.run_no+'/postref_final.mtz')

if __name__=="__main__":
  #capture starting time
  program_starts= time.time()
  logging.captureWarnings(True)
  formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
  console_handler = logging.StreamHandler()
  console_handler.setLevel(logging.ERROR)
  console_handler.setFormatter(formatter)
  logging.getLogger().addHandler(console_handler)
  logging.getLogger('py.warnings').addHandler(console_handler)
  logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', level=logging.DEBUG)
  rh = run_handler()
  rh.run(sys.argv[:1])
  print 'Elapsed times %10.1f seconds'%(time.time() - program_starts)

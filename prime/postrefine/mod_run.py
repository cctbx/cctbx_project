from __future__ import division
import time, os

class run_handler(object):
  """
  Author      : Uervirojnangkoorn, M.
  Created     : 4/15/2016
  A collection of functions used in scaling up/ parallelizing prime.
  """
  def __init__(self):
    """
    Constructor
    """

  def check_done(self, iparams, n_frames):
    #only allow the module to continue if all frames are found or time out.
    program_starts = time.time()
    while((time.time() - program_starts < iparams.timeout_seconds) and \
        (len([fname for fname in os.listdir(iparams.run_no+'/pickles/')]) < n_frames)):
      print("Found {0} frames - Elapsed times: {1} seconds".format(len([fname for fname in os.listdir(iparams.run_no+'/pickles/')]), time.time() - program_starts))
    print "%10.0f frames processed."%(len([fname for fname in os.listdir(iparams.run_no+'/pickles/')]))

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
  def count_frames(self):
    #only count files that are older than the program
    if self.iparams.queue.mode is None:
      return len([fname for fname in os.listdir(self.iparams.run_no+'/pickles/')])
    else:
      return len([fname for fname in os.listdir(self.iparams.run_no+'/pickles/') if os.path.getmtime(self.iparams.run_no+'/pickles/'+fname)-self.pg_starts > 0])

  def check_done(self, iparams, n_frames):
    self.iparams = iparams
    self.pg_starts = time.time()
    #only allow the module to continue if all frames are found or time out.
    n_frames_now = self.count_frames()
    while((time.time() - self.pg_starts < iparams.timeout_seconds) and (self.count_frames() < n_frames)):
      n_frames_tmp = self.count_frames()
      if n_frames_tmp > n_frames_now:
        n_frames_now = n_frames_tmp
        print "%10.0f frames processed %10.1f seconds"%(n_frames_now, time.time() - self.pg_starts)
    print "%10.0f frames processed."%(self.count_frames())

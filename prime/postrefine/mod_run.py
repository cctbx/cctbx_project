from __future__ import division, print_function
""" handle queue usage """
__author__ = 'Monarin Uervirojnangkoorn, monarin@gmail.com'

import time, os

class run_handler(object):

  def check_done(self, iparams):
    #only allow the module to continue if all frames are found or time out.
    program_starts = time.time()
    while((time.time() - program_starts < iparams.timeout_seconds) and \
        not os.path.isfile(os.path.join(iparams.run_no,'.done'))):
      print("Running - Elapsed times: {0:6.1f} seconds".format(time.time() - program_starts))
      time.sleep(5)
    print("Done. You can view your results and log file in ", iparams.run_no)

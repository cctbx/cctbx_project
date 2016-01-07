from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME iota.bulletproof

'''
Author      : Brewster, A.S., Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 11/06/2015
Description : Script for indexing using run_one_index_core. Designed to be
              called by easy_run.fully_buffered, protecting the calling program
              from boost errors.
'''

import sys
from xfel.phil_preferences import load_cxi_phil
from xfel.cxi.display_spots import run_one_index_core
from libtbx import easy_pickle


if __name__ == "__main__":
  # should be invoked like this: "iota.bulletproof tmppath target args"
  tmppath = sys.argv[1]
  target = sys.argv[2]
  args = sys.argv[3:]

  try:
    # index the image
    horizons_phil = load_cxi_phil(target, args)
    info = run_one_index_core(horizons_phil)
    # save specific results from the info object to be used by iota
    int_final = info.last_saved_best
    easy_pickle.dump(tmppath, int_final)
  except Exception, e:
    if hasattr(e, "classname"):
      error_message = "{}: {}".format(e.classname, e[0].replace('\n',' ')[:50])
    else:
      error_message = "{}".format(str(e).replace('\n', ' ')[:50])
    print e
    # save the error message to be picked up by iota
    easy_pickle.dump(tmppath, error_message)

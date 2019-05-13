from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.run
""" handle prime run """
from __future__ import print_function
__author__ = 'Monarin Uervirojnangkoorn, monarin@gmail.com'

from subprocess import call
from prime.postrefine.mod_run import run_handler
from prime.postrefine.mod_input import process_input
import sys

if __name__ == "__main__":
  iparams, txt_out = process_input(sys.argv[1:] if len(sys.argv) > 1 else None, flag_mkdir=False)
  if iparams.queue.mode:
    args = ["bsub","-q",iparams.queue.qname,"-n", str(iparams.queue.n_nodes), "prime.postrefine"]+sys.argv[1:] if len(sys.argv) > 1 else []
    call(args)
    print("Submitting prime job to ", iparams.queue.qname)
    runh = run_handler()
    runh.check_done(iparams)
  else:
    import prime.command_line.postrefine
    prime.command_line.postrefine.run(sys.argv[1:])

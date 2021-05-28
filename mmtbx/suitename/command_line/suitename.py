# LIBTBX_SET_DISPATCHER_NAME phenix.suitename
# LIBTBX_SET_DISPATCHER_NAME molprobity.suitename

#import libtbx.load_env
#from libtbx.utils import Usage

import os,sys,inspect
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir) 

from suitename import main  # from suitename.py in PARENT directory

def run(args):
  main()

if (__name__ == "__main__"):
  run(sys.argv[1:])

# LIBTBX_SET_DISPATCHER_NAME phenix.suitename
# LIBTBX_SET_DISPATCHER_NAME molprobity.suitename
# LIBTBX_SET_DISPATCHER_NAME cctbx.suitename

from __future__ import absolute_import, division, print_function
<<<<<<< HEAD

from iotbx.cli_parser import run_program
from mmtbx.programs import suitename

if __name__ == '__main__':
  run_program(args=sys.argv[1:], program_class=suitename.Program)
=======
import sys
from mmtbx.programs import suitename

if __name__ == '__main__':
  print("-----------------command line reporting-----------------\n")
  suitename.run(args=sys.argv[1:])
  
>>>>>>> +suitename

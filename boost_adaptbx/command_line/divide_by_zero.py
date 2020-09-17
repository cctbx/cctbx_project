from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.divide_by_zero

import boost_adaptbx.boost.python as bp
import sys

def run(args):
  assert len(args) <= 1
  if len(args) == 0:
    print("Now dividing by zero (in C++) ...")
    sys.stdout.flush()
    result = bp.ext.divide_doubles(1, 0)
    print("Result:", result)
  else:
    bp.floating_point_exceptions.division_by_zero_trapped = False
    print("Dividing by zero in C++: not gonna be caught")
    result = bp.ext.divide_doubles(1, 0)
    print("Result:", result)
    bp.floating_point_exceptions.division_by_zero_trapped = True
    print()
    print("Dividing by zero in C++: gonna crash")
    sys.stdout.flush()
    result = bp.ext.divide_doubles(1, 0)
    print("Result:", result)

if (__name__ == "__main__"):
  run(sys.argv[1:])

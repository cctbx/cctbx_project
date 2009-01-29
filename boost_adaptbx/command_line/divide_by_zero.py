# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.divide_by_zero

import boost.python
import sys

def run(args):
  assert len(args) <= 1
  if len(args) == 0:
    print "Now dividing by zero (in C++) ..."
    sys.stdout.flush()
    result = boost.python.ext.divide_doubles(1, 0)
    print "Result:", result
  else:
    boost.python.floating_point_exceptions.division_by_zero_trapped = False
    print "Dividing by zero in C++: not gonna be caught"
    result = boost.python.ext.divide_doubles(1, 0)
    print "Result:", result
    boost.python.floating_point_exceptions.division_by_zero_trapped = True
    print
    print "Dividing by zero in C++: gonna crash"
    sys.stdout.flush()
    result = boost.python.ext.divide_doubles(1, 0)
    print "Result:", result

if (__name__ == "__main__"):
  run(sys.argv[1:])

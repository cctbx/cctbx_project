from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME boost_adaptbx.nested_cpp_loops_with_check_signals

def run(args):
  assert len(args) == 2, "iterations_outer, iterations_inner"
  import boost_adaptbx.boost.python as bp
  count = bp.ext.nested_cpp_loops_with_check_signals(
    iterations_outer=int(args[0]),
    iterations_inner=int(args[1]))
  print("actual iterations outer:", count)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

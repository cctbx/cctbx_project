def run(args):
  assert len(args) == 0
  from libtbx.utils import show_string
  from libtbx import easy_run
  import libtbx.load_env
  import os
  op = os.path
  test3 = libtbx.env.under_dist(
    module_name="cbflib",
    path="pycbf/pycbf_test3.py")
  assert op.isfile(test3)
  cmd = "cbflib.python %s" % show_string(test3)
  print cmd
  easy_run.call(command=cmd)

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])

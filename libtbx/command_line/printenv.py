from __future__ import absolute_import, division, print_function
import sys, os

def show(out):
  var_names = sorted(os.environ.keys())
  for var_name in var_names:
    print("%s=%s" % (var_name, os.environ[var_name]), file=out)

def run(args):
  assert len(args) == 0
  show(out=sys.stdout)

if (__name__ == "__main__"):
  run(args=sys.argv[1:])

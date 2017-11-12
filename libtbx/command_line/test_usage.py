from __future__ import absolute_import, division, print_function

import sys

def run(args):
  from libtbx.test_utils import test_usage
  from libtbx.utils import Usage
  if not args:
    raise Usage("""\
libtbx.test_usage command

Utility to verify sensible behavior of programs when run without arguments.
""")
  cmd = args[0]
  return test_usage(cmd)

if __name__ == "__main__":
  run(sys.argv[1:])


from libtbx.test_utils import test_usage
from libtbx.utils import Usage
import sys

def run (args) :
  if (len(args) == 0) :
    raise Usage("""\
libtbx.test_usage command

Utility to verify sensible behavior of programs when run without arguments.
""")
  cmd = args[0]
  return test_usage(cmd)

if (__name__ == "__main__") :
  run(sys.argv[1:])

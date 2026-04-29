"""
tests.tst_basic
"""

from __future__ import absolute_import, division, print_function
import sys
from libtbx import easy_run

def run(args):
  ret = easy_run.call("scitbx.unicode_examples")
  ret += easy_run.call("scitbx.show_sizes")
  ret += easy_run.call("scitbx.show_exp_times 100")
  ret += easy_run.call("echo OK")
  return ret

if __name__ == "__main__":
  import sys
  run(sys.argv[1:])

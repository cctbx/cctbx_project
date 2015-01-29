# LIBTBX_SET_DISPATCHER_NAME cctbx_regression.test_nightly

from __future__ import division
from libtbx.command_line import run_tests_parallel
import sys

if (__name__ == "__main__") :
  args = [
    "module=libtbx",
    "module=boost_adaptbx",
    "module=scitbx",
    "module=cctbx",
    "module=iotbx",
    "module=mmtbx",
    "module=smtbx",
    "nproc=Auto",
  ]
  if (run_tests_parallel.run(args) > 0) :
    sys.exit(1)

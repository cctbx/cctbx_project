from __future__ import absolute_import, division, print_function
from libtbx.test_utils.pytest import discover

tst_list = [
  ["$D/test/test_cout_compile.py", "stop"],
] + discover()

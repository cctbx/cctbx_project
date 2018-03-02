from __future__ import absolute_import, division, print_function
from libtbx.test_utils.pytest import discover

tst_list = [
  "$D/tst_ext.py",
  "$D/tst_equivalence.py",
  "$D/tst_read.py",
  "$D/tst_cout.py",
  ["$D/tst_cout_compile.py", "stop"],
] + discover()

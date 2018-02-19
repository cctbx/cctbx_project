from __future__ import absolute_import, division, print_function
from libtbx.test_utils.pytest import discover

tst_list = [
  "$D/tst_ext.py",
  "$D/tst_equivalence.py",
  "$D/tst_read.py",
  "$D/tst_cout.py",
  "$D/test/tst_show_calls.py",
  "$D/test/tst_command_line.py",
  "$D/test/tst_separate_files.py",
  ["$D/tst_cout_compile.py", "stop"],
] + discover()

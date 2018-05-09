from __future__ import absolute_import, division, print_function

from libtbx.test_utils.pytest import discover

tst_list = discover()

# To write tests for dxtbx:

# 1. Test file should be named test_*.py
# 2. Test methods should be named test_*()
# 3. Nothing else needed. Rest happens by magic.

# To run dxtbx tests:

# run 'pytest' inside dxtbx directory

# For more information see:
#   https://github.com/dials/dials/wiki/pytest

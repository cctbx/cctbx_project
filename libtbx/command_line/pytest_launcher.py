# LIBTBX_SET_DISPATCHER_NAME pytest
# LIBTBX_SET_DISPATCHER_NAME py.test
from __future__ import absolute_import, division, print_function

import sys

import pytest

# modify sys.argv so the command line help shows the right executable name
sys.argv[0] = 'pytest'

sys.exit(pytest.main())

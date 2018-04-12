#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import os

import pytest

@pytest.fixture(scope="session")
def dials_regression():
  '''Return the absolute path to the dials_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    import dials_regression as dr
  except ImportError:
    pytest.skip("dials_regression required for this test")
  return os.path.dirname(dr.__file__)

from libtbx.test_utils.pytest import libtbx_collector
pytest_collect_file = libtbx_collector()

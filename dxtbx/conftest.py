#
# See https://github.com/dials/dials/wiki/pytest for documentation on how to
# write and run pytest tests, and an overview of the available features.
#

from __future__ import absolute_import, division, print_function

import libtbx.load_env
import pytest

@pytest.fixture
def dials_regression():
  '''Return the absolute path to the dials_regression module as a string.
     Skip the test if dials_regression is not installed.'''
  try:
    return libtbx.env.dist_path('dials_regression')
  except KeyError:
    pytest.skip("dials_regression required for this test")

from libtbx.test_utils.pytest import libtbx_collector
pytest_collect_file = libtbx_collector()

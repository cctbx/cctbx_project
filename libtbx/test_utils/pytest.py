from __future__ import division, absolute_import
import atexit
import libtbx.load_env
import os

def discover(module, pytestargs=None):
  '''
  pytest compatibility layer

  This function discovers pytest tests in a module directory so that they can
  be run via libtbx.run_tests_parallel without being named explicitly in the
  module's run_tests.py.

  To use this function, add
    from libtbx.test_utils.pytest import discover
    tst_list = tst_list + discover( $module )
  to your run_tests.py, substituting '$module' with the name of your module
  directory.
  You can pass further arguments to the pytest discovery process as a list in
  the pytestargs parameter.

  The function tests for the presence of pytest and mock, and generates a
  helpful warning if either is missing. You can skip the discovery of tests
  by setting the environment variable LIBTBX_SKIP_PYTEST.
  '''

  if 'LIBTBX_SKIP_PYTEST' in os.environ:
    return []

  if pytestargs is None:
    pytestargs = []

  try:
    import pytest
    import mock
  except ImportError:
    def pytest_warning():
      print "=" * 55
      print " WARNING: Skipping some tests for %s\n" % module
      print " To run all available tests you need to install pytest"
      print " and the python mocking package, e.g. by running\n"
      print "    libtbx.python -m pip install pytest mock"
      print "=" * 55
    pytest_warning()
    atexit.register(pytest_warning)
    return []

  print "Discovering pytest tests for %s:" % module
  test_list = []
  dist_dir = libtbx.env.dist_path(module)
  class TestDiscoveryPlugin:
    def pytest_itemcollected(self, item):
      test_list.append([ "libtbx.python", "-m", "pytest", '--noconftest', '--basetemp=pytest',
        '"%s"' % os.path.join(dist_dir, item.nodeid) ])
  pytest.main(['-qq', '--collect-only', '--noconftest', dist_dir] + pytestargs, plugins=[TestDiscoveryPlugin()])
  return test_list


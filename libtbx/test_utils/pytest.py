from __future__ import absolute_import, division, print_function
import atexit
import os
import sys
import itertools

import libtbx.load_env
from libtbx.test_utils.parallel import run_command

_first_pytest_collection = True
_pytest_unique_counter = itertools.count(1)

def discover(module=None, pytestargs=None):
  '''
  pytest compatibility layer

  This function discovers pytest tests in a module directory so that they can
  be run via libtbx.run_tests_parallel without being named explicitly in the
  module's run_tests.py.

  To use this function, add
    from libtbx.test_utils.pytest import discover
    tst_list = tst_list + discover()
  to your run_tests.py.
  You can pass further arguments to the pytest discovery process as a list in
  the pytestargs parameter, and you can specify the cctbx module name if the
  automatic discovery does not work for any reason.

  The function tests for the presence of pytest and mock, and generates a
  helpful warning if either is missing. You can skip the discovery of tests
  by setting the environment variable LIBTBX_SKIP_PYTEST.
  '''

  if 'LIBTBX_SKIP_PYTEST' in os.environ:
    '''If LIBTBX_SKIP_PYTEST is set then the user is running libtbx testing but
       does not want any pytests to be run.
       Alternatively the user is running pytest, and pytest wants to find all
       legacy libtbx tests. In this case we also do not attempt to recursively
       find pytests.
    '''
    # Buck stops here in either case
    return []

  if pytestargs is None:
    pytestargs = []

  try:
    import pytest
    import mock
  except ImportError:
    def pytest_warning():
      print("=" * 55)
      print(" WARNING: Skipping some tests for %s\n" % module)
      print(" To run all available tests you need to install pytest")
      print(" and the python mocking package, e.g. by running\n")
      print("    libtbx.python -m pip install pytest mock")
      print("=" * 55)
    pytest_warning()
    atexit.register(pytest_warning)
    return []

  if not module:
    # Try to determine the name of the calling module, and thus the name of the
    # cctbx module. Use exception trick to pick up the current frame.
    try:
      raise Exception()
    except Exception:
      frame = sys.exc_info()[2].tb_frame.f_back
    caller = frame.f_globals['__name__']
    if not caller.endswith('.run_tests'):
      raise RuntimeError('Only use discover() from within run_tests.py ' \
                       + 'or specify the module name manually.')
    module = caller[:-10]

  class L(list):
    """Subclass list so that it can accept additional attributes."""

  print("Discovering pytest tests for %s:" % module)
  test_list = []
  dist_dir = libtbx.env.dist_path(module)
  class TestDiscoveryPlugin:
    def pytest_itemcollected(self, item):
      global _pytest_unique_counter
      testarray = L([ "libtbx.python", "-m", "pytest", '-rsxX', '--basetemp=pytest%st%03d' % (os.path.sep, next(_pytest_unique_counter)),
        '"%s"' % (item.fspath + '::' + item.nodeid.split('::', 1)[1]) ])
      testclass = module + '.' + item.location[0].replace(os.path.sep, '.')
      if testclass.endswith('.py'):
        testclass = testclass[:-3]
      testarray.test_class = testclass
      testarray.test_name = item.name
      test_list.append(testarray)

  pytest_parameters = ['-qq', '--collect-only']

  # Only show pytest warnings during first collection
  global _first_pytest_collection
  if not _first_pytest_collection:
    pytest_parameters.append('--disable-pytest-warnings')
  _first_pytest_collection = False

  try:
    # Now set LIBTBX_SKIP_PYTEST so we can collect all pytests without pytest
    # recursively trying to find legacy libtbx tests.
    os.environ['LIBTBX_SKIP_PYTEST'] = "1"

    pytest.main(pytest_parameters + [ dist_dir ] + pytestargs, plugins=[TestDiscoveryPlugin()])
  finally:
    del os.environ['LIBTBX_SKIP_PYTEST']

  if test_list:
    # Ensure the common basetemp directory pytest/ exists
    try:
      os.mkdir('pytest')
    except OSError:
      pass

  return test_list

def libtbx_collector():
  '''
  libtbx compatibility layer:

  return a function that enables pytest to collect and run all libtbx legacy tests of a module

  To use this you need to add
    import libtbx.test_utils.pytest
    pytest_collect_file = libtbx.test_utils.pytest.libtbx_collector()
  to your conftest.py in the module root directory.
  '''

  import pytest

  class LibtbxTestException(Exception):
    '''Custom exception for error reporting.'''
    def __init__(self, stdout, stderr):
      self.stdout = stdout
      self.stderr = stderr

  class LibtbxTest(pytest.Item):
    def __init__(self, name, parent, test_command, test_parameters):
      super(LibtbxTest, self).__init__(name, parent)
      self.test_cmd = test_command
      if test_command.endswith('.py'):
        self.test_cmd = 'libtbx.python "%s"' % self.test_cmd
      self.test_parms = test_parameters
      self.full_cmd = ' '.join([self.test_cmd] + self.test_parms)
      if not hasattr(self, 'module'):
        self.module = None
      if not hasattr(self, '_fixtureinfo'):
        self._fixtureinfo = self.session._fixturemanager.getfixtureinfo(self, self.runtest, self)

    def runtest(self):
      rc = run_command(self.full_cmd)
      if rc is None:
        # run_command only returns None if CTRL+C pressed
        raise KeyboardInterrupt()
      self.add_report_section('call', 'stdout', '\n'.join(rc.stdout_lines))
      self.add_report_section('call', 'stderr', '\n'.join(rc.stderr_lines))
      if rc.stderr_lines or rc.return_code != 0:
        raise LibtbxTestException(rc.stdout_lines, rc.stderr_lines)

    def repr_failure(self, excinfo):
      '''called when self.runtest() raises an exception.'''
      if isinstance(excinfo.value, LibtbxTestException):
        return "\n".join(excinfo.value.stderr)

    def reportinfo(self):
      return self.fspath, 0, self.full_cmd

  def pytest_collect_file(path, parent):
    if 'LIBTBX_SKIP_PYTEST' in os.environ:
      '''The pytest discovery process is ran from within libtbx, so do not
         attempt to find libtbx legacy tests.'''
      return

    class LibtbxRunTestsFile(pytest.File):
      def collect(self):
        try:
          os.environ['LIBTBX_SKIP_PYTEST'] = '1'
          import importlib
          # Guess the module import path from the location of this file
          test_module = self.fspath.dirpath().basename
          # We must be directly inside the root of a configured module.
          # If this module isn't configured, then we don't want to run tests.
          if not libtbx.env.has_module(test_module):
            return
          run_tests_module = test_module + "." + self.fspath.purebasename
          run_tests = importlib.import_module(run_tests_module)
        finally:
          del os.environ['LIBTBX_SKIP_PYTEST']

        for test in run_tests.tst_list:
          from six import string_types 
          if isinstance(test, string_types):
            testfile = test
            testparams = []
            testname = 'main'
          else:
            testfile = test[0]
            testparams = [str(s) for s in test[1:]]
            testname = "_".join(str(p) for p in testparams)

          full_command = testfile.replace("$D", os.path.dirname(run_tests.__file__)). \
                                  replace("$B", libtbx.env.under_build(test_module))
          shortpath = testfile.replace("$D/", "").replace("$B/", "build/")
          pytest_file_object = pytest.File(shortpath, self.session)
          yield LibtbxTest(testname, pytest_file_object, full_command, testparams)

    if path.basename == 'run_tests.py':
      return LibtbxRunTestsFile(path, parent)

  return pytest_collect_file

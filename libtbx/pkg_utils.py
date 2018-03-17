# Python package related functions
#
# To let libtbx modules specify and satisfy their requirements and capabilities.

from __future__ import division, absolute_import, print_function

import contextlib
import os
import sys

import libtbx.load_env

try:
  import pip
  import pkg_resources

  # Don't run if pip version < 9.0.0.
  # There is no technical reason for this, the code should work still.
  # But in the interest of not upsetting legacy build systems let's be cautious.
  if not all(symbol in dir(pkg_resources) for symbol in
             ('parse_version', 'require', 'DistributionNotFound', 'VersionConflict')) \
     or pkg_resources.parse_version(pip.__version__) < pkg_resources.parse_version('9.0.0'):
    pip = None
    pkg_resources = None
except ImportError:
  pip = None
  pkg_resources = None

try:
  import setuptools
except ImportError:
  setuptools = None

def require(pkgname, version=None):
  if not pip:
    print ("\n" + "=" * 80 + "\n\n"
         + "  WARNING: Can not verify python package requirements - pip/setuptools out of date\n"
         + "  Please update pip and setuptools by running:\n\n"
         + "    libtbx.python -m pip install pip setuptools --upgrade\n\n"
         + "  or following the instructions at https://pip.pypa.io/en/stable/installing/\n\n"
         + "=" * 80 + "\n")
    return False

  if not version:
    version = ''
  requirestring = pkgname + version
  try:
    print("requires %s, has %s" % (requirestring, pkg_resources.require(requirestring)[0].version))
    return True

  except pkg_resources.DistributionNotFound:
    currentversion = '(not determined)'
    project_name = pkgname
    action = 'install'
    print("requirement %s is not currently met, package not installed" % (requirestring))

  except pkg_resources.VersionConflict:
    currentversion = pkg_resources.require(pkgname)[0].version
    project_name = pkg_resources.require(pkgname)[0].project_name
    action = 'update'
    print("requirement %s is not currently met, current version %s" % (requirestring, currentversion))

  # Check if package can be updated
  for path_item in sys.path:
    egg_link = os.path.join(path_item, project_name + '.egg-link')
    if os.path.isfile(egg_link):
      with open(egg_link, 'r') as fh:
        print (("=" * 80 + "\n\n"
             + "     WARNING: Can not update package {package} automatically.\n"
             + "It is installed as editable package for development purposes. The currently\n"
             + "installed version, {currentversion}, is too old. The required version is {requirement}.\n"
             + "Please update the package manually in its installed location:\n\n"
             + "     {packagelocation}\n\n"
             + "=" * 80 + "\n\n").format(package=pkgname, currentversion=currentversion, requirement=version, packagelocation=fh.readline().strip()))
      return False

  if not os.path.isdir(libtbx.env.under_base('.')):
    print (("=" * 80 + "\n\n"
         + "     WARNING: Can not install/update package {package} automatically.\n"
         + "You are running in a base-less installation, which disables automatic package installation\n"
         + "by convention. cf. https://github.com/cctbx/cctbx_project/issues/151\n\n"
         + "Please install/update the package manually.\n\n"
         + "=" * 80 + "\n\n").format(package=pkgname, currentversion=currentversion, requirement=version))
    return False

  print("attempting {action} of {package}...".format(action=action, package=pkgname))
  exit_code = pip.main(['install', requirestring])
  if exit_code == 0:
    print("{action} successful".format(action=action))
    return True
  else:
    print("{action} failed. please check manually".format(action=action))
    return False

@contextlib.contextmanager
def _silence():
  '''Helper context which shuts up stdout and stderr.'''
  try:
    oldstdout = os.dup(sys.stdout.fileno())
    oldstderr = os.dup(sys.stderr.fileno())
    dest_file = open(os.devnull, 'w')
    os.dup2(dest_file.fileno(), sys.stdout.fileno())
    os.dup2(dest_file.fileno(), sys.stderr.fileno())
    yield
  finally:
    if oldstdout is not None:
      os.dup2(oldstdout, sys.stdout.fileno())
    if oldstderr is not None:
      os.dup2(oldstderr, sys.stderr.fileno())
    if dest_file is not None:
      dest_file.close()

def define_entry_points(epdict, **kwargs):
  '''A function to allow non-setuptools packages (ie. libtbx modules) to use
     the setuptools entry_points mechanism. Call this function from
     libtbx_refresh.py and pass a dictionary of entry points.'''
  # Determine the name of the calling module, and thus the internal module name
  # of the run_tests file. Use exception trick to pick up the current frame.
  try:
    raise Exception()
  except Exception:
    frame = sys.exc_info()[2].tb_frame.f_back
  # now in normal python land we could just do
  # caller = frame.f_globals['__name__']
  # alas we are in the libtbx shadow world and there is no __name__. x_x
  caller = os.path.abspath(frame.f_code.co_filename)  # Get the full path of the libtbx_refresh.py file.
  refresh_file, _ = os.path.splitext(caller)
  if not refresh_file.endswith('libtbx_refresh'):
    raise RuntimeError('Entry points can only be defined from within libtbx_refresh.py')
  # the name of the parent directory of libtbx_refresh.py is the caller name
  caller = os.path.basename(os.path.dirname(refresh_file))

  # No setuptools mechanism without setuptools.
  if not setuptools:
    raise ImportError("You must install setuptools to configure package {}. Run\n  libtbx.python -m pip install setuptools".format(caller))

  # Temporarily change to build/ directory. This is where a directory named
  # libtbx.{caller}.egg-info will be created containing the entry point info.
  try:
    curdir = os.getcwd()
    os.chdir(abs(libtbx.env.build_path))
    # Now trick setuptools into thinking it is in control here.
    try:
      argv_orig = sys.argv
      sys.argv = ['setup.py', 'develop']
      # And make it run quietly
      with _silence():
        setuptools.setup(
          name='libtbx.{}'.format(caller),
          description='libtbx entry point manager for {}'.format(caller),
          entry_points=epdict,
          **kwargs
        )
    finally:
      sys.argv = argv_orig
  finally:
    os.chdir(curdir)

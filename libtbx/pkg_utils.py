# Python package related functions
#
# To let libtbx modules specify and satisfy their requirements and capabilities.

from __future__ import absolute_import, division, print_function

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
  if pip:
    if pkg_resources.parse_version(pip.__version__) >= pkg_resources.parse_version('10.0.0'):
      import pip._internal
      pip_main = pip._internal.main
    else:
      pip_main = pip.main
except ImportError:
  pip = None
  pkg_resources = None

try:
  import setuptools
except ImportError:
  setuptools = None

def _notice(*lines, **context):
  print(os.linesep + "=" * 80 + os.linesep + os.linesep +
        os.linesep.join(l.format(**context) for l in lines) + os.linesep +
        os.linesep + "=" * 80 + os.linesep)

_defined_entrypoints = set()

def require(pkgname, version=None):
  '''Ensure a package requirement is met. Install or update package as required
     and print a warning message if this can't be done due to the local
     environment, or when automatic package management is disabled by setting
     the environment variable 'LIBTBX_DISABLE_UPDATES'.
     :param pkgname: A string describing the package requirement. This will
                     generally just be a package name, but package features
                     can be specified in square brackets. Features are not
                     enforced, but will be requested during installation and
                     update.
     :param version: An optional string describing version constraints. This
                     can be a minimum version, eg. '>=1.0', a maximum version,
                     eg. '<2', or both, eg. '>=4.5,<4.6'.
     :return: True when the requirement is met, False otherwise.'''

  if not pip:
    _notice("  WARNING: Can not verify python package requirements - pip/setuptools out of date",
            "  Please update pip and setuptools by running:", "",
            "    libtbx.python -m pip install pip setuptools --upgrade", "",
            "  or following the instructions at https://pip.pypa.io/en/stable/installing/")
    return False

  if not version:
    version = ''

  # package name without feature specification
  basepkgname = pkgname.split('[')[0]

  requirestring = pkgname + version
  baserequirestring = basepkgname + version
  try:
    try:
      print("requires %s, has %s" % (requirestring, pkg_resources.require(requirestring)[0].version))
      return True
    except pkg_resources.UnknownExtra:
      print("requires %s, has %s, but without features" % (requirestring, pkg_resources.require(baserequirestring)[0].version))
      return True

  except pkg_resources.DistributionNotFound:
    currentversion = '(not determined)'
    project_name = pkgname
    action = 'install'
    print("requirement %s is not currently met, package not installed" % (requirestring))

  except pkg_resources.VersionConflict:
    currentversion = pkg_resources.require(basepkgname)[0].version
    project_name = pkg_resources.require(basepkgname)[0].project_name
    action = 'update'
    print("requirement %s is not currently met, current version %s" % (requirestring, currentversion))

  # Check if package can be updated
  for path_item in sys.path:
    egg_link = os.path.join(path_item, project_name + '.egg-link')
    if os.path.isfile(egg_link):
      with open(egg_link, 'r') as fh:
        location = fh.readline().strip()
      _notice("    WARNING: Can not update package {package} automatically.", "",
              "It is installed as editable package for development purposes. The currently",
              "installed version, {currentversion}, is too old. The required version is {requirement}.",
              "Please update the package manually in its installed location:", "",
              "    {location}",
              package=pkgname, currentversion=currentversion, requirement=version, location=location)
      return False

  if not os.path.isdir(libtbx.env.under_base('.')):
    _notice("    WARNING: Can not {action} package {package} automatically.", "",
            "You are running in a base-less installation, which disables automatic package",
            "management by convention, see https://github.com/cctbx/cctbx_project/issues/151", "",
            "Please {action} the package manually.",
            package=pkgname, currentversion=currentversion, requirement=version, action=action)
    return False

  if os.getenv('LIBTBX_DISABLE_UPDATES') and os.getenv('LIBTBX_DISABLE_UPDATES').strip() not in ('0', ''):
    _notice("    WARNING: Can not {action} package {package} automatically.", "",
            "Environment variable LIBTBX_DISABLE_UPDATES is set. Please {action} manually.",
            package=pkgname, currentversion=currentversion, requirement=version, action=action)
    return False

  print("attempting {action} of {package}...".format(action=action, package=pkgname))
  exit_code = pip_main(['install', requirestring])
  if exit_code == 0:
    print("{action} successful".format(action=action))
    return True
  else:
    print("{action} failed. please check manually".format(action=action))
    return False

@contextlib.contextmanager
def _silence():
  '''Helper context which shuts up stdout.'''
  sys.stdout.flush()
  try:
    oldstdout = os.dup(sys.stdout.fileno())
    dest_file = open(os.devnull, 'w')
    os.dup2(dest_file.fileno(), sys.stdout.fileno())
    yield
  finally:
    if oldstdout is not None:
      os.dup2(oldstdout, sys.stdout.fileno())
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
  # Extract the caller name
  caller = frame.f_globals['__name__']
  if caller == '__main__':
    # well that is not very informative, is it.
    caller = os.path.abspath(frame.f_code.co_filename)  # Get the full path of the libtbx_refresh.py file.
    refresh_file, _ = os.path.splitext(caller)
    if not refresh_file.endswith('libtbx_refresh'):
      raise RuntimeError('Entry points can only be defined from within libtbx_refresh.py')
    # the name of the parent directory of libtbx_refresh.py is the caller name
    caller = os.path.basename(os.path.dirname(refresh_file))
  else:
    if not caller.endswith('.libtbx_refresh'):
      raise RuntimeError('Entry points can only be defined from within libtbx_refresh.py')
    caller = caller[:-15]

  # No setuptools mechanism without setuptools.
  if not setuptools:
    raise ImportError("You must install setuptools to configure package {}. Run\n  libtbx.pip install setuptools".format(caller))

  if caller in _defined_entrypoints:
    raise RuntimeError("Entry points have already been defined for package %s. There must only be a single call to "
                       "define_entry_points() per calling package" % caller)
  _defined_entrypoints.add(caller)

  print("Updating entry points for {caller}".format(caller=caller))
  for ep in epdict:
    print("  {n} entries for entry point {ep}".format(ep=ep, n=len(epdict[ep])))

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

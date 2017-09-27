# Python package related functions
#
# To let libtbx modules specify and satisfy their requirements.

from __future__ import division, absolute_import

import os
import sys
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
    print "requires %s, has %s" % (requirestring, pkg_resources.require(requirestring)[0].version)
    return True

  except pkg_resources.DistributionNotFound:
    currentversion = '(not determined)'
    project_name = pkgname
    action = 'install'
    print "requirement %s is not currently met, package not installed" % (requirestring)

  except pkg_resources.VersionConflict:
    currentversion = pkg_resources.require(pkgname)[0].version
    project_name = pkg_resources.require(pkgname)[0].project_name
    action = 'update'
    print "requirement %s is not currently met, current version %s" % (requirestring, currentversion)

  # Check if package can be updated
  for path_item in sys.path:
    egg_link = os.path.join(path_item, project_name + '.egg-link')
    if os.path.isfile(egg_link):
      with open(egg_link, 'r') as fh:
        print ("=" * 80 + "\n\n"
             + "     WARNING: Can not update package {package} automatically.\n"
             + "It is installed as editable package for development purposes. The currently\n"
             + "installed version, {currentversion}, is too old. The required version is {requirement}.\n"
             + "Please update the package manually in its installed location:\n\n"
             + "     {packagelocation}\n\n"
             + "=" * 80 + "\n\n").format(package=pkgname, currentversion=currentversion, requirement=version, packagelocation=fh.readline().strip())
        return False

  print "attempting {action} of {package}...".format(action=action, package=pkgname)
  exit_code = pip.main(['install', requirestring])
  if exit_code == 0:
    print "{action} successful".format(action=action)
    return True
  else:
    print "{action} failed. please check manually".format(action=action)
    return False

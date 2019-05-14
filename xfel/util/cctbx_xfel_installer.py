
# XXX This is not a standalone installer!  It must be used as part of the
# framework in libtbx/auto_build.

"""
Installer script for cctbx.xfel packages based on automatically generated
template.  This must be moved to the proper location to work.
"""

from __future__ import absolute_import, division, print_function
import os.path
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution

class installer (install_distribution.installer) :
  # XXX most settings can be edited here
  product_name = "cctbx.xfel"
  dest_dir_prefix = "xfel"
  make_apps = []
  configure_modules = install_distribution.installer.configure_modules + \
    ['dxtbx', 'wxtbx', "gltbx", "crys3d", "xfel","dials"]
  include_gui_packages = True
  base_package_options = ['--all']
  source_packages = [ "cctbx_bundle" ] + ['cbflib','labelit','dials']
  #

  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if (__name__ == "__main__") :
  installer(sys.argv[1:])

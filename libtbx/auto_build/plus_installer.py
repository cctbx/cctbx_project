
# XXX This is not a standalone installer!  It must be used as part of the
# framework in libtbx/auto_build.

"""
Installer script for CCTBX-plus packages based on automatically generated
template.  This must be moved to the proper location to work.
"""

from __future__ import absolute_import, division, print_function

import os.path
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if (not libtbx_path in sys.path):
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution

from libtbx.auto_build.bootstrap import CCTBXBuilder

class installer(install_distribution.installer):
  # XXX most settings can be edited here
  product_name = "CCTBX"
  dest_dir_prefix = "cctbx"
  make_apps = ["cctbx.image_viewer", "phenix.data_viewer"]
  configure_modules = CCTBXBuilder.LIBTBX
  include_gui_packages = True
  base_package_options = ['--all', "--sphinx", "--ipython"]
  modules = CCTBXBuilder.HOT + CCTBXBuilder.CODEBASES

  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if (__name__ == "__main__"):
  installer(sys.argv[1:]).install()

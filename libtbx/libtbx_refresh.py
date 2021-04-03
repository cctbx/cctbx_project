from __future__ import absolute_import, division, print_function

import os
import shutil

import libtbx.load_env
from libtbx.version import create_version_files

# try creating version files again
# this should already be done by bootstrap.py, but try anyway
# this will fail if git is not available
filenames = []
try:
  filenames = create_version_files()
except Exception:
  pass

# copy version file and header to build directory
try:
  shutil.copy(filenames[0], abs(libtbx.env.build_path))
  cctbx_include = libtbx.env.under_build('include/cctbx')
  if not os.path.exists(cctbx_include):
    os.mkdir(cctbx_include)
  shutil.copy(filenames[1], cctbx_include)
except Exception:
  pass

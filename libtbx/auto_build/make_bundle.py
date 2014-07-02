
"""
Create a "bundle" (.tar.gz) of all Python modules and compiled code in a
product.  The target directory is expected to look something like this:

CCTBX-<version>/
CCTBX-<version>/build/
CCTBX-<version>/build/<mtype>/
CCTBX-<version>/build/<mtype>/base/
CCTBX-<version>/build/<mtype>/lib/
CCTBX-<version>/cctbx_project/

plus any number of module directories in the top level.  The resulting bundle
will be named bundle-<version>-<mtype>.tar.gz.
"""

from __future__ import division
from optparse import OptionParser
import os.path as op
import shutil
import time
import os
import sys
# local imports
# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from libtbx.auto_build.installer_utils import *
from libtbx.auto_build import rpath


def run (args, out=sys.stdout) :
  datestamp = time.strftime("%Y_%m_%d", time.localtime())
  parser = OptionParser()
  parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
    help="Temporary staging directory", default=os.getcwd())
  parser.add_option("--version", dest="version", action="store",
    help="Version number or code", default=datestamp)
  parser.add_option("--mtype", dest="mtype", action="store",
    help="Architecture type", default=machine_type())
  parser.add_option("--ignore", dest="ignore", action="store",
    help="Subdirectories to ignore", default="")
  parser.add_option("--remove_src", dest="remove_src", action="store",
    help="Remove compiled source files (.h, .cpp, etc.)", default=False)
  options, args = parser.parse_args(args)
  target_dir = args[0]
  assert op.isdir(target_dir), target_dir
  os.chdir(options.tmp_dir)
  build_dir = op.join(target_dir, "build", options.mtype)
  tmp_dir = op.join(options.tmp_dir, "tmp_%s" % options.version)
  assert op.isdir(build_dir), build_dir
  if os.exists(tmp_dir) :
    shutil.rmtree(tmp_dir)
  os.mkdir(tmp_dir)
  os.chdir(tmp_dir)
  # copy over non-compiled files
  print >> out, "Copying base modules..."
  ignore_dirs = options.ignore.split(",")
  for file_name in os.listdir(target_dir) :
    if (file_name == "build") or (file_name in ignore_dirs) :
      continue
    full_path = op.join(target_dir, file_name)
    if op.isdir(full_path) :
      print >> out, "  copying %s..." % file_name
      copy_tree(full_path, op.join(tmp_dir, file_name))
  # build directory
  tmp_build_dir = op.join(tmp_dir, "build", options.mtype)
  os.makedirs(tmp_build_dir)
  for dir_name in ["lib", "base"] :
    full_path = op.join(build_dir, dir_name)
    assert op.isdir(full_path)
    copy_tree(full_path, op.join(tmp_build_dir, dir_name))
  # copy over build executable directories
  for file_name in os.listdir(build_dir) :
    full_path = op.join(build_dir)
    if op.isdir(full_path) :
      module_name = file_name
      for file_name in os.listdir(full_path) :
        if (file_name == "exe") :
          copy_tree(op.join(full_path, file_name),
                    op.join(tmp_build_dir, module_name, file_name))
  # remove unnecessary base directories/files
  for dir_name in [
      "base/bin/gtk-demo",
      "base/man",
      "base/doc",
      "base/info",
      "base/share/gtk-doc",
      "base/share/locale",
    ] :
    full_path = op.join(tmp_build_dir, dir_name)
    shutil.rmtree(full_path)
  # XXX what about base/include?
  # delete unnecessary files
  find_and_delete_files(tmp_dir, file_ext=".pyc")
  find_and_delete_files(tmp_dir, file_ext=".pyo")
  find_and_delete_files(tmp_dir, file_name=".sconsign")
  find_and_delete_files(tmp_dir, file_name="CVS")
  find_and_delete_files(tmp_dir, file_name=".svn")

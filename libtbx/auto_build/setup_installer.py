#!/usr/bin/python

"""
Script to set up an installer directory tree and copy over most of the
necessary files.  We used to just keep the entire (Phenix) installer in a
separate SVN tree, but this is inconvenient when we have multiple packages
using the same system and also many third-party dependencies which need to be
kept up to date.  Note that this script provides only the bare minimum
functionality for building CCTBX installers, and other distributions will
probably need to perform additional modifications to the installer tree
before it can be tarred.
"""

from __future__ import division
from optparse import OptionParser
import os.path as op
import shutil
import time
import stat
import os
import sys
# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from package_defs import *

def run (args) :
  parser = OptionParser()
  parser.add_option("--version", dest="version", action="store",
    help="Package version", default=time.strftime("%Y_%m_%d",time.localtime()))
  parser.add_option("--pkg_dir", dest="pkg_dirs", action="append",
    help="Directory with source packages", default=None)
  parser.add_option("--basic", dest="basic", action="store_true",
    help="Only include basic prerequisite packages", default=False)
  parser.add_option("--xia2", dest="xia2", action="store_true",
    help="Include xia2 dependencies", default=False)
  parser.add_option("--dials", dest="dials", action="store_true",
    help="Include DIALS dependencies", default=False)
  parser.add_option("--gui", dest="gui", action="store_true",
    help="Include GUI dependencies", default=False)
  parser.add_option("-a", "--all", dest="all", action="store_true",
    help="Include all recommended dependencies", default=False)
  parser.add_option("--dest", dest="dest", action="store",
    help="Destination folder", default=os.getcwd())
  parser.add_option("--readme", dest="readme", action="store",
    help="Readme file", default=op.join(libtbx_path, "COPYRIGHT_2_0.txt"))
  parser.add_option("--license", dest="license", action="store",
    help="License file", default=op.join(libtbx_path, "LICENSE_2_0.txt"))
  parser.add_option("--script", dest="script", action="store",
    help="Final installation script", default=None)
  parser.add_option("--modules", dest="modules", action="store",
    help="Local modules to include", default=None)
  parser.add_option("--product_name", dest="product_name", action="store",
    help="Name of installed package", default=None)
  options, args = parser.parse_args(args)
  assert len(args) == 1
  package_name = args[-1]
  if (options.product_name is None) :
    options.product_name = package_name
  if (options.script is not None) :
    assert op.isfile(options.script)
    options.script = op.abspath(options.script)
  os.chdir(options.dest)
  # setup directory structure
  installer_dir = "%s-installer-%s" % (package_name, options.version)
  os.mkdir(installer_dir)
  os.chdir(installer_dir)
  os.mkdir("bin")
  os.mkdir("lib")
  os.mkdir("source")
  os.mkdir("dependencies")
  # copy over libtbx
  lib_dir = op.join(os.getcwd(), "lib")
  shutil.copytree(libtbx_path, op.join(lib_dir, "libtbx"))
  find_and_delete_files(lib_dir, file_ext=".pyc")
  # write VERSION
  open("VERSION", "w").write(options.version)
  if op.isfile(options.readme) :
    open("README", "w").write(open(options.readme).read())
  if op.isfile(options.license) :
    open("LICENSE", "w").write(open(options.license).read())
  fetch_all_dependencies(
    dest_dir=op.join(options.dest, installer_dir, "dependencies"),
    log=sys.stdout,
    pkg_dirs=options.pkg_dirs,
    gui_packages=(options.dials or options.gui or options.all),
    dials_packages=(options.dials or options.xia2 or options.all))
  # local packages
  fetch_package = fetch_packages(
    dest_dir=op.join(options.dest, installer_dir, "source"),
    log=sys.stdout,
    pkg_dirs=options.pkg_dirs,
    copy_files=True)
  fetch_package(
    pkg_name="cctbx_bundle_for_installer.tar.gz",
    pkg_url="http://cci.lbl.gov/build/results/current")
  if (options.modules is not None) :
    for module_name in options.modules :
      pass # TODO
  os.chdir(op.join(options.dest, installer_dir))
  # actual Python installer script
  if (options.script is not None) :
    assert op.isfile(options.script)
    open("bin/install.py", "w").write(open(options.script).read())
  else :
    # default stub.  this is pretty minimal but it will work for simple
    # packages.
    modules_list = []
    if (options.modules is not None) :
      modules_list = options.modules.split(",")
    base_package_options = []
    if (options.gui or options.all) :
      base_package_options.append("--gui")
    if (options.dials or options.all or options.xia2) :
      base_package_options.append("--dials")
    f = open("bin/install.py", "w")
    f.write("""\
import os.path
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution

class installer (install_distribution.installer) :
  product_name = "%(pkg)s"
  dest_dir_prefix = "%(product)s"
  make_apps = []
  configure_modules = install_distribution.installer.configure_modules + \\
    %(modules)s
  include_gui_packages = %(gui)s
  base_package_options = %(baseopts)s
  source_packages = [ "cctbx_bundle" ] + %(modules)s

  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

if (__name__ == "__main__") :
  installer(sys.argv[1:])
""" % { "pkg" : package_name, "product" : options.product_name,
        "modules" : modules_list, "baseopts" : base_package_options,
        "gui" : (options.gui or options.all) })
    f.close()
  # write executable Bash script wrapping Python script
  f = open("install", "w")
  f.write("""\
#!/bin/bash
#
#

UNAME=`uname`
if [ -z "$PYTHON_EXE" ]; then
  PYTHON_EXE='/usr/bin/python'
  if [ -f "/usr/bin/python2.7" ]; then
    PYTHON_EXE='/usr/bin/python2.7'
  elif [ -f "/usr/bin/python2.6" ]; then
    PYTHON_EXE='/usr/bin/python2.6'
  elif [ -f "/usr/bin/python2.5" ]; then
    PYTHON_EXE='/usr/bin/python2.5'
  elif [ -f "/usr/bin/python2" ]; then
    PYTHON_EXE='/usr/bin/python2'
  fi
fi
$PYTHON_EXE ./bin/install.py $@
""")
  f.close()
  st = os.stat("install")
  os.chmod("install", st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

if (__name__ == "__main__") :
  run(sys.argv[1:])

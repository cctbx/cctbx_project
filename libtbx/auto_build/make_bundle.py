#!/usr/bin/python

"""
Create a "bundle" (.tar.gz) of all Python modules and compiled code in a
product.  The target directory is expected to look something like this:

CCTBX-<version>/
CCTBX-<version>/base/
CCTBX-<version>/build/
CCTBX-<version>/build/lib/
CCTBX-<version>/cctbx_project/

plus any number of module directories in the top level.  The resulting bundle
will be named bundle-<version>-<mtype>.tar.gz.

Since the base modules take an especially long time to compile, they are now
part of a separate bundle.  This will allow re-use of precompiled base packages
with the latest source, which should speed up installer generation.
"""

from __future__ import absolute_import, division, print_function
from optparse import OptionParser
import os.path as op
import shutil
import time
import os
import sys
# local imports
# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(op.dirname(__file__))))
if (not libtbx_path in sys.path):
  print(libtbx_path)
  sys.path.append(libtbx_path)
from libtbx.auto_build.installer_utils import *
from libtbx.auto_build import rpath


def run(args, out=sys.stdout):
  datestamp = time.strftime("%Y_%m_%d", time.localtime())
  parser = OptionParser()
  parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
    help="Temporary staging directory", default=os.getcwd())
  parser.add_option("--version", dest="version", action="store",
    help="Version number or code", default=datestamp)
  parser.add_option("--mtype", dest="mtype", action="store",
    help="Architecture type", default=machine_type())
  parser.add_option("--ignore", dest="ignore", action="append",
    help="Subdirectories to ignore", default=[])
  parser.add_option("--remove_src", dest="remove_src", action="store_true",
    help="Remove compiled source files (.h, .cpp, etc.)", default=False)
  parser.add_option("--dest", dest="dest", action="store",
    help="Destination directory for bundle tarfiles", default=None)
  parser.add_option("--verbose", dest="verbose", action="store_true")
  options, args = parser.parse_args(args)
  target_dir = args[0]
  assert op.isdir(target_dir), target_dir
  target_dir = op.abspath(target_dir)
  if (options.dest is not None):
    assert op.isdir(options.dest)
  os.chdir(options.tmp_dir)
  pkg_dir = op.basename(target_dir)
  build_dir = op.join(target_dir, "build")
  base_dir = op.join(target_dir, "base")

  # XXX a bit of a hack - if 'modules' subdirectory exists, it is assumed that
  # all source/data packages residue there, otherwise they must be at the top
  # level
  modules_dir = op.join(target_dir, "modules")
  if (not op.isdir(modules_dir)):
    modules_dir = target_dir
  stdout_old = sys.stdout
  if (not options.verbose):
    f = open("rpath.log", "w")
    sys.stdout = f
  print("Setting rpath in base packages...", file=out)
  rpath.run([base_dir])
  print("Setting rpath in shared libraries...", file=out)
  rpath.run([build_dir])
  sys.stdout = stdout_old

  # create temp dir
  tmp_dir = op.join(options.tmp_dir, "%s_tmp" % pkg_dir)
  assert op.isdir(build_dir), build_dir
  if op.exists(tmp_dir):
    shutil.rmtree(tmp_dir)
  os.mkdir(tmp_dir)
  os.chdir(tmp_dir)

  # base and build/lib directories
  print("Copying dependencies...", file=out)
  copy_tree(op.join(target_dir, "base"), op.join(tmp_dir, "base"))

  print("Copying shared libraries...", file=out)
  tmp_build_dir = op.join(tmp_dir, "build")
  os.makedirs(tmp_build_dir)

  # save mtype information (for hypothetical future update mechanism)
  open(op.join(tmp_build_dir, "MTYPE"), "w").write(options.mtype)
  copy_tree(op.join(build_dir, "lib"), op.join(tmp_build_dir, "lib"))

  # copy over non-compiled files
  print("Copying base modules...", file=out)
  for file_name in os.listdir(modules_dir):
    if (file_name in ["build", "base"]) or (file_name in options.ignore):
      continue
    full_path = op.join(modules_dir, file_name)
    if op.isdir(full_path):
      print("  copying %s..." % file_name, file=out)
      copy_tree(full_path, op.join(tmp_dir, file_name))
      call("chmod -R a+rX %s" % op.join(tmp_dir, file_name))

  # remove unnecessary base directories/files
  for dir_name in [
      "base/bin/gtk-demo",
      "base/man",
      "base/doc",
      "base/info",
      "base/share/gtk-doc",
      "base/share/locale",
      "base/lib/python2.7/test",
    ] :
    full_path = op.join(tmp_dir, dir_name)
    if op.exists(full_path):
      shutil.rmtree(full_path)

  site_pkg_dir = op.join(tmp_dir, "base/lib/python2.7/site-packages")
  find_and_delete_files(tmp_dir, file_name="tests")
  if sys.platform.startswith("linux"):
    strip_libs(op.join(tmp_dir, "base", "lib"), log=out)

  # XXX what about base/include?

  # copy over build executable directories
  print("Copying standalone executables...", file=out)
  # for j in [i for i in os.listdir(build_dir) if os.path.isdir(os.path.join(build_dir,i,"exe"))]:
  for j in [i for i in os.listdir(build_dir) if os.path.isdir(os.path.join(build_dir, i, "exe"))]:
    print("->", op.join(build_dir, j, "exe"), file=out)
    copy_tree(op.join(build_dir, j, "exe"), op.join(tmp_build_dir, j, "exe"))

  # delete unnecessary files
  print("Deleting unnecessary files.", file=out)
  find_and_delete_files(tmp_dir, file_ext=".pyc")
  find_and_delete_files(tmp_dir, file_ext=".o")
  find_and_delete_files(tmp_dir, file_ext=".pyo")
  find_and_delete_files(tmp_dir, file_name=".sconsign")
  find_and_delete_files(tmp_dir, file_name="CVS")
  find_and_delete_files(tmp_dir, file_name=".svn")
  if (options.remove_src):
    find_and_delete_files(tmp_dir, file_ext=".cpp")
    find_and_delete_files(tmp_dir, file_ext=".hpp")
    find_and_delete_files(tmp_dir, file_ext=".cc")
    find_and_delete_files(tmp_dir, file_ext=".c")
    find_and_delete_files(tmp_dir, file_ext=".h")

  # TODO strip objects?
  os.chdir(tmp_dir)
  call("chmod -R a+rX %s" % op.join(tmp_dir, "base"))
  call("chmod -R a+rX %s" % op.join(tmp_dir, "build"))
  # create base bundle
  base_tarfile = "../base-%(version)s-%(mtype)s.tar.gz" % \
    {"version":options.version, "mtype":options.mtype}
  call("tar -czf %(tarfile)s base" %
    {"tarfile":base_tarfile}, log=out)
  shutil.rmtree("base")

  assert op.isfile(base_tarfile)
  if (options.dest is not None):
    shutil.move(base_tarfile, options.dest)
    base_tarfile = op.join(options.dest, op.basename(base_tarfile))
  print("  created base bundle %s" % base_tarfile, file=out)

  # create the product bundle
  build_tarfile = "../build-%(version)s-%(mtype)s.tar.gz" % \
    {"version":options.version, "mtype":options.mtype}

  modules_tarfile = "../modules-%(version)s-%(mtype)s.tar.gz" % \
    {"version":options.version, "mtype":options.mtype}

  call("tar -czf %(tarfile)s build" % {"tarfile":build_tarfile}, log=out)
  shutil.rmtree("build")
  call("tar -czf %(tarfile)s ." % {"tarfile":modules_tarfile}, log=out)

  assert op.isfile(build_tarfile)
  assert op.isfile(modules_tarfile)
  if (options.dest is not None):
    shutil.move(build_tarfile, options.dest)
    build_tarfile = op.join(options.dest, op.basename(build_tarfile))
    shutil.move(modules_tarfile, options.dest)
    modules_tarfile = op.join(options.dest, op.basename(modules_tarfile))
  print("  created bundle %s" % build_tarfile, file=out)
  print("  created bundle %s" % modules_tarfile, file=out)
  shutil.rmtree(tmp_dir)

if (__name__ == "__main__"):
  run(sys.argv[1:])

#!/usr/bin/python
#
# MASTER SCRIPT FOR BINARY UNIX INSTALLER ASSEMBLY
#

from __future__ import division
from optparse import OptionParser
import os.path as op
import time
import os
import sys
libtbx_path = op.abspath(op.dirname(op.dirname(op.dirname(__file__))))
if (not libtbx_path in sys.path) :
  print libtbx_path
  sys.path.append(libtbx_path)
from libtbx.auto_build import setup_installer
from libtbx.auto_build import make_bundle
from libtbx.auto_build.installer_utils import *

class installer_builder (object) :
  #---------------------------------------------------------------------
  #
  # Label for the overall package
  product_name = "CCTBX"
  # prefix for installer directory and tarball
  pkg_prefix = "cctbx"
  # script to drive actual installer
  installer_script = "cctbx_project/libtbx/auto_build/plus_installer.py"
  # optional license file
  license = "cctbx_project/libtbx/LICENSE_2_0.txt"
  # optional directory containing additional binary packages
  bin_dir = None
  # optional text files to be placed in top-level installer directory
  readme_files = [
    "cctbx_project/libtbx/COPYRIGHT_2_0.txt",
  ]
  # modules that we need to grab source for (in addition to the cctbx_bundle
  # package)
  source_modules = [
    "cbflib",
    "annlib",
    "annlib_adaptbx",
  ]
  # modules that are placed in the 'base' directory of the installer rather
  # than being bundled up with everything else.
  # XXX this is probably historical legacy rather than an actual necessity;
  # but nice to keep some flexibility in packaging
  base_modules = [
  ]
  # modules that may be part of the build, but definitely *not* desired for
  # the installer.  this can include regression directories, experimental
  # packages, or third-party code with an incompatible license.
  exclude_build_modules = [
    "phenix_regression",
    "phenix_dev",
    "chem_data",
  ]
  #
  # CODE BELOW THIS LINE SHOULD NOT NEED TO BE MODIFIED
  #---------------------------------------------------------------------
  def __init__ (self, args) :
    parser = OptionParser()
    parser.add_option("--version", dest="version", action="store",
      help="Package version",
      default=time.strftime("%Y_%m_%d", time.localtime()))
    parser.add_option("--pkg_dir", dest="pkg_dirs", action="append",
      help="Directory with source packages", default=None)
    parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
      help="Temporary staging directory", default=os.getcwd())
    parser.add_option("--mtype", dest="mtype", action="store",
      help="Architecture type", default=machine_type())
    parser.add_option("--host-tag", dest="host_tag", action="store",
      help="Host tag (OS/distribution label)", default=None)
    parser.add_option("--remove_src", dest="remove_src", action="store_true",
      help="Remove compiled source files (.h, .cpp, etc.)", default=False)
    options, args = parser.parse_args(args)
    assert len(args) == 1
    builder_dir = args[0]
    assert op.isdir(builder_dir), builder_dir
    base_dir = op.join(builder_dir, "base")
    build_dir = op.join(builder_dir, "build")
    assert op.isdir(base_dir), base_dir
    assert op.isdir(build_dir), build_dir
    modules_dir = op.join(builder_dir, "modules")
    if (not op.isdir(modules_dir)) :
      modules_dir = builder_dir
    def full_path (path_name) :
      if op.isabs(path_name) :
        return path_name
      else :
        path1 = op.join(modules_dir, path_name)
        if op.exists(path1) :
          return path1
        else :
          path2 = op.join(builder_dir, path_name)
          if op.exists(path2) :
            return path2
          else :
            raise RuntimeError("Can't find path %s" % path_name)
    os.chdir(options.tmp_dir)
    # setup basic installer directory
    setup_args = [
      "--version=%s" % options.version,
      "--binary",
      "--script=%s" % full_path(options.installer_script),
      "--product_name=%s" % self.product_name,
      "--pkg_dir=%s" % modules_dir,
    ]
    if (len(self.readme_files) > 0) :
      for readme in self.readme_files :
        setup_args.append("--readme=%s" % full_path(readme))
    if (len(self.base_modules) > 0) :
      setup_args.append("--base-modules=%s" % ",".join(self.base_modules))
    if (self.license) :
      setup_args.append("--license=%s" % full_path(self.license))
    setup_installer.run(args=setup_args + [ self.pkg_name ])
    installer_dir = self.pkg_name + "-installer-" + options.version
    assert op.isdir(installer_dir)
    # create bundles of base, build, and module directories
    os.chdir(installer_dir)
    bundle_args = [
      "--dest=%s" % os.getcwd(),
      "--version=%s" % options.version,
      "--verbose",
    ]
    if (len(self.exclude_build_modules) > 0) :
      for module in self.exclude_build_modules :
        bundle_args.append("--ignore=%s" % module)
    make_bundle.run(args=bundle_args + [ self.builder_dir ])
    # package the entire mess into the complete installer
    find_and_delete_files(os.getcwd(), file_ext=".pyc")
    os.chdir(options.tmp_dir)
    tar_prefix = installer_dir
    if (options.host_tag is not None) :
      tar_prefix += "-" + options.host_tag
    else :
      tar_prefix += "-" + options.mtype
    tar_name = tar_prefix + ".tar.gz"
    call("tar czf %s %s" % (tar_name, installer_dir))
    print "Wrote %s" % tar_name

if (__name__ == "__main__") :
  installer_builder(sys.argv[1:])

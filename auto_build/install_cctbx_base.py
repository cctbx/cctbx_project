
"""
Automated build of CCTBX and all dependencies for Linux and Mac platforms.  This
script will download and install the current Python distribution provided by LBL,
plus any packages required for various CCTBX-dependent apps to function.  In the
future this will be used as the core of the CCI nightly build system and the
Phenix installer.
"""

from __future__ import division
from cctbx_installer_utils import *
from package_defs import *
from optparse import OptionParser
import urllib2
import os
import sys

class installer (object) :
  # XXX various defaults go here, to be overridden in subclasses
  build_xia2_dependencies = False
  build_gui_dependencies = False
  build_labelit_dependencies = False

  def __init__ (self, args, log=sys.stdout) :
    self.log = log
    parser = OptionParser()
    parser.add_option("--build_dir", dest="build_dir", action="store",
      help="Build directory", default=os.path.join(os.getcwd(), "cctbx_build"))
    parser.add_option("--src_dir", dest="src_dir", action="store",
      help="Source directory", default=os.path.join(os.getcwd(), "cctbx_sources"))
    parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
      help="Temporary directory", default=os.path.join(os.getcwd(), "build_tmp"))
    parser.add_option("--nproc", dest="nproc", action="store", type="int",
      help="Number of processors", default=1)
    parser.add_option("--verbose", dest="verbose", action="store_true",
      help="Verbose output", default=False)
    parser.add_option("--pkg_dir", dest="pkg_dir", action="store",
      help="Directory with source packages", default=None)
    parser.add_option("--basic", dest="basic", action="store_true",
      help="Only install basic prerequisites", default=False)
    parser.add_option("--labelit", dest="labelit", action="store_true",
      help="Build LABELIT dependencies", default=self.build_labelit_dependencies)
    parser.add_option("--xia2", dest="xia2", action="store_true",
      help="Build xia2 dependencies", default=self.build_xia2_dependencies)
    parser.add_option("--gui", dest="gui", action="store_true",
      help="Build GUI dependencies", default=self.build_gui_dependencies)
    parser.add_option("--use_system_python", dest="use_system_python",
      action="store_true", default=False)
    parser.add_option("--all", dest="build_all", action="store_true",
      help="Build all recommended dependencies", default=False)
    options, args = parser.parse_args(args)
    self.tmp_dir = options.tmp_dir
    self.src_dir = options.src_dir
    self.build_dir = options.build_dir
    self.pkg_dir = options.pkg_dir
    self.base_dir = os.path.join(self.build_dir, "base")
    for dir_name in [self.tmp_dir, self.src_dir, self.build_dir, self.base_dir] :
      if (not os.path.isdir(dir_name)) :
        print >> log, "Creating %s" % dir_name
        os.makedirs(dir_name)
    self.nproc = options.nproc
    self.verbose = options.verbose
    self.options = options
    self.python_exe = None
    self.build_cctbx_dependencies()

  def build_cctbx_dependencies (self) :
    os.chdir(self.tmp_dir)
    options = self.options
    if (not options.use_system_python) :
      self.python_exe = self.build_python()
    else :
      pass # TODO ???
    assert os.path.exists(self.python_exe)
    if (not options.basic) :
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=NUMPY_PKG,
        pkg_name_simple="numpy",
        confirm_import_module="numpy")
    if (options.labelit) or (options.gui) or (options.build_all) :
      self.build_imaging()
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=REPORTLAB_PKG,
        pkg_name_simple="reportlab",
        confirm_import_module="reportlab")
    if (options.xia2) or (options.build_all) :
      self.build_hdf5()
    self.print_sep()
    print >> self.log, "Dependencies finished building."

  def call (self, args, log=None) :
    if (log is None) : log = self.log
    return call(args, log=log)

  def chdir (self, dir_name, log=None) :
    if (log is None) : log = self.log
    print >> log, "cd \"%s\"" % dir_name
    os.chdir(dir_name)

  def print_sep (self, char="#") :
    print >> self.log, ""
    print >> self.log, char*80
    print >> self.log, ""

  def start_building_package (self, pkg_name) :
    install_log = os.path.join(self.tmp_dir, pkg_name + "_install_log")
    self.print_sep()
    print >> self.log, "Installing %s..." % pkg_name
    print >> self.log, "  log file is %s" % install_log
    return open(install_log, "w")

  def patch_src (self, src_file, target, replace_with) :
    os.rename(src_file, src_file + ".dist")
    src_in = open(src_file + ".dist")
    src_out = open(src_file, "w")
    for line in src_in.readlines() :
      src_out.write(line.replace(target, replace_with))
    src_in.close()
    src_out.close()

  def fetch_package (self, pkg_url, pkg_name) :
    os.chdir(self.tmp_dir)
    print >> self.log, "  getting package %s..." % pkg_name
    if (self.pkg_dir is not None) :
      static_file = os.path.join(self.pkg_dir, pkg_name)
      if (os.path.exists(static_file)) :
        print >> self.log, "  %s" % static_file
        return static_file
    if (os.path.exists(pkg_name)) :
      print >> self.log, "  ./%s" % pkg_name
      return os.path.join(self.tmp_dir, pkg_name)
    else :
      full_url = "%s/%s" % (pkg_url, pkg_name)
      self.log.write("    downloading %s : " % full_url)
      f = open(pkg_name, "wb")
      data = urllib2.urlopen(full_url).read()
      self.log.write("%d bytes\n" % len(data))
      self.log.flush()
      f.write(data)
      f.close()
      return os.path.join(self.tmp_dir, pkg_name)

  def build_python (self) :
    self.chdir(self.tmp_dir)
    log = self.start_building_package("Python")
    python_tarball = self.fetch_package(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=PYTHON_PKG)
    python_dir = untar(python_tarball)
    self.chdir(python_dir, log=log)
    if (sys.platform == "darwin") :
      configure_args = [
        "--prefix=\"%s\"" % self.base_dir,
        "--enable-framework=\"%s\"" % self.base_dir,]
      self.call("./configure %s" % " ".join(configure_args), log=log)
      self.call("make -j %d" % self.nproc, log=log)
      targets = "bininstall libinstall libainstall inclinstall sharedinstall"
      self.call("make -j %d %s"% (self.nproc, targets), log=log)
      self.chdir("Mac", log=log)
      self.call("make install_Python install_pythonw install_versionedtools",
        log=log)
      self.chdir(python_dir, log=log)
      self.call("make frameworkinstallstructure frameworkinstallmaclib \
        frameworkinstallunixtools FRAMEWORKUNIXTOOLSPREFIX=\"%s\"" %
        self.base_dir, log=log)
    else :
      configure_args = ["--prefix=\"%s\"" % self.base_dir,]
      self.call("./configure %s" % " ".join(configure_args), log=log)
      self.call("make -j %d" % self.nproc, log=log)
      self.call("make install", log=log)
    log.close()
    return os.path.join(self.base_dir, "bin", "python")

  def build_python_module_simple (self,
      pkg_url,
      pkg_name,
      pkg_name_simple,
      callback_before_build=None,
      callback_after_build=None,
      confirm_import_module=None) :
    os.chdir(self.tmp_dir)
    pkg_log = self.start_building_package(pkg_name_simple)
    pkg = self.fetch_package(pkg_url, pkg_name)
    pkg_log = open("%s_install_log" % pkg_name_simple, "w")
    pkg_dir = untar(pkg, log=pkg_log)
    self.chdir(pkg_dir, log=pkg_log)
    if (callback_before_build is not None) :
      assert callback_before_build(pkg_log)
    self.call("%s setup.py install" % self.python_exe, log=pkg_log)
    if (callback_after_build is not None) :
      assert callback_after_build(pkg_log)
    os.chdir(self.tmp_dir)
    if (confirm_import_module is not None) :
      print >> self.log, "  trying to import installed module..."
      self.call("%s -c 'import %s'" % (self.python_exe, confirm_import_module))

  def build_imaging (self, patch_src=True) :
    def patch_imaging_src (out) :
      print >> out, "  patching libImaging/ZipEncode.c"
      self.patch_src(src_file="libImaging/ZipEncode.c",
                     target="Z_DEFAULT_COMPRESSION",
                     replace_with="Z_BEST_SPEED")
      return True
    callback = None
    if (patch_src) :
      callback = patch_imaging_src
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=IMAGING_PKG,
      pkg_name_simple="Imaging",
      callback_before_build=patch_imaging_src,
      confirm_import_module="Image")

  def build_hdf5 (self) :
    os.chdir(self.tmp_dir)
    pkg_log = self.start_building_package("HDF5")
    pkg = self.fetch_package(BASE_XIA_PKG_URL, HDF5_PKG)
    pkg_dir = untar(pkg, log=pkg_log)
    self.chdir(pkg_dir, log=pkg_log)
    print >> pkg_log, "Building base HDF5 library..."
    self.call("./configure --prefix=\"%s\"" % self.base_dir, log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)
    print >> pkg_log, "Building h5py..."
    pkg = self.fetch_package(BASE_XIA_PKG_URL, H5PY_PKG)
    pkg_dir = untar(pkg, log=pkg_log)
    self.chdir(pkg_dir, log=pkg_log)
    self.call("%s setup.py build --hdf5=\"%s\"" % (self.python_exe,
      self.base_dir), log=pkg_log)
    self.call("%s setup.py install" % (self.python_exe), log=pkg_log)
    self.call("%s setup.py test" % self.python_exe, log=pkg_log)

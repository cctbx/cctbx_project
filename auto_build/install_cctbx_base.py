
"""
Automated build of CCTBX dependencies for Linux and Mac platforms.
This script will download and install the current Python distribution provided
by LBL, plus any packages required for various CCTBX-dependent apps to
function.  In the future this will be used as the core of the CCI nightly build
system and the Phenix installer.
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
    print >> log, """
  ****************************************************************************
                 Automated CCTBX + dependencies build script
         This is experimental software - not guaranteed to work yet.
                report problems to cctbx-dev@cci.lbl.gov
  ****************************************************************************
"""
    parser = OptionParser()
    parser.add_option("--build_dir", dest="build_dir", action="store",
      help="Build directory", default=os.path.join(os.getcwd(), "cctbx_build"))
    parser.add_option("--src_dir", dest="src_dir", action="store",
      help="Source directory",
      default=os.path.join(os.getcwd(), "cctbx_sources"))
    parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
      help="Temporary directory",
      default=os.path.join(os.getcwd(), "build_tmp"))
    parser.add_option("--nproc", dest="nproc", action="store", type="int",
      help="Number of processors", default=1)
    parser.add_option("--verbose", dest="verbose", action="store_true",
      help="Verbose output", default=False)
    parser.add_option("--pkg_dir", dest="pkg_dir", action="store",
      help="Directory with source packages", default=None)
    parser.add_option("--basic", dest="basic", action="store_true",
      help="Only install basic prerequisites", default=False)
    parser.add_option("--labelit", dest="labelit", action="store_true",
      help="Build LABELIT dependencies",
      default=self.build_labelit_dependencies)
    parser.add_option("--xia2", dest="xia2", action="store_true",
      help="Build xia2 dependencies", default=self.build_xia2_dependencies)
    parser.add_option("--gui", dest="build_gui", action="store_true",
      help="Build GUI dependencies", default=self.build_gui_dependencies)
    #parser.add_option("--use_system_python", dest="use_system_python",
    #  action="store_true", default=False)
    parser.add_option("--all", dest="build_all", action="store_true",
      help="Build all recommended dependencies", default=False)
    parser.add_option("--scipy", dest="build_scipy", action="store_true",
      help="Build SciPy (requires Fortran compiler)", default=False)
    parser.add_option("--debug", dest="debug", action="store_true",
      help="Build in debugging mode", default=False)
    parser.add_option("--no-download", dest="no_download", action="store_true",
      help="Use only local packages (no downloads)", default=False)
    options, args = parser.parse_args(args)
    # basic setup
    self.tmp_dir = options.tmp_dir
    self.src_dir = options.src_dir
    self.build_dir = options.build_dir
    self.pkg_dir = options.pkg_dir
    self.base_dir = os.path.join(self.build_dir, "base")
    print >> log, "Setting up directories..."
    for dir_name in [self.tmp_dir,self.src_dir,self.build_dir,self.base_dir] :
      if (not os.path.isdir(dir_name)) :
        print >> log, "  creating %s" % dir_name
        os.makedirs(dir_name)
    self.nproc = options.nproc
    self.verbose = options.verbose
    self.options = options
    self.python_exe = None
    self.include_dirs = []
    self.lib_dirs = []
    self.flag_is_linux = sys.platform.startswith("linux")
    self.flag_is_mac = (sys.platform == "darwin")
    self.cppflags_start = os.environ.get("CPPFLAGS", "")
    self.ldflags_start = os.environ.get("LDFLAGS", "")
    self.build_cctbx_dependencies()

  def build_cctbx_dependencies (self) :
    print >> self.log, "*** Building dependencies first ***"
    self.print_sep()
    os.chdir(self.tmp_dir)
    options = self.options
    if (True) : #not options.use_system_python) :
      self.build_python()
    else :
      pass # TODO ???
    assert os.path.exists(self.python_exe), self.python_exe
    self.update_paths()
    if (not options.basic) :
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=NUMPY_PKG,
        pkg_name_simple="numpy",
        confirm_import_module="numpy")
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=BIOPYTHON_PKG,
        pkg_name_simple="BioPython",
        confirm_import_module="Bio")
      if (options.build_scipy) :
        # XXX SciPy requires a Fortran compiler, so it's not built by default
        self.build_python_module_simple(
          pkg_url=BASE_CCI_PKG_URL,
          pkg_name=SCIPY_PKG,
          pkg_name_simple="SciPy",
          confirm_import_module="scipy")
    if (options.labelit) or (options.build_gui) or (options.build_all) :
      self.build_imaging()
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=REPORTLAB_PKG,
        pkg_name_simple="reportlab",
        confirm_import_module="reportlab")
    if (options.xia2) or (options.build_all) :
      self.build_hdf5()
    if (options.build_gui) or (options.build_all) :
      self.build_wxpython_dependencies()
      self.build_wxpython()
      self.build_misc()
    print >> self.log, "Dependencies finished building."

  def call (self, args, log=None) :
    if (log is None) : log = self.log
    return call(args, log=log)

  def chdir (self, dir_name, log=None) :
    if (log is None) : log = self.log
    print >> log, "cd \"%s\"" % dir_name
    os.chdir(dir_name)

  def print_sep (self, char="-") :
    print >> self.log, ""
    print >> self.log, char*80
    print >> self.log, ""

  def start_building_package (self, pkg_name) :
    install_log = os.path.join(self.tmp_dir, pkg_name + "_install_log")
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
        print >> self.log, "    using %s" % static_file
        return static_file
    if (os.path.exists(pkg_name)) :
      print >> self.log, "    using ./%s" % pkg_name
      return os.path.join(self.tmp_dir, pkg_name)
    else :
      if (self.options.no_download) :
        raise RuntimeError(("Package '%s' not found on local filesystems.  ") %
          pkg_name)
      full_url = "%s/%s" % (pkg_url, pkg_name)
      self.log.write("    downloading from %s : " % pkg_url)
      f = open(pkg_name, "wb")
      data = urllib2.urlopen(full_url).read()
      assert (len(data) > 0), pkg_name
      self.log.write("%d bytes\n" % len(data))
      self.log.flush()
      f.write(data)
      f.close()
      return os.path.join(self.tmp_dir, pkg_name)

  def build_python (self) :
    log = self.start_building_package("Python")
    os.chdir(self.tmp_dir)
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
    self.python_exe = os.path.abspath(
      os.path.join(self.base_dir, "bin", "python"))
    # just an arbitrary import (with .so)
    self.verify_python_module("Python", "socket")
    self.print_sep()

  def update_paths (self) :
    os.environ["PATH"] = ("%s/bin:" % self.base_dir) + os.environ['PATH']
    lib_paths = [ os.path.join(self.base_dir, "lib") ]
    if (sys.platform == "darwin") :
      lib_paths.append("%s/base/Python.framework/Versions/Current/lib" %
        self.base_dir)
      if ("DYLD_LIBRARY_PATH" in os.environ) :
        lib_paths.append(os.environ["DYLD_LIBRARY_PATH"])
      os.environ['DYLD_LIBRARY_PATH'] = ":".join(lib_paths)
    else :
      if ("LD_LIBRARY_PATH" in os.environ) :
        lib_paths.append(os.environ["LD_LIBRARY_PATH"])
      os.environ['LD_LIBRARY_PATH'] = ":".join(lib_paths)
    inc_dir = os.path.join(self.base_dir, "include")
    if (not os.path.isdir(inc_dir)) :
      os.mkdir(inc_dir)
    self.include_dirs.append(inc_dir)
    self.lib_dirs.append(lib_paths[0])

  def set_cppflags_ldflags_tmp (self) :
    inc_paths = [ "-I\"%s\"" % p for p in self.include_dirs ]
    lib_paths = [ "-L\"%s\"" % p for p in self.lib_dirs ]
    os.environ['CPPFLAGS'] = "%s %s" % (" ".join(inc_paths),
      self.cppflags_start)
    os.environ['LDFLAGS'] = "%s %s" % (" ".join(lib_paths), self.ldflags_start)

  def restore_cppflags_ldflags (self) :
    os.environ['CPPFLAGS'] = self.cppflags_start
    os.environ['LDFLAGS'] = self.ldflags_start

  def verify_python_module (self, pkg_name_simple, module_name) :
    os.chdir(self.tmp_dir) # very important for import to work!
    self.log.write("  verifying %s installation..." % pkg_name_simple)
    self.call("%s -c 'import %s'" % (self.python_exe, module_name))
    print >> self.log, " OK"

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
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    if (callback_before_build is not None) :
      assert callback_before_build(pkg_log), pkg_name
    debug_flag = ""
    if (self.options.debug) :
      debug_flag = "--debug"
    self.call("%s setup.py build %s" % (self.python_exe, debug_flag),
      log=pkg_log)
    self.call("%s setup.py install" % self.python_exe, log=pkg_log)
    if (callback_after_build is not None) :
      assert callback_after_build(pkg_log), pkg_name
    os.chdir(self.tmp_dir)
    if (confirm_import_module is not None) :
      self.verify_python_module(pkg_name_simple, confirm_import_module)
    self.print_sep()

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

  def build_compiled_package_simple (self,
      pkg_url,
      pkg_name,
      pkg_name_simple) :
    os.chdir(self.tmp_dir)
    pkg_log = self.start_building_package(pkg_name_simple)
    pkg = self.fetch_package(pkg_url, pkg_name)
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    self.call("./configure --prefix=\"%s\"" % self.base_dir, log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)
    self.print_sep()

  def build_hdf5 (self) :
    os.chdir(self.tmp_dir)
    pkg_log = self.start_building_package("HDF5")
    pkg = self.fetch_package(BASE_XIA_PKG_URL, HDF5_PKG)
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    print >> pkg_log, "Building base HDF5 library..."
    self.call("./configure --prefix=\"%s\"" % self.base_dir, log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)
    print >> pkg_log, "Building h5py..."
    pkg = self.fetch_package(BASE_XIA_PKG_URL, H5PY_PKG)
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    self.call("%s setup.py build --hdf5=\"%s\"" % (self.python_exe,
      self.base_dir), log=pkg_log)
    self.call("%s setup.py install" % (self.python_exe), log=pkg_log)
    self.call("%s setup.py test" % self.python_exe, log=pkg_log)
    self.verify_python_module("h5py", "h5py")
    self.print_sep()

  def build_wxpython_dependencies (self) :
    self.build_compiled_package_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=FREETYPE_PKG,
      pkg_name_simple="Freetype")
    self.include_dirs.append(os.path.join(self.base_dir, "include",
      "freetype2"))
    if (sys.platform == "darwin") :
      pass
    else :
      pass # TODO add roughly a dozen more libraries here...

  def build_wxpython (self) :
    os.chdir(self.tmp_dir)
    pkg_log = self.start_building_package("wxPython")
    pkg_name = WXPYTHON_PKG
    # XXX we don't entirely trust wxPython-2.9, but it would be preferrable for
    # the future to use a single version instead
    cocoa = False
    if (self.flag_is_mac) :
      if (detect_osx_version() >= 10) :
        print >> self.log, "  running OS 10.6 or later, switching to wx 2.9"
        pkg_name = WXPYTHON_DEV_PKG
        cocoa = True
    pkg = self.fetch_package(BASE_CCI_PKG_URL, pkg_name)
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    # Stage 1: build wxWidgets libraries
    config_opts = [
      "--disable-mediactrl",
      "--with-opengl",
      "--prefix=\"%s\"" % self.base_dir,
      "--disable-unicode", # FIXME
    ]
    if (self.options.debug) :
      config_opts.extend(["--disable-optimize", "--enable-debug"])
      if (self.flag_is_linux) :
        config_opts.append("--disable-debug_gdb")
    else :
      config_opts.extend(["--enable-optimize", "--disable-debugreport"])
    if (cocoa) :
      config_opts.extend(["--with-osx_cocoa", "--enable-monolithic"])
    elif (self.flag_is_mac) :
      config_opts.extend(["--with-mac", "--enable-monolithic"])
    elif (self.flag_is_linux) :
      config_opts.extend([
        "--with-gtk",
        "--with-gtk-prefix=\"%s\"" % self.base_dir,
        "--with-gtk-exec-prefix=\"%s\"" % os.path.join(self.base_dir, "lib"),
        "--enable-graphics_ctx",
      ])
    self.set_cppflags_ldflags_tmp()
    print >> self.log, "  building wxWidgets with options:"
    for opt in config_opts :
      print >> self.log, "    %s" % opt
    self.call("./configure %s" % " ".join(config_opts), log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    if (not cocoa) : # XXX ???
      self.call("make -j %d -C contrib/src/stc" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)
    # Stage 2: build wxPython itself
    wxpy_build_opts = [
      "BUILD_GLCANVAS=1",
      "BUILD_STC=0",
      "BUILD_GIZMOS=0",
      "BUILD_DLLWIDGET=0",
      "UNICODE=0", # FIXME
    ]
    if (cocoa) :
      os.environ['CFLAGS'] = "-arch x86_64"
      wxpy_build_opts.extend(["WXPORT=osx_cocoa"])
    else :
      wxpy_build_opts.append("BUILD_OGL=0")
    self.chdir("wxPython", log=pkg_log)
    debug_flag = ""
    if (self.options.debug) :
      debug_flag = "--debug"
    print >> self.log, "  building wxPython with options:"
    for opt in wxpy_build_opts :
      print >> self.log, "    %s" % opt
    self.call("%s setup.py %s build %s" % (self.python_exe,
      " ".join(wxpy_build_opts), debug_flag), log=pkg_log)
    self.call("%s setup.py %s install" % (self.python_exe,
      " ".join(wxpy_build_opts)), log=pkg_log)
    self.restore_cppflags_ldflags()
    self.verify_python_module("wxPython", "wx")
    self.print_sep()

  def build_misc (self) :
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=PYRTF_PKG,
      pkg_name_simple="PyRTF",
      confirm_import_module="PyRTF")
    self.set_cppflags_ldflags_tmp()
    # TODO we are patching the source in the CCTBX builds (but not the Phenix
    # installer) to force it to use the correct backend - not sure if this
    # is truly necessary or if there's a cleaner way to do it
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=MATPLOTLIB_PKG,
      pkg_name_simple="Matplotlib",
      confirm_import_module="matplotlib")
    self.restore_cppflags_ldflags()

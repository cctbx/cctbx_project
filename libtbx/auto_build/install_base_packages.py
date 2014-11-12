
"""
Automated build of CCTBX dependencies for Linux and Mac platforms.
This script will download and install the current Python distribution provided
by LBL, plus any packages required for various CCTBX-dependent apps to
function.  In the future this will be used as the core of the CCI nightly build
system and the Phenix installer.
"""

from __future__ import division
from installer_utils import *
from package_defs import *
from optparse import OptionParser
import os.path as op
import sys
# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)


class installer (object) :
  # XXX various defaults go here, to be overridden in subclasses
  build_xia2_dependencies = False
  build_gui_dependencies = False
  build_labelit_dependencies = False
  build_dials_dependencies = False

  def __init__ (self, args=None, packages=None, log=sys.stdout) :
    assert (sys.platform in ["linux2", "linux3", "darwin"])
    check_python_version()
    self.log = log
    print >> log, """
  ****************************************************************************
                 Automated CCTBX dependencies build script
                 report problems to cctbx-dev@cci.lbl.gov
  ****************************************************************************
"""
    dist_dir = op.dirname(op.dirname(op.dirname(__file__)))
    parser = OptionParser()
    parser.add_option("--build_dir", dest="build_dir", action="store",
      help="Build directory", default=os.getcwd())
    parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
      help="Temporary directory",
      default=op.join(os.getcwd(), "build_tmp"))
    parser.add_option("-p", "--nproc", dest="nproc", action="store",
      type="int", help="Number of processors", default=1)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
      help="Verbose output", default=False)
    parser.add_option("--pkg_dir", dest="pkg_dirs", action="append",
      help="Directory with source packages", default=None)
    parser.add_option("--basic", dest="basic", action="store_true",
      help="Only install basic prerequisites", default=False)
    parser.add_option("--labelit", dest="labelit", action="store_true",
      help="Build LABELIT dependencies",
      default=self.build_labelit_dependencies)
    parser.add_option("--xia2", dest="xia2", action="store_true",
      help="Build xia2 dependencies", default=self.build_xia2_dependencies)
    parser.add_option("--dials", dest="dials", action="store_true",
      help="Build DIALS dependencies", default=self.build_dials_dependencies)
    parser.add_option("--gui", dest="build_gui", action="store_true",
      help="Build GUI dependencies", default=self.build_gui_dependencies)
    parser.add_option("--with-python", dest="with_python",
      help="Use specified Python interpreter")
    parser.add_option("--with-system-python", dest="with_system_python",
      help="Use the system Python interpreter", action="store_true")
    parser.add_option("-a", "--all", dest="build_all", action="store_true",
      help="Build all recommended dependencies", default=False)
    parser.add_option("--scipy", dest="build_scipy", action="store_true",
      help="Build SciPy (requires Fortran compiler)", default=False)
    parser.add_option("--sphinx", dest="build_sphinx", action="store_true",
      help="Build Sphinx (and numpydoc) for generating documentation")
    parser.add_option("--ipython", dest="build_ipython", action="store_true",
      help="Build IPython", default=False)
    parser.add_option("--pyopengl", dest="build_pyopengl", action="store_true",
      help="Build PyOpenGL", default=False)
    parser.add_option("-g", "--debug", dest="debug", action="store_true",
      help="Build in debugging mode", default=False)
    parser.add_option("--no-download", dest="no_download", action="store_true",
      help="Use only local packages (no downloads)", default=False)
    parser.add_option("--python-shared", dest="python_shared",
      action="store_true", default=False,
      help="Compile Python as shared library (Linux only)")
    parser.add_option("--skip-if-exists", action="store_true",
      help="Exit if build_dir/base exists; used with automated builds.")
    options, args = parser.parse_args(args)
    # basic setup
    self.tmp_dir = options.tmp_dir
    #self.src_dir = options.src_dir
    self.build_dir = options.build_dir
    self.pkg_dirs = options.pkg_dirs
    self.base_dir = op.join(self.build_dir, "base")
    if options.skip_if_exists and os.path.exists(self.base_dir):
      print >> log, "Base directory already exists and --skip-if-exists set; exiting."
      return
    print >> log, "Setting up directories..."
    for dir_name in [self.tmp_dir,self.build_dir,self.base_dir] : #self.src_dir
      if (not op.isdir(dir_name)) :
        print >> log, "  creating %s" % dir_name
        os.makedirs(dir_name)
    self.nproc = options.nproc
    self.verbose = options.verbose
    self.options = options
    self.have_freetype = False
    self.include_dirs = []
    self.lib_dirs = []
    self.flag_is_linux = sys.platform.startswith("linux")
    self.flag_is_mac = (sys.platform == "darwin")
    self.cppflags_start = os.environ.get("CPPFLAGS", "")
    self.ldflags_start = os.environ.get("LDFLAGS", "")
    # Which Python interpreter:
    self.python_exe = options.with_python or None
    if options.with_system_python:
      self.python_exe = sys.executable
    # configure package download
    self.fetch_package = fetch_packages(
      dest_dir=self.tmp_dir,
      log=log,
      pkg_dirs=options.pkg_dirs,
      no_download=options.no_download)

    # Use a specific Python interpreter if provided.
    packages = []
    if self.python_exe:
      self.set_python(self.python_exe)
    else:
      packages += ['python']
    # Always build hdf5.
    packages += ['hdf5']
    # Add specified packages.
    if not options.basic:
      for env_var in ["BLAS","ATLAS","LAPACK"]:
        os.environ[env_var] = "None"
      packages += ['numpy', 'BioPython']
    if options.build_scipy:
      packages += ['scipy']
    if options.build_gui or options.build_all or options.labelit:
      if sys.platform == "darwin":
        packages += ['py2app']
      packages += ['Imaging', 'reportlab']
    if options.build_gui or options.build_all:
      packages += ['gui']
    if options.build_sphinx:
      packages += ['sphinx']
    if options.build_ipython:
      packages += ['ipython']
    if options.build_pyopengl:
      packages += ['PyOpenGL']
    # Or use specified packages if provided.
    if args:
      args = args[1:]
    packages = set(args or packages)
    self.build_dependencies(packages=packages)
    # On Mac OS X all of the Python-related executables located in base/bin
    # are actually symlinks to absolute paths inside the Python.framework, so
    # we replace them with symlinks to relative paths.
    if (sys.platform == "darwin") :
      print >> log, "Regenerating symlinks with relative paths..."
      regenerate_relative_symlinks(op.join(self.base_dir, "bin"), log=log)

  def build_dependencies(self, packages=None):
    packages = packages or []
    print >> self.log, "Building dependencies: %s"%(", ".join(packages))
    self.print_sep()
    os.chdir(self.tmp_dir)
    options = self.options
    if 'python' in packages:
      self.build_python()
    if 'numpy' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=NUMPY_PKG,
        pkg_name_label="numpy",
        confirm_import_module="numpy")
    if 'BioPython' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=BIOPYTHON_PKG,
        pkg_name_label="BioPython",
        confirm_import_module="Bio")
    if 'scipy' in packages:
      # requires Fortran compiler.
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=SCIPY_PKG,
        pkg_name_label="SciPy",
        confirm_import_module="scipy")
    if 'py2app' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=PY2APP_PKG,
        pkg_name_label="py2app",
        confirm_import_module="py2app")
    if 'Imaging' in packages:
      self.build_imaging()
    if 'reportlab' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=REPORTLAB_PKG,
        pkg_name_label="reportlab",
        confirm_import_module="reportlab")
    if 'hdf5' in packages:
      self.build_hdf5()
    if 'gui' in packages:
      self.build_wxpython_dependencies()
      self.build_wxpython()
      self.build_misc()
    if 'sphinx' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=SPHINX_PKG,
        pkg_name_label="Sphinx",
        confirm_import_module="sphinx")
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=NUMPYDOC_PKG,
        pkg_name_label="numpydoc",
        confirm_import_module=None)
    if 'ipython' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=IPYTHON_PKG,
        pkg_name_label="IPython",
        confirm_import_module="IPython")
    if 'PyOpenGL' in packages:
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=PYOPENGL_PKG,
        pkg_name_label="PyOpenGL",
        confirm_import_module="OpenGL")
    print >> self.log, "Dependencies finished building."

  def call (self, args, log=None, **kwargs) :
    if (log is None) : log = self.log
    return call(args, log=log, **kwargs)

  def configure_and_build (self, config_args=(), log=None, make_args=()) :
    self.call("./configure %s" % " ".join(list(config_args)), log=log)
    self.call("make -j %d %s" % (self.nproc, " ".join(list(make_args))), log=log)
    self.call("make install", log=log)

  def chdir (self, dir_name, log=None) :
    if (log is None) : log = self.log
    print >> log, "cd \"%s\"" % dir_name
    os.chdir(dir_name)

  def touch_file (self, file_name) :
    f = open(file_name, "w")
    f.write("")
    f.close()
    assert op.isfile(file_name)

  def print_sep (self, char="-") :
    print >> self.log, ""
    print >> self.log, char*80
    print >> self.log, ""

  def start_building_package (self, pkg_name) :
    os.chdir(self.tmp_dir)
    install_log = op.join(self.tmp_dir, pkg_name + "_install_log")
    print >> self.log, "Installing %s..." % pkg_name
    print >> self.log, "  log file is %s" % install_log
    return open(install_log, "w")

  def patch_src (self, src_file, target, replace_with, output_file=None) :
    if isinstance(target, str) :
      assert isinstance(replace_with, str)
      target = [ target ]
      replace_with = [ replace_with ]
    assert len(target) == len(replace_with)
    in_file = src_file
    if (output_file is None) :
      in_file += ".dist"
      os.rename(src_file, in_file)
    src_in = open(in_file)
    if (output_file is None) :
      output_file = src_file
    src_out = open(output_file, "w")
    for line in src_in.readlines() :
      for target_str, replacement_str in zip(target, replace_with) :
        line = line.replace(target_str, replacement_str)
      src_out.write(line)
    src_in.close()
    src_out.close()

  def fetch_untar_and_chdir (self, pkg_name, pkg_url=None, log=None) :
    if (log is None) : log = self.log
    pkg = self.fetch_package(pkg_name=pkg_name, pkg_url=pkg_url)
    pkg_dir = untar(pkg, log=log)
    os.chdir(pkg_dir)

  def build_python (self) :
    log = self.start_building_package("Python")
    os.chdir(self.tmp_dir)
    python_tarball = self.fetch_package(PYTHON_PKG)
    python_dir = untar(python_tarball)
    self.chdir(python_dir, log=log)
    if (sys.platform == "darwin") :
      configure_args = [
        "--prefix=\"%s\"" % self.base_dir,
        "--enable-framework=\"%s\"" % self.base_dir,]
      self.call("./configure %s" % " ".join(configure_args), log=log)
      self.call("make -j %d" % self.nproc, log=log)
      targets = "bininstall libinstall libainstall inclinstall sharedinstall"
      # XXX don't parallelize here - sometimes breaks on Mac
      self.call("make %s"% (targets), log=log)
      self.chdir("Mac", log=log)
      self.call("make install_Python install_pythonw install_versionedtools",
        log=log)
      self.chdir(python_dir, log=log)
      self.call("make frameworkinstallstructure frameworkinstallmaclib \
        frameworkinstallunixtools FRAMEWORKUNIXTOOLSPREFIX=\"%s\"" %
        self.base_dir, log=log)
    else :
      # Linux
      # Ian: Setup build to use rpath $ORIGIN to find libpython.so.
      # Also note that I'm using shell=False and passing args as list
      #   ... to minimize opportunities for mucking up the "\$$"
      configure_args = ["--prefix", self.base_dir]
      if (self.options.python_shared):
        configure_args.append("--enable-shared")
        configure_args.append("LDFLAGS=-Wl,-rpath=\$$ORIGIN/../lib")
      self.call([os.path.join(python_dir, 'configure')] + configure_args,
        log=log,
        cwd=python_dir,
        shell=False)
      self.call('make -j %s install'%(self.nproc), log=log, cwd=python_dir)
      self.call('make install', log=log, cwd=python_dir)
    python_exe = op.abspath(op.join(self.base_dir, "bin", "python"))
    self.set_python(op.abspath(python_exe))
    self.print_sep()
    log.close()

  def set_python(self, python_exe):
    print >> self.log, "Using Python: %s"%python_exe
    self.python_exe = python_exe
    # Just an arbitrary import (with .so)
    self.verify_python_module("Python", "socket")

    # Check python version >= 2.5.
    python_version = check_output([self.python_exe, '-c', 'import sys; print "%s.%s.%s"%(sys.version_info[0], sys.version_info[1], sys.version_info[2])'])
    python_version = python_version.strip()
    try:
      check_output([self.python_exe, '-c', 'import sys; assert sys.version_info[0] == 2; assert sys.version_info[1] >= 5'])
    except Exception, e:
      print >> self.log, """
Error: Python 2.5 or higher required. Python 3 is not supported.
Found Python version:
  %s
"""%python_version
      raise e

    # Check that we have write access to site-packages dir
    # by creating a temporary file with write permissions.
    # Open with tempfile to auto-handle unlinking.
    import tempfile
    site_packages = check_output([self.python_exe, '-c', 'from distutils.sysconfig import get_python_lib; print(get_python_lib())'])
    site_packages = site_packages.strip()
    # print >> self.log,  "Checking for write permissions:", site_packages
    try:
      with tempfile.TemporaryFile(dir=site_packages) as f:
        pass
    except OSError, e:
      print >> self.log, """
Error: You don't appear to have write access to
the Python site-packages directory:
  %s
Installation of Python packages may fail.
      """%site_packages
      raise e
    # Update paths.
    self.update_paths()

  def update_paths(self):
    os.environ["PATH"] = ("%s/bin:" % self.base_dir) + os.environ['PATH']
    lib_paths = [ op.join(self.base_dir, "lib") ]
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
    inc_dir = op.join(self.base_dir, "include")
    if (not op.isdir(inc_dir)) :
      os.mkdir(inc_dir)
    self.include_dirs.append(inc_dir)
    self.lib_dirs.append(lib_paths[0])

  def set_cppflags_ldflags_tmp (self) :
    # XXX ideally we would like to quote the paths to allow for spaces, but
    # the compiler doesn't like this
    inc_paths = [ "-I%s" % p for p in self.include_dirs ]
    lib_paths = [ "-L%s" % p for p in self.lib_dirs ]
    os.environ['CPPFLAGS'] = "%s %s" % (" ".join(inc_paths),
      self.cppflags_start)
    os.environ['LDFLAGS'] = "%s %s" % (" ".join(lib_paths), self.ldflags_start)

  def restore_cppflags_ldflags (self) :
    os.environ['CPPFLAGS'] = self.cppflags_start
    os.environ['LDFLAGS'] = self.ldflags_start

  def verify_python_module (self, pkg_name_label, module_name) :
    os.chdir(self.tmp_dir) # very important for import to work!
    self.log.write("  verifying %s installation..." % pkg_name_label)
    self.call("%s -c 'import %s'" % (self.python_exe, module_name))
    print >> self.log, " OK"

  def build_python_module_simple (self,
      pkg_url,
      pkg_name,
      pkg_name_label,
      callback_before_build=None,
      callback_after_build=None,
      confirm_import_module=None) :
    pkg_log = self.start_building_package(pkg_name_label)
    pkg = self.fetch_untar_and_chdir(pkg_name=pkg_name, pkg_url=pkg_url,
      log=pkg_log)
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
      self.verify_python_module(pkg_name_label, confirm_import_module)
    self.print_sep()

  def build_imaging (self, patch_src=True) :
    self.build_freetype()
    def patch_imaging_src (out) :
      print >> out, "  patching libImaging/ZipEncode.c"
      self.patch_src(src_file="libImaging/ZipEncode.c",
                     target="Z_DEFAULT_COMPRESSION",
                     replace_with="Z_BEST_SPEED")
      # point to freetype installation
      site_file = open("setup_site.py", "w")
      site_file.write('FREETYPE_ROOT = "%s", "%s"\n' %
        (op.join(self.base_dir, "lib"), op.join(self.base_dir, "include")))
      return True
    callback = None
    if (patch_src) :
      callback = patch_imaging_src
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=IMAGING_PKG,
      pkg_name_label="Imaging",
      callback_before_build=patch_imaging_src,
      confirm_import_module="Image")

  def build_compiled_package_simple (self,
      pkg_name,
      pkg_name_label,
      pkg_url=None) :
    pkg_log = self.start_building_package(pkg_name_label)
    pkg = self.fetch_untar_and_chdir(pkg_name=pkg_name, pkg_url=pkg_url,
      log=pkg_log)
    self.configure_and_build(
      config_args=["--prefix=\"%s\"" % self.base_dir,],
      log=pkg_log)
    self.print_sep()

  def build_hdf5 (self) :
    pkg_log = self.start_building_package("HDF5")
    self.fetch_untar_and_chdir(pkg_name=HDF5_PKG, pkg_url=BASE_XIA_PKG_URL,
      log=pkg_log)
    print >> pkg_log, "Building base HDF5 library..."
    make_args = []
    # XXX the HDF5 library uses '//' for comments, which will break if the
    # compiler doesn't support C99 by default.  for some bizarre reason this
    # includes certain (relatively new) versions of gcc...
    if sys.platform.startswith("linux") :
      make_args.append("CFLAGS=\"--std=c99\"")
    self.configure_and_build(
      config_args=["--prefix=\"%s\"" % self.base_dir,],
      log=pkg_log,
      make_args=make_args)
    print >> pkg_log, "Building h5py..."
    os.chdir(self.tmp_dir)
    self.fetch_untar_and_chdir(pkg_url=BASE_XIA_PKG_URL, pkg_name=H5PY_PKG,
      log=pkg_log)
    self.call("%s setup.py build --hdf5=\"%s\"" % (self.python_exe,
      self.base_dir), log=pkg_log)
    self.call("%s setup.py install" % (self.python_exe), log=pkg_log)
    # FIXME this appears to fail on CentOS 5 (gcc 4.1.2)
    #self.call("%s setup.py test" % self.python_exe, log=pkg_log)
    self.verify_python_module("h5py", "h5py")
    self.print_sep()

  def build_freetype (self) :
    assert (not self.have_freetype)
    self.build_compiled_package_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=FREETYPE_PKG,
      pkg_name_label="Freetype")
    self.include_dirs.append(op.join(self.base_dir, "include",
      "freetype2"))
    self.have_freetype = True

  def build_wxpython_dependencies (self) :
    if (not self.have_freetype) :
      self.build_freetype()
    self.build_compiled_package_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=LIBPNG_PKG,
      pkg_name_label="libpng")
    if (sys.platform == "darwin") :
      # XXX nuke libpng dylibs - I forget why we need to do this
      lib_dir = op.join(self.base_dir, "lib")
      for file_name in os.listdir(lib_dir) :
        if ("png" in file_name.lower()) and (file_name.endswith(".dylib")) :
          full_path = op.join(lib_dir, file_name)
          os.remove(full_path)
    else :
      self.build_wxpython_dependencies_linux()

  def build_wxpython_dependencies_linux (self) :
    pkg_config_dir = op.join(self.base_dir, "lib", "pkgconfig")
    if (not op.isdir(pkg_config_dir)) :
      os.makedirs(pkg_config_dir)
    pkg_config_paths = [pkg_config_dir] + os.environ.get("PKG_CONFIG_PATH", [])
    os.environ['PKG_CONFIG_PATH'] = ":".join(pkg_config_paths)
    prefix_arg = "--prefix=\"%s\"" % self.base_dir
    # gettext
    pkg_log = self.start_building_package("gettext")
    self.fetch_untar_and_chdir(pkg_name=GETTEXT_PKG, log=pkg_log)
    os.chdir("gettext-runtime")
    self.set_cppflags_ldflags_tmp()
    gettext_conf_args = [ prefix_arg, "--disable-java", "--disable-csharp",
      "--disable-intl-java", "--disable-gcj" ]
    self.configure_and_build(config_args=gettext_conf_args, log=pkg_log)
    self.print_sep()
    # glib
    glib_log = self.start_building_package("glib")
    self.fetch_untar_and_chdir(pkg_name=GLIB_PKG, log=glib_log)
    msgfmt_bin = op.join(self.base_dir, "bin", "msgfmt")
    gettext_bin = op.join(self.base_dir, "bin", "xgettext")
    def touch_bin_files () :
      self.touch_file(msgfmt_bin)
      self.touch_file(gettext_bin)
      self.call("chmod 744 \"%s\"" % msgfmt_bin)
      self.call("chmod 744 \"%s\"" % gettext_bin)
    touch_bin_files()
    make_args = []
    if (self.options.debug) : make_args.append("CFLAGS=-g")
    self.configure_and_build(
      config_args=[prefix_arg],
      log=glib_log,
      make_args=make_args)
    os.remove(msgfmt_bin)
    os.remove(gettext_bin)
    self.print_sep()
    # expat
    expat_log = self.start_building_package("expat")
    self.fetch_untar_and_chdir(pkg_name=EXPAT_PKG, log=expat_log)
    self.configure_and_build(config_args=[prefix_arg], log=expat_log)
    header_files = ["./lib/expat_external.h", "./lib/expat.h"]
    for header in header_files :
      self.call("./conftools/install-sh -c -m 644 %s \"%s\"" % (header,
        op.join(self.base_dir, "include")), log=expat_log)
    self.print_sep()
    # fontconfig
    fc_log = self.start_building_package("fontconfig")
    self.fetch_untar_and_chdir(pkg_name=FONTCONFIG_PKG, log=fc_log)
    if (not op.isdir(op.join(self.base_dir, "share", "fonts"))) :
      os.makedirs(op.join(self.base_dir, "share", "fonts"))
    if (not op.isdir(op.join(self.base_dir,"etc","fonts"))) :
      os.makedirs(op.join(self.base_dir,"etc","fonts"))
    fc_config_args = [ prefix_arg,
      "--disable-docs",
      "--with-expat-includes=\"%s\"" % op.join(self.base_dir, "include"),
      "--with-expat-lib=\"%s\"" % op.join(self.base_dir, "lib"),
      "--with-add-fonts=\"%s\"" % op.join(self.base_dir,"share","fonts"),
      "--with-confdir=\"%s\"" % op.join(self.base_dir,"etc","fonts"),
      "--with-docdir=\"%s\"" % op.join(self.base_dir, "doc"),
      "--with-freetype-config=freetype-config", ]
    self.configure_and_build(config_args=fc_config_args, log=fc_log)
    self.print_sep()
    # render, xrender, xft
    for pkg, name in zip([RENDER_PKG,XRENDER_PKG,XFT_PKG],
                         ["render","Xrender","Xft"]) :
      self.build_compiled_package_simple(pkg_name=pkg, pkg_name_label=name)
    # pixman
    pix_log = self.start_building_package("pixman")
    self.fetch_untar_and_chdir(pkg_name=PIXMAN_PKG, log=pix_log)
    self.configure_and_build(
      config_args=[ prefix_arg, "--disable-gtk" ],
      log=pix_log,
      make_args=make_args) # re-used from above
    self.print_sep()
    # cairo, pango, atk
    for pkg, name in zip([CAIRO_PKG, PANGO_PKG, ATK_PKG],
                         ["cairo", "pango", "atk"]) :
      self.build_compiled_package_simple(pkg_name=pkg, pkg_name_label=name)
    # tiff
    tiff_log = self.start_building_package("tiff")
    self.fetch_untar_and_chdir(pkg_name=TIFF_PKG, log=tiff_log)
    os.environ['MANSCHEME'] = "bsd-source-cat"
    os.environ['DIR_MAN'] = op.join(self.base_dir, "man")
    config_args = [ prefix_arg, "--noninteractive", "--with-LIBGL=no",
                    "--with-LIBIMAGE=no" ]
    self.configure_and_build(
      config_args=config_args,
      log=tiff_log,
      make_args=make_args)
    self.print_sep()
    # gtk+
    touch_bin_files()
    gtk_log = self.start_building_package("gtk+")
    self.fetch_untar_and_chdir(pkg_name=GTK_PKG, log=gtk_log)
    gtk_config_args = [
      prefix_arg,
      "--disable-cups",
      "--without-libjpeg",
    ]
    self.call("./configure %s" % " ".join(gtk_config_args), log=gtk_log)
    self.call("make -j %d SRC_SUBDIRS='gdk-pixbuf gdk gtk modules'" %
      self.nproc, log=gtk_log)
    self.call("make install SRC_SUBDIRS='gdk-pixbuf gdk gtk modules'",
      log=gtk_log)
    os.remove(msgfmt_bin)
    os.remove(gettext_bin)
    self.print_sep()
    # gtk-engine
    self.build_compiled_package_simple(pkg_name=GTK_ENGINE_PKG,
      pkg_name_label="gtk-engine")
    # fonts
    fonts_log = self.start_building_package("fonts")
    share_dir = op.join(self.base_dir, "share")
    if (not op.isdir(share_dir)) :
      os.makedirs(share_dir)
    pkg = self.fetch_package(pkg_name=FONT_PKG)
    os.chdir(share_dir)
    untar(pkg, log=fonts_log, verbose=True)
    self.print_sep()
    os.chdir(self.tmp_dir)

  def build_wxpython (self) :
    pkg_log = self.start_building_package("wxPython")
    pkg_name = WXPYTHON_PKG
    # XXX we don't entirely trust wxPython-2.9, but it would be preferrable for
    # the future to use a single version instead
    cocoa = False
    if (self.flag_is_mac) :
      if (detect_osx_version() >= 10) :
        print >> self.log, "  running OS 10.6 or later, switching to wx 3.0"
        pkg_name = WXPYTHON_DEV_PKG
        cocoa = True
    pkg = self.fetch_package(pkg_name)
    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    if (self.flag_is_mac and get_os_version() == "10.10") :
      # Workaround wxwidgets 3.0.2 compilation error on Yosemite
      # See:
      #   http://trac.wxwidgets.org/ticket/16329
      #   http://goharsha.com/blog/compiling-wxwidgets-3-0-2-mac-os-x-yosemite/
      print >> self.log, "  patching src/osx/webview_webkit.mm"
      self.patch_src(src_file="src/osx/webview_webkit.mm",
                     target=("#include <WebKit/WebKit.h>",),
                     replace_with=("#include <WebKit/WebKitLegacy.h>",))
    # Stage 1: build wxWidgets libraries
    config_opts = [
      "--disable-mediactrl",
      "--with-opengl",
      "--prefix=\"%s\"" % self.base_dir,
    ]
    if (self.options.debug) :
      config_opts.extend(["--disable-optimize", "--enable-debug"])
      if (self.flag_is_linux) :
        config_opts.append("--disable-debug_gdb")
    else :
      config_opts.extend(["--enable-optimize", "--disable-debugreport"])
    if (cocoa) :
      config_opts.extend(["--with-osx_cocoa", "--enable-monolithic",
        "--with-macosx-version-min=10.6", "--enable-unicode",])
      os.environ["CXXFLAGS"] = \
        "-DMAC_OS_X_VERSION_MIN_REQUIRED=MAC_OS_X_VERSION_10_6"
    elif (self.flag_is_mac) :
      config_opts.extend(["--with-mac", "--enable-monolithic"])
    elif (self.flag_is_linux) :
      config_opts.extend([
        "--with-gtk",
        "--with-gtk-prefix=\"%s\"" % self.base_dir,
        "--with-gtk-exec-prefix=\"%s\"" % op.join(self.base_dir, "lib"),
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
    ]
    if (cocoa) :
      os.environ['CFLAGS'] = "-arch x86_64"
      wxpy_build_opts.extend(["WXPORT=osx_cocoa", "UNICODE=1",])
    else :
      wxpy_build_opts.extend(["BUILD_OGL=0", "UNICODE=0"])
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
      pkg_name_label="PyRTF",
      confirm_import_module="PyRTF")
    self.set_cppflags_ldflags_tmp()
    # TODO we are patching the source to force it to use the correct backend.
    # I'm not sure if this is truly necessary or if there's a cleaner way...
    def patch_matplotlib_src (out) :
      print >> out, "  patching setup.cfg"
      self.patch_src(src_file="setup.cfg.template",
                     output_file="setup.cfg",
                     target=("#backend = Agg", "#basedirlist = /usr"),
                     replace_with=("backend = WXAgg",
                      "basedirlist = /usr, %s" % self.base_dir))
      return True
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=MATPLOTLIB_PKG,
      pkg_name_label="Matplotlib",
      callback_before_build=patch_matplotlib_src,
      confirm_import_module="matplotlib")
    self.restore_cppflags_ldflags()
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=SEND2TRASH_PKG,
      pkg_name_label="send2trash",
      )#confirm_import_module="send2trash")

  # TODO
  def write_dispatcher_include (self) :
    raise NotImplementedError()

def check_wxpython_build_dependencies (log=sys.stderr) :
  try :
    call(["pkg-config", "--version"], log=log)
  except RuntimeError :
    return """
ERROR The GUI components require pkg-config to build
Please install pkg-config to compile these components or use the --no-gui
option to disable compilation of GUI components.
"""
  # TODO pkg-config version check
  try :
    call(["bash", "--version"], log=log)
  except RuntimeError :
    return """
ERROR: The GUI requires bash to be available to build
Please install bash to compile these components, or use the --no-gui option
to disable GUI compilation.
"""
  if (not op.exists("/usr/include/X11/X.h") and
      not op.exists("/usr/X11R6/include/X11/X.h")) :
    return """
ERROR: The X-windows headers appear to be missing
Please install the X11 development packages to compile the GUI components,
or use the --no-gui option to disable GUI compilation.
"""
  return None

if __name__ == "__main__":
  installer(args=sys.argv, log=sys.stdout)

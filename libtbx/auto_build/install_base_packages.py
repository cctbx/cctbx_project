
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
import time
import zipfile

# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)

python_dependencies = {"_ssl" : "Secure Socket Library",
                       "zlib" : "XCode command line tools",
                     }

class installer (object) :
  def __init__ (self, args=None, packages=None, log=sys.stdout) :
    #assert (sys.platform in ["linux2", "linux3", "darwin"])
    check_python_version() # (2,6)
    self.log = log
    print >> log, """
  ****************************************************************************
                 Automated CCTBX dependencies build script
                 report problems to cctbx-dev@cci.lbl.gov
  ****************************************************************************
"""
    dist_dir = op.dirname(op.dirname(op.dirname(__file__)))
    parser = OptionParser()
    # Basic options
    parser.add_option("--build_dir", dest="build_dir", action="store",
      help="Build directory", default=os.getcwd())
    parser.add_option("--tmp_dir", dest="tmp_dir", action="store",
      help="Temporary directory",
      default=op.join(os.getcwd(), "base_tmp"))
    parser.add_option("--pkg_dir", dest="pkg_dirs", action="append",
      help="Directory with source packages", default=None)
    parser.add_option("-p", "--nproc", dest="nproc", action="store",
      type="int", help="Number of processors", default=1)
    parser.add_option("-v", "--verbose", dest="verbose", action="store_true",
      help="Verbose output", default=False)
    parser.add_option("--skip-if-exists", action="store_true",
      help="Exit if build_dir/base exists; used with automated builds.")
    parser.add_option("--no-download", dest="no_download", action="store_true",
      help="Use only local packages (no downloads)", default=False)
    parser.add_option("--download-only", dest="download_only", action="store_true",
      help="Only download missing packages, do not compile", default=False)
    parser.add_option("--skip-base", dest="skip_base", action="store",default="")
    parser.add_option("--python-shared", dest="python_shared",
      action="store_true", default=False,
      help="Compile Python as shared library (Linux only)")
    parser.add_option("--with-python", dest="with_python",
      help="Use specified Python interpreter")
    parser.add_option("--with-system-python", dest="with_system_python",
      help="Use the system Python interpreter", action="store_true")
    parser.add_option("-g", "--debug", dest="debug", action="store_true",
      help="Build in debugging mode", default=False)
    # Package set options.
    parser.add_option("--cctbx", dest="cctbx", action="store_true",
      help="Build CCTBX dependencies")
    parser.add_option("--phenix", dest="phenix", action="store_true",
      help="Build PHENIX dependencies")
    parser.add_option("--labelit", dest="labelit", action="store_true",
      help="Build LABELIT dependencies")
    parser.add_option("--xia2", dest="xia2", action="store_true",
      help="Build xia2 dependencies")
    parser.add_option("--dials", dest="dials", action="store_true",
      help="Build DIALS dependencies")
    parser.add_option("--gui", dest="build_gui", action="store_true",
      help="Build GUI dependencies")
    parser.add_option("-a", "--all", dest="build_all", action="store_true",
      help="Build all recommended dependencies", default=False)
    # Specific add-on packages.
    parser.add_option("--scipy", dest="build_scipy", action="store_true",
      help="Build SciPy (requires Fortran compiler)", default=False)
    parser.add_option("--ipython", dest="build_ipython", action="store_true",
      help="Build IPython", default=False)
    parser.add_option("--wxpython3", dest="use_wxpython3", action="store_true",
      help="Use wxPython3", default=False)

    # Basic setup.
    options, args = parser.parse_args(args)
    self.nproc = options.nproc
    self.verbose = options.verbose
    self.options = options
    self.include_dirs = []
    self.lib_dirs = []
    self.flag_is_linux = sys.platform.startswith("linux")
    self.flag_is_mac = (sys.platform == "darwin")
    self.cppflags_start = os.environ.get("CPPFLAGS", "")
    self.ldflags_start = os.environ.get("LDFLAGS", "")

    # Directory setup.
    self.tmp_dir = options.tmp_dir
    self.build_dir = options.build_dir
    self.pkg_dirs = options.pkg_dirs
    self.base_dir = op.join(self.build_dir, "base")
    self.prefix = "--prefix=\"%s\""%self.base_dir
    if options.skip_if_exists and os.path.exists(self.base_dir):
      print >> log, "Base directory already exists and --skip-if-exists set; exiting."
      return
    print >> log, "Setting up directories..."
    for dir_name in [self.tmp_dir,self.build_dir,self.base_dir]:
      if (not op.isdir(dir_name)) :
        print >> log, "  creating %s" % dir_name
        os.makedirs(dir_name)
    self.check_python_dependencies()

    # Configure package download
    self.fetch_package = fetch_packages(
      dest_dir=self.tmp_dir,
      log=log,
      pkg_dirs=options.pkg_dirs,
      no_download=options.no_download)
    # Shortcut: Extract python for Windows bundled with all preinstalled modules
    if sys.platform == "win32":
      import platform
      if platform.architecture()[0] == '64bit':
        winpythonpkg = WIN64PYTHON_PKG
      else:
        winpythonpkg = WIN32PYTHON_PKG
      self.fetch_package(pkg_name=winpythonpkg, pkg_url=BASE_CCI_PKG_URL)
      winpython = zipfile.ZipFile(os.path.join(self.tmp_dir, winpythonpkg), 'r')
      members = winpython.namelist()
      for zipinfo in members:
        print >> self.log, "extracting", zipinfo
        winpython.extract(zipinfo, path=os.path.join(self.base_dir,'bin'))
      winpython.close()
      # On Windows quit now as all required modules are in the precompiled python package
      return

    # Which Python interpreter:
    self.python_exe = None
    if os.path.exists(os.path.join(self.build_dir, 'base', 'bin', 'python')):
      self.python_exe = os.path.join(self.build_dir, 'base', 'bin', 'python')
    elif options.with_python:
      self.python_exe = options.with_python
    elif options.with_system_python:
      self.python_exe = sys.executable

    if self.python_exe:
      print >> log, "Using Python interpreter: %s" % self.python_exe

    import platform

    if not self.python_exe and 'SuSE' in platform.platform():
      if 'CONFIG_SITE' in os.environ:
        print >> log, 'SuSE detected; clobbering CONFIG_SITE in environ'
        try:
          del(os.environ['CONFIG_SITE'])
        except: # intentional
          pass

    # Set package config.
    pkg_config_dir = op.join(self.base_dir, "lib", "pkgconfig")
    if (not op.isdir(pkg_config_dir) and not options.download_only):
      os.makedirs(pkg_config_dir)
    pkg_config_paths = [pkg_config_dir] + os.environ.get("PKG_CONFIG_PATH", "").split(":")
    os.environ['PKG_CONFIG_PATH'] = ":".join(pkg_config_paths)

    # Set BLAS/ATLAS/LAPACK
    for env_var in ["BLAS","ATLAS","LAPACK"]:
      os.environ[env_var] = "None"

    # Select packages to build.
    packages = self.configure_packages(options)
    # Override and specified packages if provided.
    if len(args) > 1:
      packages = args[1:]

    # Do the work!
    self.check_dependencies(packages=packages)
    self.build_dependencies(packages=packages)

    # On Mac OS X all of the Python-related executables located in base/bin
    # are actually symlinks to absolute paths inside the Python.framework, so
    # we replace them with symlinks to relative paths.
    if self.flag_is_mac and not self.options.download_only:
      print >> log, "Regenerating symlinks with relative paths..."
      regenerate_relative_symlinks(op.join(self.base_dir, "bin"), log=log)

  def configure_packages(self, options):
    packages = []
    # Package groups.
    if options.phenix:
      options.build_gui = True
      options.build_all = True
    if options.dials:
      options.build_gui = True
      options.build_all = True
      packages += ['imaging']
    if options.labelit:
      options.build_gui = True
      options.build_all = True
      packages += ['imaging', 'reportlab']

    # Use a specific Python interpreter if provided.
    if self.python_exe:
      self.set_python(self.python_exe)
    else:
      packages += ['python']

    # Always build hdf5 and numpy.
    packages += ['cython', 'hdf5', 'numpy', 'setuptools', 'pip', 'docutils']
    packages += ["libsvm"]
    # GUI packages.
    if options.build_gui or options.build_all:
      packages += [
        'png',
        'freetype',
        'matplotlib',
        'pyopengl',
        'wxpython',
      ]
      if self.flag_is_mac:
        # Use system libpng.
        packages += ['py2app']
      if self.flag_is_linux:
        packages += [
          'tiff',
          'gettext',
          'glib',
          'expat',
          'fontconfig',
          'render',
          'pixman',
          'cairo',
          'gtk',
          'fonts',
        ]

    # Additional recommended dependencies.
    if options.build_all:
      packages += ['biopython', 'misc', 'sphinx']

    # Non-all packages.
    # Scipy
    if options.build_scipy:
      packages += ['scipy']

    # IPython
    if options.build_ipython:
      packages += ['ipython']

    # Package dependencies.
    if 'imaging' in packages:
      packages += ['freetype']

    return set(packages)

  def call (self, args, log=None, **kwargs) :
    if (log is None) : log = self.log
    return call(args, log=log, verbose=self.verbose, **kwargs)

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

  def check_download_only (self, pkg_name) :
    if self.options.download_only:
      print >> self.log, "  skipping installation of %s (--download-only)" % pkg_name
      return True
    else:
      print >> self.log, "  installing %s..." % pkg_name
      return False

  def patch_src (self, src_file, target, replace_with, output_file=None) :
    from shutil import copymode
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
    copymode(in_file, output_file)

  def untar_and_chdir (self, pkg, log=None) :
    if (log is None) : log = self.log
    pkg_dir = untar(pkg, log=log)
    os.chdir(pkg_dir)

  def check_python_version(self, major=2, minor=6):
    python_version = check_output([
      self.python_exe,
      '-c',
      'import sys; print "%s.%s.%s"%(sys.version_info[0], sys.version_info[1], sys.version_info[2])',
      ])
    python_version = python_version.strip()
    try:
      check_output([
        self.python_exe,
        '-c',
        'import sys; assert sys.version_info[0] == %d; assert sys.version_info[1] >= %d' % (major, minor),
        ])
    except (OSError, RuntimeError), e:
      print >> self.log, """
Error: Python 2.6 or higher required. Python 3 is not supported.
Found Python version:
  %s
"""%python_version
      raise e

  def check_python_dependencies(self):
    for module, desc in python_dependencies.items():
      try:
        self.verify_python_module(module, module)
      except RuntimeError, e:
        print '\n\n\n%s\n\nThis python does not have %s installed\n\n%s\n\n\n' % (
          "*"*80,
          "%s - %s" % (module, desc),
          "*"*80,
        )
        raise e

  def set_python(self, python_exe):
    print >> self.log, "Using Python: %s"%python_exe
    self.python_exe = python_exe
    # Just an arbitrary import (with .so)
    self.verify_python_module("Python", "socket")
    # check that certain modules are avaliable in python
    self.check_python_dependencies()
    # Check python version >= 2.6.
    self.check_python_version()

    # Check that we have write access to site-packages dir
    # by creating a temporary file with write permissions.
    # Open with tempfile to auto-handle unlinking.
    import tempfile
    site_packages = check_output([self.python_exe, '-c', 'from distutils.sysconfig import get_python_lib; print(get_python_lib())'])
    site_packages = site_packages.strip()
    print >> self.log,  "Checking for write permissions:", site_packages
    try:
      f = tempfile.TemporaryFile(dir=site_packages)
    except (OSError, RuntimeError), e:
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
    if self.flag_is_linux:
      if ("LD_LIBRARY_PATH" in os.environ) :
        lib_paths.append(os.environ["LD_LIBRARY_PATH"])
      os.environ['LD_LIBRARY_PATH'] = ":".join(lib_paths)
    inc_dir = op.join(self.base_dir, "include")
    if (not op.isdir(inc_dir)) :
      os.mkdir(inc_dir)
    self.include_dirs.append(inc_dir)
    self.lib_dirs.append(lib_paths[0])

  def set_cppflags_ldflags(self):
    # XXX ideally we would like to quote the paths to allow for spaces, but
    # the compiler doesn't like this
    inc_paths = ["-I%s" % p for p in self.include_dirs]
    lib_paths = ["-L%s" % p for p in self.lib_dirs]
    os.environ['CPPFLAGS'] = "%s %s"%(" ".join(inc_paths), self.cppflags_start)
    os.environ['LDFLAGS'] = "%s %s"%(" ".join(lib_paths), self.ldflags_start)

  def verify_python_module (self, pkg_name_label, module_name) :
    os.chdir(self.tmp_dir) # very important for import to work!
    if hasattr(self, "python_exe"): python_exe = self.python_exe
    else:
      python_exe = sys.executable
    self.log.write("  verifying %s installation in %s" % (
      pkg_name_label,
      python_exe,
    ))
    if sys.platform == "win32": # Windows is picky about using double quotes rather than single quotes
      self.call('%s -c "import %s"' % (python_exe, module_name))
    else:
      self.call("%s -c 'import %s'" % (python_exe, module_name))
    print >> self.log, " OK"

  def workarounds(self):
    '''Look at the directory I am in at the moment and the platform I am
    running on and (perhaps) mess with things to make the build work.'''

    import os
    import platform

    # Ubuntu & Xrender - libtool is broken - replace with system one if
    # installed...

    if 'xrender' in os.path.split(os.getcwd())[-1] and \
      'Ubuntu' in platform.platform() and os.path.exists('libtool'):
      if os.path.exists(os.path.join('/', 'usr', 'bin', 'libtool')):
        self.log.write('Removing xrender libtool; replace with system\n')
        os.remove('libtool')
        os.symlink(os.path.join('/', 'usr', 'bin', 'libtool'), 'libtool')
      else:
        self.log.write('Cannot removing xrender libtool; not installed\n')

    return

  def configure_and_build (self, config_args=(), log=None, make_args=()) :
    self.call("./configure %s" % " ".join(list(config_args)), log=log)
    self.workarounds()
    self.call("make -j %d %s" % (self.nproc, " ".join(list(make_args))), log=log)
    self.call("make install", log=log)

  def build_compiled_package_simple (self,
      pkg_name,
      pkg_name_label,
      pkg_url=None) :
    pkg_log = self.start_building_package(pkg_name_label)
    pkg = self.fetch_package(pkg_name=pkg_name, pkg_url=pkg_url)
    if self.check_download_only(pkg_name): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    self.configure_and_build(
      config_args=[self.prefix],
      log=pkg_log)

  def build_python_module_simple (self,
      pkg_url,
      pkg_name,
      pkg_name_label,
      callback_before_build=None,
      callback_after_build=None,
      confirm_import_module=None) :
    pkg_log = self.start_building_package(pkg_name_label)
    pkg, size = self.fetch_package(pkg_name=pkg_name,
                                   pkg_url=pkg_url,
                                   return_file_and_status=True,
      )
    if self.check_download_only(pkg_name): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
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

  def check_dependencies(self, packages=None):
    packages = packages or []
    if 'scipy' in packages:
      compilers = ['gfortran', 'g77', 'g95', 'f77', 'f95']
      found = []
      for compiler in compilers:
        try:
          check_output(compiler)
          found.append(compiler)
        except RuntimeError, e:
          # Found compiler, but error like "no input files"
          found.append(compiler)
        except OSError, e:
          # Command not found
          pass
      if not found:
        raise Exception("No Fortran compiler found for Scipy. Requires one of: %s"%(", ".join(compilers)))

  def build_dependencies(self, packages=None):
    # Build in the correct dependency order.
    packages = packages or []
    order = [
      'python',
      'numpy',
      'cython',
      'hdf5',
      'setuptools',
      'libsvm',
      'pip',
      'biopython',
      'reportlab',
      'sphinx',
      'docutils',
      'ipython',
      'pyopengl',
      'scipy',
      'py2app',
      'misc',
      # ...
      'png',
      'freetype',
      'matplotlib',
      'imaging',
      # START GUI PACKAGES
      'gettext',
      'glib',
      'expat',
      'fontconfig',
      'render',
      'pixman',
      'tiff',
      'cairo',
      'gtk',
      'fonts',
      'wxpython',
    ]
    self.options.skip_base = self.options.skip_base.split(",")
    packages_order = []
    for i in order:
      if i in packages:
        if i in self.options.skip_base: continue
        packages_order.append(i)

    if self.options.download_only:
      print >> self.log, "Downloading dependencies: %s"%(" ".join(packages_order))
    else:
      print >> self.log, "Building dependencies: %s"%(" ".join(packages_order))

    os.chdir(self.tmp_dir)
    for i in packages_order:
      self.set_cppflags_ldflags() # up-to-date LDFLAGS/CPPFLAGS
      self.print_sep()
      if (i == 'wxpython' and self.options.use_wxpython3):
        i = 'wxpython3'
      if ( hasattr(self, "built_%s" % i) and
           getattr(self, "built_%s" % i)()
         ):
        print >> self.log, "Installation of %s is OK - skipping" % i
        continue
      t0=time.time()
      getattr(self, 'build_%s'%i)()
      print >> self.log, "  package %s took %0.1fs to install" % (
        i,
        time.time()-t0,
        )

    if self.options.download_only:
      print >> self.log, "Dependencies finished downloading."
    else:
      print >> self.log, "Dependencies finished building."

  #######################################################
  ##### Build Individual Packages #######################
  #######################################################

  def build_python (self) :
    log = self.start_building_package("Python")
    os.chdir(self.tmp_dir)
    python_tarball = self.fetch_package(PYTHON_PKG)
    if self.check_download_only(PYTHON_PKG): return
    python_dir = untar(python_tarball)
    self.chdir(python_dir, log=log)
    if self.flag_is_mac:
      configure_args = [
        self.prefix,
        "--enable-framework=\"%s\"" % self.base_dir,]
      self.call("./configure %s" % " ".join(configure_args), log=log)
      # Parallel build of python 2.7 needs to be limited to cope with race conditions
      self.call("make -j %d" % min(self.nproc, 3), log=log)
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
    else:
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
      # Parallel build of python 2.7 needs to be limited to cope with race conditions
      self.call('make -j %s install'% min(3, self.nproc), log=log, cwd=python_dir)
      self.call('make install', log=log, cwd=python_dir)
    python_exe = op.abspath(op.join(self.base_dir, "bin", "python"))
    self.set_python(op.abspath(python_exe))
    log.close()

  def simple_log_parse_test(self, log_filename, line):
    if not os.path.exists(log_filename): return False
    f = file(log_filename, "rb")
    log_lines = f.read()
    f.close()
    if log_lines.find(line)>-1:
      return True
    return False

  def build_setuptools(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=SETUPTOOLS_PKG,
      pkg_name_label="setuptools")

  def build_libsvm(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=LIBSVM_PKG,
      pkg_name_label="libsvm")

  def build_pip(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=PIP_PKG,
      pkg_name_label="pip")

  def build_numpy(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=NUMPY_PKG,
      pkg_name_label="numpy",
      confirm_import_module="numpy")

  def build_docutils(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=DOCUTILS_PKG,
      pkg_name_label="docutils",
      confirm_import_module="docutils")

  def build_biopython(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=BIOPYTHON_PKG,
      pkg_name_label="biopython",
      confirm_import_module="Bio")

  def build_imaging (self, patch_src=True) :
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
      pkg_name_label="imaging",
      callback_before_build=patch_imaging_src,
      confirm_import_module="Image")

  def build_scipy(self):
    # requires Fortran compiler.
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=SCIPY_PKG,
      pkg_name_label="SciPy",
      confirm_import_module="scipy")

  def build_py2app(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=PY2APP_PKG,
      pkg_name_label="py2app",
      confirm_import_module="py2app")

  def build_reportlab(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=REPORTLAB_PKG,
      pkg_name_label="reportlab",
      confirm_import_module="reportlab")

  def build_sphinx(self):
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

  def build_ipython(self):
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=IPYTHON_PKG,
      pkg_name_label="IPython",
      confirm_import_module="IPython")

  def build_pyopengl(self):
      self.build_python_module_simple(
        pkg_url=BASE_CCI_PKG_URL,
        pkg_name=PYOPENGL_PKG,
        pkg_name_label="pyopengl",
        confirm_import_module="OpenGL")

  def built_cython(self):
    log_filename = "cython_install_log"
    line = "Finished processing dependencies for Cython"
    return self.simple_log_parse_test(log_filename, line)

  def build_cython(self):
    self.build_python_module_simple(
      pkg_url=BASE_CYTHON_PKG_URL,
      pkg_name=CYTHON_PKG,
      pkg_name_label="cython",
      confirm_import_module="Cython")

  def built_hdf5(self):
    log_filename = "HDF5_install_log"
    line = "Finished processing dependencies for h5py"
    return self.simple_log_parse_test(log_filename, line)

  def build_hdf5 (self):
    pkg_log = self.start_building_package("HDF5")
    hdf5pkg = self.fetch_package(pkg_name=HDF5_PKG, pkg_url=BASE_HDF5_PKG_URL)
    h5pypkg = self.fetch_package(pkg_name=H5PY_PKG, pkg_url=BASE_H5PY_PKG_URL)
    if ( self.check_download_only(HDF5_PKG) and
         self.check_download_only(H5PY_PKG)
         ):
      return
    self.untar_and_chdir(pkg=hdf5pkg, log=pkg_log)
    print >> pkg_log, "Building base HDF5 library..."
    make_args = []
    # XXX the HDF5 library uses '//' for comments, which will break if the
    # compiler doesn't support C99 by default.  for some bizarre reason this
    # includes certain (relatively new) versions of gcc...
    if self.flag_is_linux:
      make_args.append("CFLAGS=\"--std=c99\"")
    self.configure_and_build(
      config_args=[self.prefix],
      log=pkg_log,
      make_args=make_args)

    # not strictly necessary, but prints 'installing' message
    if self.check_download_only(H5PY_PKG): return
    os.chdir(self.tmp_dir)

    self.untar_and_chdir(pkg=h5pypkg, log=pkg_log)
    self.call("%s setup.py configure --hdf5=\"%s\"" % (self.python_exe,
      self.base_dir), log=pkg_log)
    self.call("%s setup.py build" % (self.python_exe), log=pkg_log)
    self.call("%s setup.py install" % (self.python_exe), log=pkg_log)
    # FIXME this appears to fail on CentOS 5 (gcc 4.1.2)
    #self.call("%s setup.py test" % self.python_exe, log=pkg_log)
    self.verify_python_module("h5py", "h5py")

  def build_freetype (self) :
    self.build_compiled_package_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=FREETYPE_PKG,
      pkg_name_label="Freetype")
    self.include_dirs.append(op.join(self.base_dir, "include", "freetype2"))

  def build_png(self):
    self.build_compiled_package_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=LIBPNG_PKG,
      pkg_name_label="libpng")

  def build_gettext(self):
    # gettext
    pkg_log = self.start_building_package("gettext")
    pkg = self.fetch_package(pkg_name=GETTEXT_PKG)
    if self.check_download_only(GETTEXT_PKG): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    os.chdir("gettext-runtime")
    gettext_conf_args = [self.prefix, "--disable-java", "--disable-csharp",
      "--disable-intl-java", "--disable-gcj"]
    self.configure_and_build(config_args=gettext_conf_args, log=pkg_log)

  def build_glib(self):
    # glib
    pkg_log = self.start_building_package("glib")
    pkg = self.fetch_package(pkg_name=GLIB_PKG)
    if self.check_download_only(GLIB_PKG): return

    # Mock executables.
    if (not op.isdir(op.join(self.base_dir, "bin"))) :
      os.makedirs(op.join(self.base_dir, "bin"))
    msgfmt_bin = op.join(self.base_dir, "bin", "msgfmt")
    gettext_bin = op.join(self.base_dir, "bin", "xgettext")
    self.touch_file(msgfmt_bin)
    self.touch_file(gettext_bin)
    self.call("chmod 744 \"%s\""%msgfmt_bin)
    self.call("chmod 744 \"%s\""%gettext_bin)
    # glib
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    self.configure_and_build(
      config_args=[self.prefix],
      log=pkg_log)

  def build_expat(self):
    # expat
    pkg_log = self.start_building_package("expat")
    pkg = self.fetch_package(pkg_name=EXPAT_PKG)
    if self.check_download_only(EXPAT_PKG): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    self.configure_and_build(config_args=[self.prefix], log=pkg_log)
    header_files = ["./lib/expat_external.h", "./lib/expat.h"]
    for header in header_files :
      self.call("./conftools/install-sh -c -m 644 %s \"%s\"" % (header,
        op.join(self.base_dir, "include")), log=pkg_log)

  def build_fontconfig(self):
    # fontconfig
    pkg_log = self.start_building_package("fontconfig")
    pkg = self.fetch_package(pkg_name=FONTCONFIG_PKG)
    if self.check_download_only(FONTCONFIG_PKG): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    # Create font directories.
    if (not op.isdir(op.join(self.base_dir, "share", "fonts"))) :
      os.makedirs(op.join(self.base_dir, "share", "fonts"))
    if (not op.isdir(op.join(self.base_dir,"etc","fonts"))) :
      os.makedirs(op.join(self.base_dir,"etc","fonts"))
    fc_config_args = [ self.prefix,
      "--disable-docs",
      "--with-expat-includes=\"%s\"" % op.join(self.base_dir, "include"),
      "--with-expat-lib=\"%s\"" % op.join(self.base_dir, "lib"),
      "--with-add-fonts=\"%s\"" % op.join(self.base_dir,"share","fonts"),
      "--with-confdir=\"%s\"" % op.join(self.base_dir,"etc","fonts"),
      "--with-docdir=\"%s\"" % op.join(self.base_dir, "doc"),
      "--with-freetype-config=freetype-config", ]
    self.configure_and_build(config_args=fc_config_args, log=pkg_log)

  def build_render(self):
    # render, xrender, xft
    for pkg, name in zip([RENDER_PKG,XRENDER_PKG,XFT_PKG], ["render","Xrender","Xft"]):
      self.build_compiled_package_simple(pkg_name=pkg, pkg_name_label=name)

  def build_pixman(self):
    # pixman
    pkg_log = self.start_building_package("pixman")
    pkg = self.fetch_package(pkg_name=PIXMAN_PKG)
    if self.check_download_only(PIXMAN_PKG): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    self.configure_and_build(
      config_args=[self.prefix, "--disable-gtk" ],
      log=pkg_log)

  def build_cairo(self):
    # cairo, pango, atk
    for pkg, name in zip([CAIRO_PKG, PANGO_PKG, ATK_PKG], ["cairo", "pango", "atk"]) :
      self.build_compiled_package_simple(pkg_name=pkg, pkg_name_label=name)

  def build_tiff(self):
    # tiff
    pkg_log = self.start_building_package("tiff")
    pkg = self.fetch_package(pkg_name=TIFF_PKG)
    if self.check_download_only(TIFF_PKG): return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    os.environ['MANSCHEME'] = "bsd-source-cat"
    os.environ['DIR_MAN'] = op.join(self.base_dir, "man")
    config_args = [self.prefix, "--noninteractive", "--with-LIBGL=no", "--with-LIBIMAGE=no" ]
    self.configure_and_build(
      config_args=config_args,
      log=pkg_log)

  def build_gtk(self):
    # gtk+
    pkg_log = self.start_building_package("gtk+")
    pkg, size = self.fetch_package(pkg_name=GTK_PKG, return_file_and_status=True)
    self.fetch_package(pkg_name=GTK_ENGINE_PKG)
    if self.check_download_only(GTK_PKG) and self.check_download_only(GTK_ENGINE_PKG):
      return
    self.untar_and_chdir(pkg=pkg, log=pkg_log)
    gtk_config_args = [
      self.prefix,
      "--disable-cups",
      "--without-libjpeg",
    ]

    # cups breaks the install; --disable-cups not respected; kludge by renaming
    # the cups thing (why are we including printer drivers here anyway?)

    import shutil
    shutil.copyfile('configure', 'configure.save')

    configure = open('configure.save', 'r').read().replace(
      'cups-config', 'junk-config')
    open('configure', 'w').write(configure)
    os.system('chmod +x ./configure')

    # end kludge

    self.call("./configure %s" % " ".join(gtk_config_args), log=pkg_log)
    self.call("make -j %d SRC_SUBDIRS='gdk-pixbuf gdk gtk modules'" %
      self.nproc, log=pkg_log)
    self.call("make install SRC_SUBDIRS='gdk-pixbuf gdk gtk modules'",
      log=pkg_log)

    # not strictly necessary, but prints 'installing' message
    if self.check_download_only(GTK_ENGINE_PKG): return
    # gtk-engine
    self.build_compiled_package_simple(pkg_name=GTK_ENGINE_PKG,
      pkg_name_label="gtk-engine")

  def build_fonts(self):
    # fonts
    fonts_log = self.start_building_package("fonts")
    pkg = self.fetch_package(pkg_name=FONT_PKG)
    if self.check_download_only("fonts"): return

    share_dir = op.join(self.base_dir, "share")
    if (not op.isdir(share_dir)) :
      os.makedirs(share_dir)
    os.chdir(share_dir)
    untar(pkg, log=fonts_log, verbose=True)
    os.chdir(self.tmp_dir)

  def build_wxpython (self) :
    pkg_log = self.start_building_package("wxPython")
    pkg_name = WXPYTHON_PKG
    # XXX we don't entirely trust wxPython-2.9, but it would be preferrable for
    # the future to use a single version instead
    if (self.flag_is_mac or self.options.use_wxpython3) :
      pkg_name = WXPYTHON_DEV_PKG
    pkg = self.fetch_package(pkg_name)
    if self.check_download_only(pkg_name): return

    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)

    # Unconditionally append Debian i386/x86_64 multilib directories
    # to wxPython's list of library search paths.
    line = "SEARCH_LIB=\"`echo \"$SEARCH_INCLUDE\" | " \
           "sed s@include@$wx_cv_std_libpath@g` /usr/$wx_cv_std_libpath"
    self.patch_src(src_file="configure",
                   target=(line, ),
                   replace_with=(line +
                                 " /usr/lib/i386-linux-gnu" +
                                 " /usr/lib/x86_64-linux-gnu", ))

    if (self.flag_is_mac and get_os_version() == "10.10") :
      # Workaround wxwidgets 3.0.2 compilation error on Yosemite
      # This will be fixed in 3.0.3.
      # See:
      #   http://trac.wxwidgets.org/ticket/16329
      #   http://goharsha.com/blog/compiling-wxwidgets-3-0-2-mac-os-x-yosemite/
      print >> self.log, "  patching src/osx/webview_webkit.mm"
      self.patch_src(src_file="src/osx/webview_webkit.mm",
                     target=("#include <WebKit/WebKit.h>",),
                     replace_with=("#include <WebKit/WebKitLegacy.h>",))

    # Stage 1: build wxWidgets libraries
    config_opts = [
      self.prefix,
      "--with-opengl",
      "--enable-unicode"
    ]
    if (self.options.debug) :
      config_opts.extend(["--disable-optimize",
                          "--enable-debug"])
      if (self.flag_is_linux) :
        config_opts.append("--disable-debug_gdb")
    else :
      config_opts.extend(["--enable-optimize",
                          "--disable-debugreport"])

    # if (cocoa) :
    if (self.flag_is_mac) :
      config_opts.extend([
        "--with-osx_cocoa",
        #"--with-macosx-version-min=%s" % get_os_version(),
        "--with-mac",
        "--enable-monolithic",
        "--disable-mediactrl",
        "--without-liblzma"
      ])
    elif (self.flag_is_linux) :
      config_opts.extend([
        "--with-gtk",
        "--with-gtk-prefix=\"%s\"" % self.base_dir,
        "--with-gtk-exec-prefix=\"%s\"" % op.join(self.base_dir, "lib"),
        "--enable-graphics_ctx",
        "--enable-mediactrl",
        "--enable-display",
      ])

    install_gizmos = True
    print >> self.log, "  building wxWidgets with options:"
    for opt in config_opts :
      print >> self.log, "    %s" % opt
    self.call("./configure %s" % " ".join(config_opts), log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    if (not self.flag_is_mac) : # XXX ???
      self.call("make -j %d -C contrib/src/stc" % self.nproc, log=pkg_log)
      if install_gizmos:
        self.call("make -j %d -C contrib/src/gizmos" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)
    if (not self.flag_is_mac) : # XXX ???
      self.call("make -C contrib/src/stc install", log=pkg_log)
      if install_gizmos:
        self.call("make -C contrib/src/gizmos install", log=pkg_log)

    # Stage 2: build wxPython itself
    wxpy_build_opts = [
      "BUILD_GLCANVAS=1",
      "BUILD_DLLWIDGET=0",
      "UNICODE=1"
    ]
    if install_gizmos:
      wxpy_build_opts.append("BUILD_GIZMOS=1")
    if self.flag_is_mac:
      os.environ['CFLAGS'] = os.environ.get('CFLAGS', '') + " -arch x86_64"
      wxpy_build_opts.extend(["BUILD_STC=1",
                              "WXPORT=osx_cocoa"])
    else :
      wxpy_build_opts.extend(["BUILD_STC=0",
                              "BUILD_OGL=0",
                              "WX_CONFIG=%s/bin/wx-config" %self.base_dir])
    self.chdir("wxPython", log=pkg_log)
    debug_flag = ""
    if (self.options.debug) :
      debug_flag = "--debug"
    print >> self.log, "  building wxPython with options:"
    for opt in wxpy_build_opts :
      print >> self.log, "    %s" % opt
    self.call("%s setup.py %s build_ext %s" % (self.python_exe,
      " ".join(wxpy_build_opts), debug_flag), log=pkg_log)
    self.call("%s setup.py %s install" % (self.python_exe,
      " ".join(wxpy_build_opts)), log=pkg_log)
    self.verify_python_module("wxPython", "wx")

  def build_wxpython3(self):
    assert 0
    pkg_log = self.start_building_package("wxPython")
    pkg_name = WXPYTHON_DEV_PKG
    pkg = self.fetch_package(pkg_name)
    if self.check_download_only(pkg_name):
      return

    pkg_dir = untar(pkg, log=pkg_log)
    os.chdir(pkg_dir)
    if (self.flag_is_mac and get_os_version() == "10.10") :
      # Workaround wxwidgets 3.0.2 compilation error on Yosemite
      print >> self.log, "  patching src/osx/webview_webkit.mm"
      self.patch_src(src_file="src/osx/webview_webkit.mm",
                     target=("#include <WebKit/WebKit.h>",),
                     replace_with=("#include <WebKit/WebKitLegacy.h>",))

    # Stage 1: build wxWidgets libraries
    config_opts = [
      self.prefix,
      "--disable-mediactrl",
      "--with-opengl",
    ]
    if (self.flag_is_mac) :
      config_opts.extend([
        "--with-mac",
        "--with-osx_cocoa",
        "--with-macosx-version-min=10.6",
        "--enable-monolithic",
        "--enable-unicode"
        "--enable-monolithic"
      ])
    elif (self.flag_is_linux) :
      config_opts.extend([
        "--with-gtk",
        "--with-gtk-prefix=\"%s\"" % self.base_dir,
        "--with-gtk-exec-prefix=\"%s\"" % op.join(self.base_dir, "lib"),
        "--enable-graphics_ctx",
      ])

    print >> self.log, "  configure wxWidgets with options:"
    for opt in config_opts :
      print >> self.log, "    %s" % opt
    self.call("./configure %s" % " ".join(config_opts), log=pkg_log)
    self.call("make -j %d" % self.nproc, log=pkg_log)
    self.call("make install", log=pkg_log)

    os.chdir(os.path.join(pkg_dir, 'wxPython'))
    self.call([
      self.python_exe,
      'build-wxpython.py',
      '--no_wxbuild',
      '--osx_cocoa',
      '--prefix', self.base_dir,
      '--install'
      ], log=pkg_log)
    self.verify_python_module("wxPython", "wx")


  def build_matplotlib(self):
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

  def build_misc (self) :
    self.build_python_module_simple(
      pkg_url=BASE_CCI_PKG_URL,
      pkg_name=PYRTF_PKG,
      pkg_name_label="PyRTF",
      confirm_import_module="PyRTF")
    # TODO we are patching the source to force it to use the correct backend.
    # I'm not sure if this is truly necessary or if there's a cleaner way...
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

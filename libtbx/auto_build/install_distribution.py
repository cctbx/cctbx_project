#!/usr/bin/python

"""
Template for building Unix/Linux installers for CCTBX-based packages, e.g.
Phenix or LABELIT.  This is derived from the Phenix installer, with options
specific to Phenix removed.

The installer makes several assumptions about the distribution bundle, and
will break if it does not encounter the expected files and directories.  The
recommended way to create an installer package is to start with the script
setup_installer.py in the same directory.  The expected directory structure is:

INSTALLER_DIR/bin/
INSTALLER_DIR/bin/install   # main executable, subclasses 'installer' class
INSTALLER_DIR/lib
INSTALLER_DIR/lib/libtbx/   # complete libtbx distribution
INSTALLER_DIR/bundles/      # compressed binary distribution
INSTALLER_DIR/source/       # source tarballs
INSTALLER_DIR/dependencies  # third-party software

The last four directories contain the actual distribution; depending on whether
the installer is source or binary, they may not all be present.  The 'bundles'
directory will have one or more tar.gz files containing all compiled files
(and many Python modules).  The 'sources' directory will contain a set of
tar.gz files for the distributed modules, one of which will be a CCTBX source
bundle (plus various add-ons like ccp4io, cbflib, etc.).  'dependencies' is for
third-party packages such as Python that need to be compiled for source
installation.

In addition to the required contents, the installer may also contain files
named README, LICENSE, etc.; these will automatically be propagated to the
installed package.  Information describing the version may also be provided
but the implementation of this is left up to each product.  Installers may also
have additional directories for non-compiled packages such as geometry
restraints or documentation.

To implement an installer for distribution, an install script should be written
along these lines:

  from libtbx.auto_build import install_distribution
  import os.path
  import sys

  class installer (install_distribution.installer) :
    product_name = "PHENIX"
    dest_dir_prefix = "phenix"
    make_apps = ["phenix"]
    installer_dir = os.path.dirname(os.path.dirname(__file__))

  if (__name__ == "__main__") :
    installer(sys.argv[1:])

Note that the location of the install script is not mandatory, but it must be
able to correctly detect the base installer directory.
"""

from __future__ import division
from optparse import OptionParser
import os.path as op
import shutil
import time
import os
import sys
# set nproc automatically if possible
try :
  import multiprocessing
  NPROC_DEFAULT = min(24, multiprocessing.cpu_count())
except Exception :
  NPROC_DEFAULT = 1
# local imports
# XXX HACK
libtbx_path = op.abspath(op.dirname(op.dirname(__file__)))
if (not libtbx_path in sys.path) :
  sys.path.append(libtbx_path)
from libtbx.auto_build import write_gui_dispatcher_include
from libtbx.auto_build import regenerate_module_files
from libtbx.auto_build import install_base_packages
from libtbx.auto_build import create_mac_app
from libtbx.auto_build import package_defs
from libtbx.auto_build.installer_utils import *

class InstallerError (Exception) :
  pass

class installer (object) :
  """
  Base class for installers.  Some methods and class attributes must be
  re-implemented in subclasses!
  """
  # some basic configuration variables - override in subclasses
  product_name = "CCTBX"
  destination = "/usr/local"
  dest_dir_prefix = "cctbx"
  installer_dir = os.environ.get(product_name + "_INSTALLER", None)
  include_gui_packages = True
  remove_sources_default = False
  # options passed to install_base_packages.py
  base_package_options = []
  source_packages = [
    "cctbx_bundle",
  ]
  # modules that need to be configured for use in the final installation.
  # this will automatically include dependencies.
  configure_modules = ["mmtbx", "smtbx"]
  # programs to make graphical .app launchers for (Mac only)
  make_apps = []
  # architectures supported for a particular distribution.  those listed here
  # will be detected automatically.
  supported_mtypes = [
    "intel-linux-2.6",
    "intel-linux-2.6-x86_64",
    "mac-intel-osx",
    "mac-intel-osx-x86_64",
  ]

  def __init__ (self, args, out=sys.stdout) :
    self.args = args
    self.out = out
    self.apps_built = False # flag for Mac OS
    check_python_version()
    self.parse_options()
    self.basic_setup()
    self.determine_installation_type()
    self.setup_directories()
    if (self.binary_install) :
      self.install_from_binary()
    else :
      self.install_from_source()
    self.install_finalize()

  def parse_options (self) :
    """
    Process command-line options.  These may be supplemented by subclasses that
    override the method add_product_specific_options.
    """
    parser = OptionParser(
      description=("Command-line installer for %s crystallography "+
        "software on Mac and Linux systems.  Please note that Mac users "+
        "may also use the graphical installer, but this is limited to single "+
        "workstations.") % self.product_name)
    parser.add_option("--prefix", dest="prefix", action="store",
      help="Directory where %s is to be installed" % self.product_name,
      default=self.destination)
    parser.add_option("--nproc", dest="nproc", action="store", type="int",
      help="Number of processors for source install", default=NPROC_DEFAULT)
    if (self.include_gui_packages) :
      parser.add_option("--no-gui", dest="no_gui", action="store_true",
        help="Disable building of GUI dependencies (source install only)",
        default=False)
    parser.add_option("--base-only", dest="base_only", action="store_true",
      help="Only install base libraries and files", default=False)
    parser.add_option("--binary", dest="binary", action="store_true",
      help="Use pre-compiled binary install", default=None)
    parser.add_option("--source", dest="source", action="store_true",
      help="Force source installation", default=None)
    parser.add_option("--openmp", dest="openmp", action="store_true",
      help="Enable OpenMP compilation if possible", default=False)
    parser.add_option("--no-opt", dest="no_opt", action="store_true",
      help="Disable optimization during compilation", default=False)
    parser.add_option("--no-app", dest="no_app", action="store_true",
      help="Disable generation of Mac app launcher(s)", default=False)
    parser.add_option("--debug", dest="debug", action="store_true",
      help="Turn on debugging during compilation", default=False)
    parser.add_option("--top-level-sources", action="store_true",
      dest="top_level_sources",
      help="Keep modules at top-level directory instead of 'modules' "+
           "subdirectory (source installer only)", default=False)
    if (self.remove_sources_default) :
      parser.add_option("--no-compact", dest="compact", action="store_false",
        help="Don't remove unnecessary files such as compiled sources")
    else :
      parser.add_option("--compact", dest="compact", action="store_true",
        help="Remove unnecessary files such as compiled sources to save space")
    parser.add_option("--try-unsupported", dest="try_unsupported",
      action="store_true", default=False,
      help="Attempt source install on unsupported platform")
    parser.add_option("--python_static", default=False, action="store_true",
      help="Compile Python as static executable and library (Linux only)")
    parser.add_option("--makedirs", default=False, action="store_true",
      help="Create installation path prefix if not already present")
    self.add_product_specific_options(parser)
    self.options, args = parser.parse_args(self.args)

  def machine_type (self) :
    """
    Determine the mtype string.  The four pre-defined mtypes will be detected
    automatically (with Linux kernel 3.x treated as 2.6), but additional mtypes
    can be used if this method is re-implemented.
    """
    return machine_type()

  def basic_setup (self) :
    out = self.out
    self.mtype = self.machine_type()
    self.version = self.get_version()
    assert (self.version is not None)
    self.src_dir = op.join(self.installer_dir, "source")
    self.bundle_dir = op.join(self.installer_dir, "bundles")
    self.dependencies_dir = op.join(self.installer_dir, "dependencies")
    print >> out, """
  ==========================================================================
                        %(product)s Installation

                        version: %(version)s
                   machine type: %(mtype)s
                     OS version: %(os)s
                     user shell: %(shell)s
                    destination: %(dest)s
                # of processors: %(nproc)s
  =========================================================================
  """ % { "product" : self.product_name,
          "version" : self.version,
          "mtype"   : self.mtype,
          "os"      : get_os_version(),
          "shell"   : os.environ['SHELL'],
          "dest"    : self.options.prefix,
          "nproc"   : self.options.nproc, }
    if (self.mtype is None) :
      raise InstallerError("Machine type not recognized")
    elif ((not self.mtype in self.supported_mtypes) and
          (not self.options.try_unsupported)) :
      raise InstallerError("""
  %(mtype)s is not a supported platform, installation aborted
    use the --try-unsupported option to attempt installation on this platform
    use the --no-gui option for a core %(product)s installation
  """ % { "mtype" : self.mtype, "product" : self.product_name})
    if (not op.isdir(self.options.prefix)) :
      if (self.options.makedirs) :
        pass
      else :
        raise InstallerError(
          "Destination directory does not exist:\n" +
          ("  %s" % self.options.prefix) +
          "If you are sure you want to install here, re-run with the extra\n"+
          "argument --makedirs.")
    self.flag_build_gui = False
    if (self.include_gui_packages) :
      self.flag_build_gui = not self.options.no_gui
    # environment variables required by other scripts
    os.environ["%s_MTYPE" % self.product_name] = self.mtype
    os.environ["%s_INSTALLER" % self.product_name] = self.installer_dir

  def determine_installation_type (self) :
    """
    Inspect the installer contents to determine whether this is a binary
    distribution or source packages.  If both are present, binary is preferred.
    """
    out = self.out
    source_install = binary_install = None
    if (op.isdir(self.src_dir)) :
      source_install = True
    bundle_dir = op.join(self.installer_dir, "bundles")
    self.base_bundle = op.join(self.bundle_dir, "base-%s-%s.tar.gz" %
      (self.version, self.mtype))
    self.build_bundle_file = op.join(self.bundle_dir, "build-%s-%s.tar.gz" %
      (self.version, self.mtype))
    self.modules_bundle_file = op.join(self.bundle_dir, "modules-%s-%s.tar.gz" %
      (self.version, self.mtype))
    self.alias_mtype = self.mtype
    if (op.isfile(self.build_bundle_file)) :
      binary_install = True
    elif (self.mtype.endswith("x86_64")) :
      mtype_32bit = re.sub("-x86_64", "", self.mtype)
      bundle32 = op.join(bundle_dir, "build-%s-%s.tar.gz" % (self.version,
        mtype_32bit))
      if (op.isfile(bundle32)) :
        print >> out, "This is the 32-bit build for your platform- "
        print >> out, "automatically switching machine type to %s." % \
          self.alias_mtype
        self.alias_mtype = mtype_32bit
        binary_install = True
        self.build_bundle_file = bundle32
        self.modules_bundle_file = op.join(self.bundle_dir,
          "modules-%s-%s.tar.gz" % (self.version, self.alias_mtype))
        self.base_bundle = op.join(self.bundle_dir, "base-%s-%s.tar.gz" %
          (self.version, self.alias_mtype))
    if (not binary_install) :
      print >> out, "No binary bundles found for %s" % self.mtype
      if (op.isdir(bundle_dir)) :
        print >> out, "bundles directory contents:"
        for file_name in os.listdir(bundle_dir) :
          print >> out,"  %s" % file_name
        print >> out, "expected:"
        print >> out, op.basename(self.build_bundle_file)
      elif (source_install) :
        print "Okay, this must be the source-only installer."
    if (not binary_install) and (not source_install) :
      raise InstallerError("no binary installer found for %s, aborting." %
        self.mtype)
    elif (source_install and binary_install) :
      if (self.options.source) :
        print >> out, "performing source installation"
        binary_install = False
      else :
        print >> out, "source and binary both available, defaulting to "+\
                      "binary install\n"
        source_install = False
    # check for wxPython compilation dependencies on Linux
    if ((source_install) and (self.flag_build_gui) and
        (sys.platform != "darwin")) :
      tmp_log = open("/dev/null", "w")
      error = install_base_packages.check_wxpython_build_dependencies(
        log=tmp_log)
      if (error is not None) :
        raise InstallerError(error)
    self.source_install = source_install
    self.binary_install = binary_install
    self.base_bundle = op.join(self.bundle_dir, "base-%s-%s.tar.gz" %
      (self.version, self.alias_mtype))
    self.base_binary_install = op.isfile(self.base_bundle)
    return True

  def setup_directories (self) :
    out = self.out
    hostname = os.uname()[1].split(".")[0]
    self.tmp_dir = op.join(self.installer_dir, "tmp", self.mtype, hostname)
    if (not op.isdir(self.tmp_dir)) :
      os.makedirs(self.tmp_dir)
    self.dest_dir = op.join(self.options.prefix, "%s-%s" %
      (self.dest_dir_prefix, self.version))
    if (not op.isdir(self.dest_dir)) :
      try :
        os.makedirs(self.dest_dir)
      except OSError, e :
        if (e.errno == 13) :
          raise InstallerError("you do not have write permissions to %s" %
            self.options.prefix)
        else :
          raise InstallerError("can't create '%s': %s" % (self.dest_dir, e))
    else :
      print >> out, "Warning: %s already exists, will be overwritten" % \
        self.dest_dir
      if (not os.access(self.dest_dir, os.W_OK)) :
        raise InstallerError("you do not have write permissions to %s" %
          self.dest_dir)
    self.build_dir = op.join(self.dest_dir, "build")
    self.base_dir = op.join(self.dest_dir, "base")
    if (not op.exists(self.build_dir)) :
      os.makedirs(self.build_dir)
    if (not self.options.top_level_sources) :
      self.modules_dir = op.join(self.dest_dir, "modules")
      if (not op.exists(self.modules_dir)) :
        os.makedirs(self.modules_dir)
    else :
      self.modules_dir = op.join(self.dest_dir)
    # environment variables required by other scripts
    os.environ["%s_LOC" % self.product_name] = self.dest_dir
    os.environ["%s_BUILD" % self.product_name] = self.build_dir

  def show_installation_paths (self) :
    out = self.out
    print >> out, """
%(product)s installation target directory <%(product)s_LOC> set to:
   %(dest_dir)s
%(product)s installation source directory set to:
   %(inst_dir)s
%(product)s installation build directory set to:
   %(build_dir)s
%(product)s temporary build directory set to:
   %(tmp_dir)s
""" % { "product" : self.product_name,
        "dest_dir" : self.dest_dir,
        "build_dir" : self.build_dir,
        "inst_dir" : self.installer_dir,
        "tmp_dir" : self.tmp_dir, }

  #---------------------------------------------------------------------
  # SOURCE INSTALL
  #
  def install_from_source (self) :
    out = self.out
    self.show_installation_paths()
    # PART 1: dependencies
    # install Python 2.7 and other dependencies (wxPython, etc.)
    if (not op.exists(self.base_dir)) :
      os.makedirs(self.base_dir)
    if (self.base_binary_install) :
      os.chdir(self.dest_dir)
      print >> out, "Using precompiled base bundle..."
      untar(self.base_bundle, log=out, verbose=True,
        change_ownership=False, check_output_path=False)
      os.chdir(self.installer_dir)
    else :
      base_args = [
        # XXX will this allow spaces in path names?
        "--build_dir=%s" % self.dest_dir,
        "--tmp_dir=%s" % self.tmp_dir,
        "--pkg_dir=%s" % self.src_dir,
        "--pkg_dir=%s" % self.dependencies_dir,
        "--nproc=%s" % self.options.nproc,
        "--no-download",
      ] + self.base_package_options
      # XXX wxPython (et al.) installation is optional depending on the
      # tarfile actually being present
      wxpython_file = op.join(self.dependencies_dir, package_defs.WXPYTHON_PKG)
      if (self.flag_build_gui) and op.isfile(wxpython_file) :
        base_args.append("--gui")
      elif ("--all" in base_args) :
        base_args.remove("--all")
      if (self.options.debug) :
        base_args.append("--debug")
      if (not self.options.python_static) :
        base_args.append("--python-shared")
      install_base_packages.installer(
        args=base_args,
        log=out)
    python_bin = op.join(self.base_dir, "bin", "python")
    assert op.isfile(python_bin)
    print >> out, ""
    print >> out, "******************************"
    print >> out, "* BASE INSTALLATION COMPLETE *"
    print >> out, "******************************"
    print >> out, ""
    print >> out, "-" * 80
    print >> out, ""
    # PART 2: product packages
    print >> out, "Building core %s components" % self.product_name
    if (self.options.base_only) :
      print >> out, "--base-only was specified, so skipping %s libraries" % \
        self.product_name
      return
    for pkg_name in self.source_packages :
      pkg_file = op.join(self.src_dir, pkg_name + ".tar.gz")
      if (not op.isfile(pkg_file)) :
        raise RuntimeError("Can't find %s" % pkg_file)
      print >> out, "  unwrapping %s..." % op.basename(pkg_file)
      log_file = op.join(self.tmp_dir, "%s.log" % pkg_name)
      log = open(log_file, "w")
      os.chdir(self.modules_dir)
      untar(pkg_file, log=log, verbose=True, change_ownership=False,
        check_output_path=(not "cctbx_bundle" in pkg_file))
      log.close()
      if (pkg_file.startswith("cctbx_bundle")) :
        if (not op.isdir("cctbx_project")) :
          raise RuntimeError("Directory 'cctbx_project' not found!")
    tag_file = op.join(self.modules_dir, "TAG")
    if op.isfile(tag_file) :
      os.rename(tag_file, op.join(self.modules_dir, "cctbx_bundle_TAG"))
    self.product_specific_setup_before_compile(log=out)
    log_file = op.join(self.tmp_dir, "compile_%s.log" % self.dest_dir_prefix)
    log = open(log_file, "w")
    os.chdir(self.build_dir)
    if op.isfile("libtbx_refresh_is_completed") :
      os.remove("libtbx_refresh_is_completed")
    config_path = op.join(self.modules_dir, "cctbx_project", "libtbx",
      "configure.py")
    config_args = [
      "--scan_boost",
      "--command_version_suffix=%s" % self.version,
      "--current_working_directory=\"%s\"" % self.build_dir,
    ] + self.configure_modules
    if self.options.openmp :
      config_args.append("--enable-openmp-if-possible=True")
    if (self.options.debug) :
      print >> out, ("configuring %s components with debugging..." % \
        self.product_name),
      config_args.insert(0, "--build=debug")
    elif (self.options.no_opt) :
      print >> out, ("configuring %s components with no optimization..." % \
        self.product_name),
      config_args.insert(0, "--build=quick")
    else :
      print >> out, ("configuring %s components with optimization..." % \
        self.product_name),
      config_args.insert(0, "--build=release")
    call(" ".join([python_bin, config_path] + config_args), log=log)
    print >> out, "ok"
    # XXX in the original shell script, I had to run 'libtbx.configure reel'
    # to get REEL actually configured - need to double-check this
    scons_args = ["."]
    if (not os.environ.get("LIBTBX_FULL_TESTING") in [None, ""]) :
      scons_args.append("-k")
    scons_args.extend(["-j", str(self.options.nproc)])
    print >> out, ("compiling %s components - this might take an hour..." % \
      self.product_name),
    scons_bin = op.join(self.build_dir, "bin", "libtbx.scons")
    assert op.isfile(scons_bin)
    call(" ".join([scons_bin] + scons_args), log=log)
    print >> out, "done"
    self.product_specific_source_install(log=log)
    print >> out, ""
    msg = "%s INSTALLATION COMPLETE" % self.product_name
    print >> out, "*" * (len(msg) + 4)
    print >> out, "* " + msg + " *"
    print >> out, "*" * (len(msg) + 4)
    print >> out, ""

  #---------------------------------------------------------------------
  # BINARY INSTALL
  #
  def install_from_binary (self) :
    """
    Unpackage the binary bundle in the destination directory.
    """
    out = self.out
    tmp_bin_dir = op.join(self.tmp_dir, "binary")
    log_file = op.join(self.tmp_dir, "binary.log")
    log = open(log_file, "w")
    if (op.exists(tmp_bin_dir)) :
      print >> out, "removing existing temporary directory...",
      shutil.rmtree(tmp_bin_dir)
      print >> out, "ok"
    os.makedirs(tmp_bin_dir)
    # check for existing installations
    print >> out, "finding existing installations..."
    lib_dir = op.join(self.build_dir, "lib")
    if (op.exists(lib_dir)) :
      print >> out, "Removing out-of-date lib directory...",
      shutil.rmtree(lib_dir)
      print >> out, "ok"
    if (op.exists(self.base_dir)) :
      print >> out, "Removing out-of-date base software directory...",
      shutil.rmtree(self.base_dir)
      print >> out, "ok"
    # XXX there used to be a bunch of version check stuff in here, but this
    # is no longer necessary
    # unwrap the tarball
    print >> out, "Installing new binary package...",
    os.chdir(self.dest_dir)
    if (self.base_binary_install) :
      untar(self.base_bundle, log=log, verbose=True,
        change_ownership=False, check_output_path=False)
    untar(self.build_bundle_file, log=log, verbose=True,
      change_ownership=False,
      check_output_path=False)
    os.chdir(self.modules_dir)
    untar(self.modules_bundle_file, log=log, verbose=True,
      change_ownership=False,
      check_output_path=False)
    if (not op.isdir(lib_dir)) or (not op.isdir(self.base_dir)) :
      raise RuntimeError("One or more missing directories:\n  %s\n  %s" %
        (lib_dir, self.base_dir))
    if (sys.platform != "darwin") :
      os.environ["LD_LIBRARY_PATH"] = op.join(self.base_dir, "lib")
    print >> out, "ok"
    self.product_specific_binary_install(log=log)
    os.chdir(self.build_dir)
    print >> out, "Configuring %s components..." % (self.product_name),
    self.reconfigure(log=log)
    print >> out, "ok"
    print >> out, ""
    msg = "%s BINARY INSTALLATION COMPLETE" % self.product_name
    print >> out, "*" * (len(msg) + 4)
    print >> out, "* " + msg + " *"
    print >> out, "*" * (len(msg) + 4)
    print >> out, ""

  #---------------------------------------------------------------------
  # NON-COMPILED COMPONENTS AND FINAL SETUP
  #
  def reconfigure (self, log) :
    """
    Run libtbx/configure.py to configure the build in the new location.
    """
    os.chdir(self.build_dir)
    args = [
      "\"%s/bin/python\"" % self.base_dir,
      "\"%s/cctbx_project/libtbx/configure.py\"" % self.modules_dir,
      "--current_working_directory=%s" % self.build_dir,
      "--command_version_suffix=%s" % self.version,
    ] + self.configure_modules
    try :
      call(args=" ".join(args), log=log)
    except RuntimeError :
      raise InstallerError("configuration step incomplete!  See the log file "+
        "for detailed error messages.")

  def install_finalize (self) :
    """
    Set up dispatchers and assorted shell scripts, create app bundles, etc.
    """
    out = self.out
    print >> out, ""
    print >> out, "*" * 72
    print >> out, "FINALIZING %s INSTALLATION" % self.product_name
    log_path = op.join(self.tmp_dir, "install_finalize.log")
    print >> out, "  log file is %s" % log_path
    log = open(log_path, "w")
    if (self.flag_build_gui) and (sys.platform != "darwin") :
      os.environ["LD_LIBRARY_PATH"] = op.join(self.base_dir, "lib")
      regenerate_module_files.run(args=["--build_dir=%s" % self.dest_dir],
        out=out)
    # write dispatcher_include file
    print >> out, "generating %s environment additions for dispatchers" % \
      self.product_name
    print >> out, ""
    dispatcher = op.join(self.build_dir, "dispatcher_include_%s.sh" %
      self.dest_dir_prefix)
    if (op.isfile(dispatcher)) :
      os.remove(dispatcher)
    env_prefix = self.product_name.upper() # e.g. "Phenix" -> "PHENIX"
    prologue = "\n".join([
      "export %s=\"%s\"" % (env_prefix, self.dest_dir),
      "export %s_VERSION=%s" % (env_prefix, self.version),
      "export %s_ENVIRONMENT=1" % env_prefix,
      "export %s_MTYPE=%s" % (env_prefix, self.mtype),
    ] + self.product_specific_dispatcher_prologue())
    epilogue = "\n".join(self.product_specific_dispatcher_epilogue())
    dispatcher_opts = [
      "--build_dir=%s" % self.build_dir,
      "--base_dir=%s" % self.base_dir,
      "--suffix=%s" % self.dest_dir_prefix,
      "--gtk_version=2.10.0", # XXX this can change!
      "--quiet",
    ]
    if (not self.flag_build_gui) :
      dispatcher_opts.append("--ignore_missing_dirs")
    # FIXME this will happen regardless of whether the GUI modules are being
    # distributed or not - will this be problematic?
    write_gui_dispatcher_include.run(
      args=dispatcher_opts,
      prologue=prologue,
      epilogue=epilogue,
      out=out)
    assert op.isfile(dispatcher)
    self.write_environment_files(out)
    # run configure.py to generate dispatchers
    print >> out, "configuring %s components..." % self.product_name
    os.chdir(self.build_dir)
    if (op.exists("libtbx_refresh_is_completed")) :
      os.remove("libtbx_refresh_is_completed")
    self.reconfigure(log=log)
    os.chdir(self.build_dir)
    assert op.isfile("setpaths.sh")
    bin_dir = op.join(self.build_dir, "bin")
    os.environ["PATH"] = "%s:%s" % (bin_dir, os.environ["PATH"])
    # compile .py files
    print >> out, "precompiling .py files...",
    os.chdir(self.modules_dir)
    call(args="libtbx.py_compile_all", log=log)
    print >> out, "ok"
    # copy README et al.
    for file_name in ["CHANGES","LICENSE","README","README-DEV","SOURCES"] :
      src_file = op.join(self.installer_dir, file_name)
      if op.exists(src_file) :
        # XXX use our own implementation instead of shutil.copyfile
        copy_file(src_file, op.join(self.dest_dir, file_name))
    # generate .app (Mac only)
    if ((sys.platform == "darwin") and (len(self.make_apps) > 0) and
        (not self.options.no_app) and self.flag_build_gui) :
      os.chdir(self.build_dir)
      for app_name in self.make_apps :
        args = [
          "libtbx.create_mac_app",
          app_name,
          "--app_name=%s-%s" % (app_name, self.version),
          "--dest=%s" % self.dest_dir,
        ]
        if (self.mtype == "mac-intel-osx-x86_64") :
          args.append("--alias_build")
          #args.append("--python_interpreter=/usr/bin/python")
        print >> out, ("Generating Mac app launcher for %s..." % app_name),
        try :
          call(args=" ".join(args), log=log)
        except RuntimeError, e :
          print ""
          print "  ERROR:"
          print "  " + str(e)
          print "  installation will continue anyway."
        else :
          app_file = op.join(self.dest_dir, "%s-%s.app" %
            (app_name, self.version))
          if (not op.exists(app_file)) :
            print >> out, " failed."
            app_file = None
          else :
            print >> out, "ok"
            self.apps_built = True
    # run custom finalization
    self.product_specific_finalize_install(log)
    # remove source files if desired
    if (self.options.compact) :
      self.reduce_installation_size(out)
    # reconfigure one last time (possibly unnecessary)
    self.display_final_message()
    if self.apps_built and (not "SSH_CLIENT" in os.environ) :
      call(args=["open", self.dest_dir], log=out)

  def write_environment_files (self, out) :
    """
    Generate shell scripts in the top-level installation directory that can
    be used to set up the user environment to run the software.  This is
    implemented as a separate method because some products (e.g. Phenix) may
    have their own needs.
    """
    # csh/tcsh environment setup file
    print >> out, "generating %s environment setup scripts:" % \
      self.product_name
    csh_file = op.join(self.dest_dir, "%s_env.csh" % self.dest_dir_prefix)
    print >> out, "  csh: %s/%s_env.csh" % (self.dest_dir, self.dest_dir_prefix)
    env_prefix = self.product_name.upper() # e.g. "Phenix" -> "PHENIX"
    env_csh = open(csh_file, "w")
    env_csh.write("#!/bin/csh -f\n")
    env_csh.write("#\n")
    env_csh.write("setenv %s \"%s\"\n" % (env_prefix, self.dest_dir))
    env_csh.write("setenv %s_VERSION %s\n" % (env_prefix, self.version))
    env_csh.write("source $%s/build/setpaths.csh\n" % (env_prefix))
    env_csh.close()
    # phenix_env.sh
    sh_file = op.join(self.dest_dir, "%s_env.sh" % self.dest_dir_prefix)
    print >> out, "  sh:  %s" % sh_file
    env_sh = open(sh_file, "w")
    env_sh.write("#!/bin/sh\n")
    env_sh.write("#\n")
    env_sh.write("export %s=\"%s\"\n" % (env_prefix, self.dest_dir))
    env_sh.write("export %s_VERSION=%s\n" % (env_prefix, self.version))
    env_sh.write(". $%s/build/setpaths.sh\n" % (env_prefix))
    env_sh.close()

  def get_version (self) :
    """
    Determine the version suffix (if any) for the destination directory.  This
    must be implemented by subclasses unless a file named VERSION is present
    at the top level of the installer folder.
    """
    version_file = op.join(self.installer_dir, "VERSION")
    if op.isfile(version_file) :
      return open(version_file).read().strip()
    return NotImplementedError()

  def reduce_installation_size (self, out) :
    """
    Remove all files not required for program execution to save disk space,
    such as any C++ sources.  This can potentially save well over 100MB and is
    recommended for purely user-facing packages, but may make it more difficult
    to use packages for development purposes.
    """
    print >> out, "Removing unnecessary files to reduce disk usage...",
    # XXX should this include .o files?
    remove_extensions = [".cpp", ".cc", ".hpp", ".h", ".hh"]
    n_deleted = 0
    for dirname, dirnames, filenames in os.walk(self.modules_dir) :
      for file_name in filenames :
        for ext in remove_extensions :
          if file_name.endswith(ext) :
            full_path = op.join(dirname, file_name)
            try :
              os.remove(full_path)
            except Exception, e :
              print >> out, "  WARNING: error removing %s" % full_path
              print >> out, str(e)
            else :
              n_deleted += 1
    n_deleted_other = self.product_specific_reduce_installation_size(out)
    if (n_deleted_other is not None) :
      n_deleted += n_deleted_other
    print >> out, "%d files deleted" % n_deleted

  #---------------------------------------------------------------------
  # STUBS FOR SUBCLASSABLE METHODS
  def add_product_specific_options (self, parser) :
    """
    Add command-line options specific to the distributed package.
    """
    pass

  def product_specific_binary_install (self, log) :
    """
    Perform additional actions required for the binary installation of a
    specific product.
    """
    pass

  def product_specific_setup_before_compile (self, log) :
    """
    Perform any necessary modifications to the sources prior to compilation.
    """
    pass

  def product_specific_source_install (self, log) :
    """
    Build additional sources, e.g. Rosetta
    """
    pass

  def product_specific_dispatcher_prologue (self) :
    """
    Environment modifications to be included near the start of the dispatchers.
    """
    return []

  def product_specific_dispatcher_epilogue (self) :
    """
    Environment modifications to be included at the end of the dispatchers.
    """
    return []

  def product_specific_finalize_install (self, log) :
    """
    Additional installation setup, file cleanup, more add-ons, etc.
    """
    pass

  def product_specific_reduce_installation_size (self, log) :
    """
    Remove unused files specific to this product, and return the number deleted.
    """
    return 0

  def display_final_message (self) :
    """
    Final instructions for user, etc.
    """
    pass

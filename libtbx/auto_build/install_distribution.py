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
from libtbx.auto_build.installer_utils import *

class InstallerError(Exception):
  pass

class installer(object):
  """
  Base class for installers. Some methods and class attributes must be
  re-implemented in subclasses!
  """
  # Basic configuration variables - override in subclasses
  organization = "gov.lbl.cci"
  product_name = "CCTBX"
  destination = "/usr/local"
  dest_dir_prefix = "cctbx"
  installer_dir = os.environ.get(product_name + "_INSTALLER", None)
  include_gui_packages = True
  # Options passed to install_base_packages.py
  base_package_options = []
  base_modules = []
  modules = [
    "cctbx_project",
  ]
  # Modules that need to be configured for use in the final installation.
  #   this will automatically include dependencies.
  configure_modules = ["mmtbx", "smtbx"]
  # Programs to make graphical .app launchers for (Mac only)
  make_apps = []
  # Architectures supported for a particular distribution. Those listed here
  #   will be detected automatically.
  supported_mtypes = [
    "intel-linux-2.6",
    "intel-linux-2.6-x86_64",
    "mac-intel-osx",
    "mac-intel-osx-x86_64",
    "Windows32bit",
    "Windows64bit",
  ]

  def __init__ (self, args=None, out=sys.stdout) :
    self.args = args or []
    self.out = out
    self.parse_options()

  def parse_options (self):
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
    parser.add_option("--no-app", dest="no_app", action="store_true",
      help="Disable generation of Mac app launcher(s)", default=False)
    parser.add_option("--compact", dest="compact", action="store_true",
      help="Remove unnecessary files such as compiled sources to save space")
    parser.add_option("--try-unsupported", dest="try_unsupported",
      action="store_true", default=False,
      help="Attempt source install on unsupported platform")
    # Source only options
    parser.add_option("--no-gui", dest="no_gui", action="store_true",
      help="Disable building of GUI dependencies (source install only)",
      default=False)
    parser.add_option("--base-only", dest="base_only", action="store_true",
      help="Only install base libraries and files", default=False)
    parser.add_option("--source", dest="source", action="store_true",
      help="Force source installation", default=None)
    parser.add_option("--openmp", dest="openmp", action="store_true",
      help="Enable OpenMP compilation if possible", default=False)
    parser.add_option("--no-opt", dest="no_opt", action="store_true",
      help="Disable optimization during compilation", default=False)
    parser.add_option("--debug", dest="debug", action="store_true",
      help="Turn on debugging during compilation", default=False)
    parser.add_option("--python_static", default=False, action="store_true",
      help="Compile Python as static executable and library (Linux only)")
    # Deprecated
    parser.add_option("--makedirs", default=False, action="store_true",
      help="Create installation path prefix if not already present")
    parser.add_option("--binary", dest="binary", action="store_true",
      help="Use pre-compiled binary install", default=None)
    parser.add_option("--nopycompile", dest="nopycompile", action="store_true",
      help="Do not precompile python files", default=False)
    self.add_product_specific_options(parser)
    self.options, args = parser.parse_args(self.args)

  def install(self):
    check_python_version()
    self.parse_options()
    self.basic_setup()
    self.check_directories()
    if sys.platform != "win32":
      self.print_banner()
      if not os.path.exists('build'):
        print >> self.out, "No build directory exists; trying source installation."
        self.options.source = True
      if self.options.source:
        print >> self.out, "Source installation specified."
        self.install_from_source()
      else:
        self.install_from_binary()
    self.install_finalize()
    self.print_header('Installation complete!')

  def basic_setup (self) :
    # Check version
    self.version = self.get_version()
    assert (self.version is not None)
    # GUI Flag
    self.flag_build_gui = False
    if (self.include_gui_packages) :
      self.flag_build_gui = not self.options.no_gui

    # Check this is a supported architecture.
    self.mtype = self.machine_type()
    if (self.mtype is None) :
      raise InstallerError("Machine type not recognized")
    elif ((not self.mtype in self.supported_mtypes) and
          (not self.options.try_unsupported)) :
      raise InstallerError("""
  %(mtype)s is not a supported platform, installation aborted
    use the --try-unsupported option to attempt installation on this platform
    use the --no-gui option for a core %(product)s installation
  """%{ "mtype" : self.mtype, "product" : self.product_name})

  def check_directories(self):
    if sys.platform == "win32":
      self.dest_dir = os.getcwd()
      self.tmp_dir = self.dest_dir
      self.build_dir = op.join(self.dest_dir, "build")
      self.base_dir = op.join(self.dest_dir, "base")
      self.modules_dir = op.join(self.dest_dir, "modules")
      return
    # The default behavior for nearly all program's --prefix options
    # is to create the directory, so I don't think the --makedirs option
    # is necessary.

    # Do not overwrite an existing installation.
    self.dest_dir = op.abspath(op.join(self.options.prefix, "%s-%s"%(self.dest_dir_prefix, self.version)))
    if os.path.exists(self.dest_dir):
      raise InstallerError("""
  Installation directory already exists:
    %s
  Please remove this directory and try again.
    """%self.dest_dir)

    # Other useful directories.
    self.tmp_dir = op.join(self.installer_dir, "base_tmp")
    self.build_dir = op.join(self.dest_dir, "build")
    self.base_dir = op.join(self.dest_dir, "base")
    self.modules_dir = op.join(self.dest_dir, "modules")
    for i in [self.dest_dir, self.tmp_dir]:
      if not os.path.exists(i):
        os.makedirs(i)

    # Environment variables required by other scripts
    os.environ["%s_MTYPE" % self.product_name] = self.mtype
    os.environ["%s_INSTALLER" % self.product_name] = self.installer_dir
    os.environ["%s_LOC" % self.product_name] = self.dest_dir
    os.environ["%s_BUILD" % self.product_name] = self.build_dir

  def print_banner(self):
    print >> self.out, """
    ==========================================================================
                          %(product)s Installation

                          version: %(version)s
                     machine type: %(mtype)s
                       OS version: %(os)s
                      destination: %(dest)s
                  # of processors: %(nproc)s
    =========================================================================
    """ % { "product" : self.product_name,
            "version" : self.version,
            "mtype"   : self.mtype,
            "os"      : get_os_version(),
            "dest"    : self.dest_dir,
            "nproc"   : self.options.nproc,
      }

  def machine_type(self):
    """
    Determine the mtype string.  The four pre-defined mtypes will be detected
    automatically (with Linux kernel 3.x treated as 2.6), but additional mtypes
    can be used if this method is re-implemented.
    """
    return machine_type()

  #---------------------------------------------------------------------
  # BINARY INSTALL
  #
  def install_from_binary (self) :
    """
    Unpackage the binary bundle in the destination directory.
    """
    # Copy base, build, and modules
    for i in ['base', 'build', 'modules', 'doc']:
      if os.path.exists(os.path.join(self.installer_dir, i)):
        copy_tree(os.path.join(self.installer_dir, i), os.path.join(self.dest_dir, i))

    # Reconfigure
    log = open(os.path.join(self.tmp_dir, "binary.log"), "w")
    if not os.path.exists(self.tmp_dir):
      os.makedirs(self.tmp_dir)
    if (sys.platform != "darwin") :
      os.environ["LD_LIBRARY_PATH"] = op.join(self.base_dir, "lib")
    self.product_specific_binary_install(log=log)

    os.chdir(self.build_dir)
    print >> self.out, "Configuring %s components..."%(self.product_name)
    self.reconfigure(log=log)

  #---------------------------------------------------------------------
  # SOURCE INSTALL
  #
  def install_from_source(self):
    log = self.out # open(os.path.join(self.tmp_dir, "source.log"), "w")
    call([
      'python',
      os.path.join('modules', 'cctbx_project', 'libtbx', 'auto_build', 'bootstrap.py'),
      '--builder', self.dest_dir_prefix,
      'base',
      'build'
    ], log=log)
    self.product_specific_source_install(log=log)
    self.install_from_binary()

    # raise NotImplementedError("Source installer returning soon.")
    # self.show_installation_paths()
    #
    # dependencies_dir = op.join(self.installer_dir, "dependencies")
    # bundle_dir = op.join(self.installer_dir, "bundles")
    # base_bundle = op.join(bundle_dir, "base.tar.gz")
    # build_bundle_file = op.join(bundle_dir, "build.tar.gz")
    # modules_bundle_file = op.join(bundle_dir, "modules.tar.gz")
    #
    # # PART 1: dependencies
    # # install Python 2.7 and other dependencies (wxPython, etc.)
    # if os.path.isfile(base_bundle):
    #   os.chdir(self.dest_dir)
    #   print >> self.out, "Using precompiled base bundle..."
    #   untar(base_bundle, log=self.out, verbose=True, change_ownership=False, check_output_path=False)
    # else :
    #   base_args = [
    #     # XXX will this allow spaces in path names?
    #     "--build_dir=%s" % self.dest_dir,
    #     "--tmp_dir=%s" % tmp_dir,
    #     "--pkg_dir=%s" % src_dir,
    #     "--pkg_dir=%s" % dependencies_dir,
    #     "--nproc=%s" % self.options.nproc,
    #     "--no-download",
    #   ] + self.base_package_options
    #   if (self.flag_build_gui):
    #     base_args.append("--gui")
    #   if (self.options.debug):
    #     base_args.append("--debug")
    #   if (not self.options.python_static):
    #     base_args.append("--python-shared")
    #   install_base_packages.installer(args=base_args, log=self.out)
    #
    # python_bin = op.join(self.base_dir, "bin", "python")
    # assert op.isfile(python_bin)
    # self.print_header('Base installation complete')
    #
    # # PART 2: product packages
    # os.chdir(self.installer_dir)
    # print >> self.out, "Building core %s components" % self.product_name
    # if (self.options.base_only) :
    #   print >> self.out, "--base-only was specified, so skipping %s libraries" % \
    #     self.product_name
    #   return
    #
    # for pkg_name in self.modules :
    #   pkg_file = op.join(self.src_dir, pkg_name + ".tar.gz")
    #   if (not op.isfile(pkg_file)) :
    #     raise RuntimeError("Can't find %s" % pkg_file)
    #   print >> self.out, "  unwrapping %s..." % op.basename(pkg_file)
    #   log_file = op.join(self.tmp_dir, "%s.log" % pkg_name)
    #   log = open(log_file, "w")
    #   os.chdir(self.modules_dir)
    #   untar(pkg_file, log=log, verbose=True, change_ownership=False,
    #     check_output_path=(not "cctbx_bundle" in pkg_file))
    #   log.close()
    #   if (pkg_file.startswith("cctbx_bundle")) :
    #     if (not op.isdir("cctbx_project")) :
    #       raise RuntimeError("Directory 'cctbx_project' not found!")
    #
    # tag_file = op.join(self.modules_dir, "TAG")
    # if op.isfile(tag_file) :
    #   os.rename(tag_file, op.join(self.modules_dir, "cctbx_bundle_TAG"))
    # self.product_specific_setup_before_compile(log=self.out)
    # log_file = op.join(self.tmp_dir, "compile_%s.log" % self.dest_dir_prefix)
    # log = open(log_file, "w")
    # os.chdir(self.build_dir)
    # if op.isfile("libtbx_refresh_is_completed") :
    #   os.remove("libtbx_refresh_is_completed")
    # config_path = op.join(self.modules_dir, "cctbx_project", "libtbx",
    #   "configure.py")
    # config_args = [
    #   "--scan_boost",
    #   "--command_version_suffix=%s" % self.version,
    #   "--current_working_directory=\"%s\"" % self.build_dir,
    # ] + self.configure_modules
    # if self.options.openmp :
    #   config_args.append("--enable-openmp-if-possible=True")
    # if (self.options.debug) :
    #   print >> self.out, ("configuring %s components with debugging..." % \
    #     self.product_name),
    #   config_args.insert(0, "--build=debug")
    # elif (self.options.no_opt) :
    #   print >> self.out, ("configuring %s components with no optimization..." % \
    #     self.product_name),
    #   config_args.insert(0, "--build=quick")
    # else :
    #   print >> self.out, ("configuring %s components with optimization..." % \
    #     self.product_name),
    #   config_args.insert(0, "--build=release")
    # call(" ".join([python_bin, config_path] + config_args), log=log)
    # print >> self.out, "ok"
    # # XXX in the original shell script, I had to run 'libtbx.configure reel'
    # # to get REEL actually configured - need to double-check this
    # scons_args = ["."]
    # if (not os.environ.get("LIBTBX_FULL_TESTING") in [None, ""]) :
    #   scons_args.append("-k")
    # scons_args.extend(["-j", str(self.options.nproc)])
    # print >> self.out, ("compiling %s components - this might take an hour..." % \
    #   self.product_name),
    # scons_bin = op.join(self.build_dir, "bin", "libtbx.scons")
    # assert op.isfile(scons_bin)
    # call(" ".join([scons_bin] + scons_args), log=log)
    # print >> self.out, "done"
    # self.product_specific_source_install(log=log)

  def show_installation_paths (self) :
    print >> self.out, """
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
  # NON-COMPILED COMPONENTS AND FINAL SETUP
  #
  def reconfigure (self, log) :
    """
    Run libtbx/configure.py to configure the build in the new location.
    """
    os.chdir(self.build_dir)
    args = [
      os.path.join(self.base_dir, 'bin', 'python'),
      os.path.join(self.modules_dir, 'cctbx_project', 'libtbx', 'configure.py'),
      "--current_working_directory", self.build_dir,
      "--command_version_suffix", self.version,
    ] + self.configure_modules
    if 'win32'==sys.platform:
      args = [
        os.path.join(self.base_dir, 'bin', 'python', 'python.exe'),
        os.path.join(self.modules_dir, 'cctbx_project', 'libtbx', 'configure.py'),
        "--command_version_suffix", self.version,
      ] + self.configure_modules

    print self.build_dir
    print args

    if 1: #try :
      call(args=args, log=log)
    else: #except RuntimeError :
      raise InstallerError("Configuration step incomplete!  See the log file for detailed error messages.")

  def install_finalize (self) :
    """
    Set up dispatchers and assorted shell scripts, create app bundles, etc.
    """
    self.print_header('Finalizing %s installation'%self.product_name)
    log_path = op.join(self.tmp_dir, "install_finalize.log")
    print >> self.out, "Log file: %s"%log_path
    log = open(log_path, "w")

    # Write environment files.
    self.write_environment_files()

    # Regenerate module files.
    if (self.flag_build_gui) and (sys.platform != "darwin") :
      os.environ["LD_LIBRARY_PATH"] = op.join(self.base_dir, "lib")
      regenerate_module_files.run(
        args=["--build_dir=%s" % self.dest_dir],
        out=self.out)

    # Write dispatcher_include file.
    print >> self.out, "Generating %s environment additions for dispatchers..." % \
      self.product_name
    fnsuffix = '.sh'
    envcmd = "export"
    if sys.platform == "win32":
      fnsuffix = '.bat'
      envcmd = "set"
    dispatcher = op.join(self.build_dir, "dispatcher_include_%s%s" %
      (self.dest_dir_prefix, fnsuffix))
    if (op.isfile(dispatcher)) :
      os.remove(dispatcher)
    env_prefix = self.product_name.upper() # e.g. "Phenix" -> "PHENIX"
    prologue = "\n".join([
      "%s %s=\"%s\"" % (envcmd, env_prefix, self.dest_dir),
      "%s %s_VERSION=%s" % (envcmd, env_prefix, self.version),
      "%s %s_ENVIRONMENT=1" % (envcmd, env_prefix),
      "%s %s_MTYPE=%s" % (envcmd, env_prefix, self.mtype),
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
    print 'Calling write_gui_dispatcher_include'
    print '  args %s' % dispatcher_opts
    print '  prologue %s' % prologue
    print '  epilogue %s' % epilogue
    write_gui_dispatcher_include.run(
      args=dispatcher_opts,
      prologue=prologue,
      epilogue=epilogue,
      out=self.out)
    assert op.isfile(dispatcher)

    # Run configure.py to generate dispatchers
    print >> self.out, "Configuring %s components..." % self.product_name
    os.chdir(self.build_dir)
    # ???
    if (op.exists("libtbx_refresh_is_completed")) :
      os.remove("libtbx_refresh_is_completed")
    self.reconfigure(log=log)
    os.chdir(self.build_dir)
    assert op.isfile("setpaths%s" %fnsuffix)
    if sys.platform != "win32":
      os.environ["PATH"] = "%s:%s" % (op.join(self.build_dir, "bin"), os.environ["PATH"])
    else:
      os.environ["PATH"] = "%s;%s" % (op.join(self.build_dir, "bin"), os.environ["PATH"])

    if not self.options.nopycompile:
      # Compile .py files
      print >> self.out, "Precompiling .py files..."
      os.chdir(self.modules_dir)
      call(args="libtbx.py_compile_all", log=log)

    # Copy README et al.
    for file_name in ["CHANGES", "LICENSE", "README", "README-DEV", "SOURCES"] :
      src_file = op.join(self.installer_dir, file_name)
      if op.exists(src_file) :
        dest_file = op.join(self.dest_dir, file_name)
        # XXX use our own implementation instead of shutil.copyfile
        if sys.platform == "win32" and src_file==dest_file:
          # writing to the same file on Windows renders it empty
          continue
        copy_file(src_file, dest_file)

    # generate .app (Mac only)
    apps_built = False
    if ((sys.platform == "darwin") and (len(self.make_apps) > 0) and
        (not self.options.no_app) and self.flag_build_gui) :
      os.chdir(self.build_dir)
      for app_name in self.make_apps :
        args = [
          "libtbx.create_mac_app",
          app_name,
          "--app_name=%s-%s" % (app_name, self.version),
          "--dest=%s" % self.dest_dir,
          "--alias_build"
        ]
        print >> self.out, "Generating Mac app launcher for %s..."%app_name
        try :
          call(args=" ".join(args), log=log)
        except RuntimeError, e :
          print "  ERROR:"
          print "  " + str(e)
          print "  installation will continue anyway."
        else :
          app_file = op.join(self.dest_dir, "%s-%s.app" %
            (app_name, self.version))
          if (not op.exists(app_file)) :
            print >> self.out, " failed."
            app_file = None
          else :
            apps_built = True

    # run custom finalization
    self.product_specific_finalize_install(log)

    if sys.platform == "win32":
      return

    # remove source files if desired
    if (self.options.compact):
      self.reduce_installation_size()

    self.display_final_message()

    # Fix permissions
    call([
      'chmod',
      '-R',
      'u+rw,a+rX',
      self.dest_dir
      ])

    # Show the app.
    if apps_built and (not "SSH_CLIENT" in os.environ):
      try:
        call(args=["open", self.dest_dir], log=self.out)
      except Exception, e:
        # Will fail in non-interactive environments.
        pass

  #---------------------------------------------------------------------
  # STUBS FOR SUBCLASSABLE METHODS

  def write_environment_files (self) :
    """
    Generate shell scripts in the top-level installation directory that can
    be used to set up the user environment to run the software.  This is
    implemented as a separate method because some products (e.g. Phenix) may
    have their own needs.
    """
    # csh/tcsh environment setup file
    print >> self.out, "Generating %s environment setup scripts..."%self.product_name
    env_prefix = self.product_name.upper() # e.g. "Phenix" -> "PHENIX"
    if sys.platform == "win32":
      f = open(os.path.join(self.dest_dir, '%s_env.bat'%self.dest_dir_prefix), 'w')
      f.write("set %s=%s\n" % (env_prefix, self.dest_dir))
      f.write("set %s_VERSION=%s\n" % (env_prefix, self.version))
      f.write("call %%%s%%\\build\\setpaths.bat\n" % (env_prefix))
      f.close()
      return

    f = open(os.path.join(self.dest_dir, '%s_env.csh'%self.dest_dir_prefix), 'w')
    f.write("#!/bin/csh -f\n")
    f.write("setenv %s \"%s\"\n" % (env_prefix, self.dest_dir))
    f.write("setenv %s_VERSION %s\n" % (env_prefix, self.version))
    f.write("source $%s/build/setpaths.csh\n" % (env_prefix))
    f.close()

    f = open(os.path.join(self.dest_dir, '%s_env.sh'%self.dest_dir_prefix), 'w')
    f.write("#!/bin/sh\n")
    f.write("#\n")
    f.write("export %s=\"%s\"\n" % (env_prefix, self.dest_dir))
    f.write("export %s_VERSION=%s\n" % (env_prefix, self.version))
    f.write(". $%s/build/setpaths.sh\n" % (env_prefix))
    f.close()

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

  def reduce_installation_size (self) :
    """
    Remove all files not required for program execution to save disk space,
    such as any C++ sources.  This can potentially save well over 100MB and is
    recommended for purely user-facing packages, but may make it more difficult
    to use packages for development purposes.
    """
    print >> self.out, "Removing unnecessary files to reduce disk usage...",
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
              print >> self.out, "  WARNING: error removing %s" % full_path
              print >> self.out, str(e)
            else :
              n_deleted += 1
    n_deleted_other = self.product_specific_reduce_installation_size()
    if (n_deleted_other is not None) :
      n_deleted += n_deleted_other
    print >> self.out, "%d files deleted" % n_deleted

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

  def print_header(self, msg):
    print >> self.out, ""
    print >> self.out, "*"*(len(msg) + 4)
    print >> self.out, "* %s *"%msg
    print >> self.out, "*"*(len(msg) + 4)
    print >> self.out, ""

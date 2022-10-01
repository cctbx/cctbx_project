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

  class installer(install_distribution.installer):
    product_name = "PHENIX"
    dest_dir_prefix = "phenix"
    make_apps = ["phenix"]
    installer_dir = os.path.dirname(os.path.dirname(__file__))

  if (__name__ == "__main__"):
    installer(sys.argv[1:])

Note that the location of the install script is not mandatory, but it must be
able to correctly detect the base installer directory.
"""

from __future__ import absolute_import, division, print_function
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
if (not libtbx_path in sys.path):
  sys.path.append(libtbx_path)
from libtbx.auto_build import write_gui_dispatcher_include
from libtbx.auto_build import regenerate_module_files
from libtbx.auto_build import install_base_packages
from libtbx.auto_build import create_mac_app
from libtbx.auto_build.installer_utils import *

class InstallerError(Exception):
  __orig_module__ = __module__
  # trick to get just "Sorry" instead of "libtbx.utils.Sorry"
  __module__ = Exception.__module__

  def reset_module(self):
    """
    Reset the class module on an instance to libtbx.utils.
    """
    self.__class__.__module__ = self.__class__.__orig_module__

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
    "intel-windows",
    "intel-windows-x86_64",
  ]
  # Boolean flags to control installer behaviour in a fine-grained manner
  flags = [
    'create_versioned_dispatchers'
  ]

  def __init__(self, args=None, out=sys.stdout):
    self.args = args or []
    self.out = out
    self.parse_options()

  def parse_options(self):
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
    parser.add_option("--verbose", "-v", dest="verbose", action="store_true",
      help="Provide more detailed output during installation", default=False)
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
    parser.add_option("--use-conda", default=False, action="store_true",
      help="Install conda dependencies (source install only)")
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
    self.product_specific_preinstallation_hook()
    self.check_directories()
    if sys.platform != "win32":
      self.print_banner()
      if not os.path.exists('build'):
        print("No build directory exists; trying source installation.", file=self.out)
        self.options.source = True
      if self.options.source:
        print("Source installation specified.", file=self.out)
        self.install_from_source()
      else:
        self.install_from_binary()
    self.install_finalize()
    self.print_header('Installation complete!')

  def basic_setup(self):
    # Check version
    self.version = self.get_version()
    assert (self.version is not None)
    # GUI Flag
    self.flag_build_gui = False
    if (self.include_gui_packages):
      self.flag_build_gui = not self.options.no_gui

    # Check this is a supported architecture.
    self.mtype = self.machine_type()
    if (self.mtype is None):
      raise InstallerError("Machine type not recognized")
    elif ((not self.mtype in self.supported_mtypes) and
          (not self.options.try_unsupported)):
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
      # check for conda
      if os.path.isdir(op.join(self.installer_dir, "conda_base")) or \
         os.path.exists(op.join(self.installer_dir, "conda_base.tar")):
        self.base_dir = op.join(self.dest_dir, "conda_base")
      return
    # The default behavior for nearly all program's --prefix options
    # is to create the directory, so I don't think the --makedirs option
    # is necessary.

    if not os.access(self.options.prefix, os.W_OK):
      if not os.path.exists(self.options.prefix):
        os.makedirs(self.options.prefix)
      else:
        raise InstallerError("""
  Installation directory not writeable:
    %s
  Please specify an alternative directory using the --prefix option."""
      %self.options.prefix)

    # Do not overwrite an existing installation.
    self.dest_dir = op.abspath(op.join(
      self.options.prefix, "%s-%s"%(self.dest_dir_prefix, self.version)))
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

    # check for conda
    if os.path.isdir(op.join(self.installer_dir, "conda_base")) or \
       os.path.exists(op.join(self.installer_dir, "conda_base.tar")):
      self.base_dir = op.join(self.dest_dir, "conda_base")

    # Environment variables required by other scripts
    os.environ["%s_MTYPE" % self.product_name] = self.mtype
    os.environ["%s_INSTALLER" % self.product_name] = self.installer_dir
    os.environ["%s_LOC" % self.product_name] = self.dest_dir
    os.environ["%s_BUILD" % self.product_name] = self.build_dir

  def print_banner(self):
    print("""
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
      }, file=self.out)

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
  def install_from_binary(self):
    """
    Unpackage the binary bundle in the destination directory.
    """
    # Copy base, build, and modules
    for i in ['base', 'conda_base', 'build', 'modules', 'doc']:
      if os.path.exists(os.path.join(self.installer_dir, i)):
        copy_tree(os.path.join(self.installer_dir, i), os.path.join(self.dest_dir, i))

    # Copy conda_base packaged with conda-pack if available
    # only this file or the conda_base directory will exist in the installer
    conda_base_tarfile = os.path.join(self.installer_dir, 'conda_base.tar')
    if os.path.exists(conda_base_tarfile):
      import subprocess
      import tarfile
      tarball = tarfile.open(conda_base_tarfile)
      dest_dir = os.path.join(self.dest_dir, 'conda_base')
      tarball.extractall(path=dest_dir)
      tarball.close()
      cwd = os.getcwd()
      os.chdir(dest_dir)
      unpack_cmd = [sys.executable, os.path.join('.', 'bin', 'conda-unpack')]
      if sys.platform == 'win32':
        unpack_cmd = [os.path.join('.', 'Scripts', 'conda-unpack.exe')]
      subprocess.check_call(unpack_cmd)
      os.chdir(cwd)

    # Reconfigure
    log = open(os.path.join(self.tmp_dir, "binary.log"), "w")
    if not os.path.exists(self.tmp_dir):
      os.makedirs(self.tmp_dir)
    if (sys.platform != "darwin"):
      os.environ["LD_LIBRARY_PATH"] = "lib:%s:%s" % \
        (op.join(self.base_dir, "lib"), op.join(self.base_dir, "lib64"))
    self.product_specific_binary_install(log=log)

    os.chdir(self.build_dir)
    print("Configuring %s components..."%(self.product_name), file=self.out)
    self.reconfigure(log=log)

  #---------------------------------------------------------------------
  # SOURCE INSTALL
  #
  def install_from_source(self):
    log = self.out # open(os.path.join(self.tmp_dir, "source.log"), "w")
    cmd = [
      sys.executable,
      os.path.join('modules', 'cctbx_project', 'libtbx', 'auto_build', 'bootstrap.py'),
      'base',
      'build',
      '--builder', self.dest_dir_prefix,
      '--nproc', str(self.options.nproc),
    ]
    if self.options.use_conda:
      cmd.append('--use-conda')
      self.base_dir = op.join(self.dest_dir, "conda_base")
    call(cmd, log=log)
    self.product_specific_source_install(log=log)
    self.install_from_binary()

  def show_installation_paths(self):
    print("""
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
        "tmp_dir" : self.tmp_dir, }, file=self.out)

  #---------------------------------------------------------------------
  # NON-COMPILED COMPONENTS AND FINAL SETUP
  #
  def reconfigure(self, log):
    """
    Run libtbx/configure.py to configure the build in the new location.
    """
    os.chdir(self.build_dir)

    base_python = os.path.join(self.base_dir, 'bin', 'python')
    if 'win32' in sys.platform:
      base_python = os.path.join(self.base_dir, 'bin', 'python', 'python.exe')

    # check for conda
    if self.base_dir.endswith('conda_base'):
      if 'darwin' in sys.platform:
        base_python = os.path.join(self.base_dir, 'python.app', 'Contents',
                                   'MacOS', 'python')
      elif 'win32' in sys.platform:
        base_python = os.path.join(self.base_dir, 'python.exe')

    args = [
      base_python,
      os.path.join(self.modules_dir, 'cctbx_project', 'libtbx', 'configure.py'),
      "--current_working_directory", self.build_dir
    ]
    if 'win32'==sys.platform:
      args = [
        base_python,
        os.path.join(self.modules_dir, 'cctbx_project', 'libtbx', 'configure.py'),
      ]
    if 'create_versioned_dispatchers' in self.flags:
      args += [ "--command_version_suffix", self.version ]
    args += self.configure_modules

    # check for conda
    if self.base_dir.endswith('conda_base'):
      args += ['--use_conda']

    if self.options.verbose:
      print(self.build_dir)
      print(args)

    if 1: #try :
      call(args=args, log=log, verbose=self.options.verbose)
    else: #except RuntimeError :
      raise InstallerError("Configuration step incomplete!  See the log file for detailed error messages.")

  def install_finalize(self):
    """
    Set up dispatchers and assorted shell scripts, create app bundles, etc.
    """
    self.print_header('Finalizing %s installation'%self.product_name)
    log_path = op.join(self.tmp_dir, "install_finalize.log")
    print("Log file: %s"%log_path, file=self.out)
    log = open(log_path, "w")

    # Write environment files.
    self.write_environment_files()

    # Regenerate module files.
    if (self.flag_build_gui) and (sys.platform != "darwin") and \
      (not self.base_dir.endswith('conda_base')):
      os.environ["LD_LIBRARY_PATH"] = "lib:%s:%s" % \
        (op.join(self.base_dir, "lib"), op.join(self.base_dir, "lib64"))
      regenerate_module_files.run(
        os.path.join(self.dest_dir, 'base'),
        out=self.out)

    # Write dispatcher_include file.
    print("Generating %s environment additions for dispatchers..." % \
      self.product_name, file=self.out)
    fnsuffix = '.sh'
    envcmd = "export"
    if sys.platform == "win32":
      fnsuffix = '.bat'
      envcmd = "set"
    dispatcher = op.join(self.build_dir, "dispatcher_include_%s%s" %
      (self.dest_dir_prefix, fnsuffix))
    if (op.isfile(dispatcher)):
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
    if (not self.flag_build_gui):
      dispatcher_opts.append("--ignore_missing_dirs")
    # check for conda
    if (self.base_dir.endswith('conda_base')):
      dispatcher_opts += ["--use_conda", "--ignore_missing_dirs"]
    # FIXME this will happen regardless of whether the GUI modules are being
    # distributed or not - will this be problematic?
    print('Calling write_gui_dispatcher_include')
    print('  args %s' % dispatcher_opts)
    print('  prologue %s' % prologue)
    print('  epilogue %s' % epilogue)
    write_gui_dispatcher_include.run(
      args=dispatcher_opts,
      prologue=prologue,
      epilogue=epilogue,
      out=self.out)
    assert op.isfile(dispatcher)

    # Run configure.py to generate dispatchers
    print("Configuring %s components..." % self.product_name, file=self.out)
    os.chdir(self.build_dir)
    # ???
    if (op.exists("libtbx_refresh_is_completed")):
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
      print("Precompiling .py files...", file=self.out)
      os.chdir(self.modules_dir)
      call(args="libtbx.py_compile_all -i", log=log)

    # Copy README et al.
    for file_name in ["CHANGES", "LICENSE", "README", "README-DEV", "SOURCES"] :
      src_file = op.join(self.installer_dir, file_name)
      if op.exists(src_file):
        dest_file = op.join(self.dest_dir, file_name)
        # XXX use our own implementation instead of shutil.copyfile
        if sys.platform == "win32" and src_file==dest_file:
          # writing to the same file on Windows renders it empty
          continue
        copy_file(src_file, dest_file)

    # generate .app (Mac only)
    apps_built = False
    if ((sys.platform == "darwin") and (len(self.make_apps) > 0) and
        (not self.options.no_app) and self.flag_build_gui):
      os.chdir(self.build_dir)
      for app_name in self.make_apps :
        args = [
          "libtbx.create_mac_app",
          app_name,
          "--app_name=%s-%s" % (app_name, self.version),
          "--dest=%s" % self.dest_dir,
          "--alias_build"
        ]
        print("Generating Mac app launcher for %s..."%app_name, file=self.out)
        try :
          call(args=" ".join(args), log=log)
        except RuntimeError as e:
          print("  ERROR:")
          print("  " + str(e))
          print("  installation will continue anyway.")
        else :
          app_file = op.join(self.dest_dir, "%s-%s.app" %
            (app_name, self.version))
          if (not op.exists(app_file)):
            print(" failed.", file=self.out)
            app_file = None
          else :
            apps_built = True

    # run custom finalization
    self.product_specific_finalize_install(log)

    # remove source files if desired
    if (self.options.compact):
      self.reduce_installation_size()

    self.display_final_message()

    if sys.platform == "win32":
      return

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
      except Exception:
        # Will fail in non-interactive environments.
        pass

  #---------------------------------------------------------------------
  # STUBS FOR SUBCLASSABLE METHODS

  def write_environment_files(self):
    """
    Generate shell scripts in the top-level installation directory that can
    be used to set up the user environment to run the software.  This is
    implemented as a separate method because some products (e.g. Phenix) may
    have their own needs.
    """
    # csh/tcsh/bat environment setup file
    print("Generating %s environment setup scripts..."%self.product_name, file=self.out)
    env_prefix = self.product_name.upper() # e.g. "Phenix" -> "PHENIX"
    if sys.platform == "win32":
      f = open(os.path.join(self.dest_dir, '%s_env.bat'%self.dest_dir_prefix), 'w')
      f.write("@echo off\n")
      # Use the %~dp0 alias for specifying the full path to where phenix_env.bat will be located.
      # This presumes phenix_env.bat will always reside where it has originally been installed to.
      f.write("set %s=%%~dp0\n" %env_prefix)
      f.write("set %s_VERSION=%s\n" % (env_prefix, self.version))
      f.write("call \"%%%s%%\\build\\setpaths.bat\"\n" % (env_prefix))
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

  def get_version(self):
    """
    Determine the version suffix (if any) for the destination directory.  This
    must be implemented by subclasses unless a file named VERSION is present
    at the top level of the installer folder.
    """
    version_file = op.join(self.installer_dir, "VERSION")
    if op.isfile(version_file):
      return open(version_file).read().strip()
    return NotImplementedError()

  def reduce_installation_size(self):
    """
    Remove all files not required for program execution to save disk space,
    such as any C++ sources.  This can potentially save well over 100MB and is
    recommended for purely user-facing packages, but may make it more difficult
    to use packages for development purposes.
    """
    print("Removing unnecessary files to reduce disk usage...", end=' ', file=self.out)
    # XXX should this include .o files?
    remove_extensions = [".cpp", ".cc", ".hpp", ".h", ".hh"]
    n_deleted = 0
    for dirname, dirnames, filenames in os.walk(self.modules_dir):
      for file_name in filenames :
        for ext in remove_extensions :
          if file_name.endswith(ext):
            full_path = op.join(dirname, file_name)
            try :
              os.remove(full_path)
            except Exception as e:
              print("  WARNING: error removing %s" % full_path, file=self.out)
              print(str(e), file=self.out)
            else :
              n_deleted += 1
    n_deleted_other = self.product_specific_reduce_installation_size()
    if (n_deleted_other is not None):
      n_deleted += n_deleted_other
    print("%d files deleted" % n_deleted, file=self.out)

  def add_product_specific_options(self, parser):
    """
    Add command-line options specific to the distributed package.
    """
    pass

  def product_specific_prepackage_hook(self, directory):
    """
    Modify files, etc. before the installer package is created.

    :param directory: base directory of the installer package
    """
    pass

  def product_specific_preinstallation_hook(self):
    """
    Perform additional checks on parsed command line options or the
    destination machine before any installation action takes place.
    """
    pass

  def product_specific_binary_install(self, log):
    """
    Perform additional actions required for the binary installation of a
    specific product.
    """
    pass

  def product_specific_setup_before_compile(self, log):
    """
    Perform any necessary modifications to the sources prior to compilation.
    """
    pass

  def product_specific_source_install(self, log):
    """
    Build additional sources, e.g. Rosetta
    """
    pass

  def product_specific_dispatcher_prologue(self):
    """
    Environment modifications to be included near the start of the dispatchers.
    """
    return []

  def product_specific_dispatcher_epilogue(self):
    """
    Environment modifications to be included at the end of the dispatchers.
    """
    return []

  def product_specific_finalize_install(self, log):
    """
    Additional installation setup, file cleanup, more add-ons, etc.
    """
    pass

  def product_specific_reduce_installation_size(self, log):
    """
    Remove unused files specific to this product, and return the number deleted.
    """
    return 0

  def display_final_message(self):
    """
    Final instructions for user, etc.
    """
    pass

  def print_header(self, msg):
    print("", file=self.out)
    print("*"*(len(msg) + 4), file=self.out)
    print("* %s *"%msg, file=self.out)
    print("*"*(len(msg) + 4), file=self.out)
    print("", file=self.out)

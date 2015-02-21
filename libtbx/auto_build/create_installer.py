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
import shutil
import time
import stat
import os
import sys
import subprocess
import imp

import libtbx.auto_build.rpath

# XXX HACK
libtbx_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
if (not libtbx_path in sys.path):
  sys.path.append(libtbx_path)

INSTALL_SH = """\
#!/bin/bash
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
"""

BASHRC = """\
#!/bin/bash
export %(env_prefix)s=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
export %(env_prefix)s_VERSION=%(version)s
source $%(env_prefix)s/build/setpaths.sh
"""

CSHRC = """\
#!/bin/csh -f
set testpath=($_)
if ("$testpath" != "") then
    set testpath=$testpath[2]
else
    set testpath=$0
endif
set rootdir=`dirname $testpath`
set fullpath=`cd $rootdir && pwd -P`
setenv LIBTBX_BUILD_RELOCATION_HINT $fullpath
setenv %(env_prefix)s $fullpath
setenv %(env_prefix)s_VERSION %(version)s
source $%(env_prefix)s/build/setpaths.csh
"""

def makedirs(path):
  try:
    os.makedirs(path)
  except Exception, e:
    print "Directory already exists: %s"%path

def archive(source, destination, tarfile=None):
  if source == destination:
    print "Source and destination are the same, skipping: %s"%source
    return
  assert not os.path.exists(destination), "File exists: %s"%destination
  print "Copying: %s -> %s"%(source, destination)
  if not os.path.exists(source):
    print "Warning: source does not exist! Skipping: %s"%source
    return
  shutil.copytree(
    source,
    destination,
    ignore=shutil.ignore_patterns('*.pyc', '*.pyo', '.svn', '.git', '.swp', '.sconsign', '.o'),
    symlinks=True
    )

def tar(source, tarfile, cwd=None):
  assert not os.path.exists(tarfile), "File exists: %s"%tarfile
  print "Archiving: %s -> %s"%(source, tarfile)
  subprocess.check_call([
      'tar',
      '-cz',
      '-f', tarfile,
      source
    ],
    cwd=cwd)

class SetupInstaller(object):
  def __init__(self, **kwargs):
    self.install_script = kwargs.get('install_script')
    self.version = kwargs.get('version')
    self.script = kwargs.get('script')
    #
    self.dest_dir = os.path.abspath(kwargs.get('dest_dir'))
    self.root_dir = os.path.abspath(kwargs.get('root_dir') or os.getcwd())
    self.dist_dir = os.path.abspath(kwargs.get('dist_dir') or os.path.join(self.root_dir, 'dist'))

    self.license = kwargs.get('license')
    if self.license:
      self.license = os.path.abspath(self.license)
    self.readme = kwargs.get('readme') or [os.path.join(libtbx_path, 'COPYRIGHT_2_0.txt')]
    if self.readme:
      self.readme = [os.path.abspath(i) for i in self.readme]

    # Is this a source or binary installer?
    self.binary = kwargs.get('binary')

    # Load the installer class, get the list of modules.
    assert os.path.isfile(self.install_script)
    installer_module = imp.load_source('install_script', self.install_script)
    self.installer = installer_module.installer()

  def run(self):
    # Setup directory structure
    print "Installer will be %s"%self.dest_dir
    assert not os.path.exists(self.dest_dir), "Installer dir exists: %s"%self.dest_dir
    makedirs(self.dest_dir)
    for i in ['bin', 'lib']:
      makedirs(os.path.join(self.dest_dir, i))
    self.copy_info()
    self.write_environment_files()
    self.copy_libtbx()
    self.copy_doc()
    self.copy_modules()
    if self.binary:
      self.copy_dependencies()
      self.copy_build()
    self.fix_permissions()
    self.make_dist()
    if self.binary and sys.platform == "darwin":
      self.make_dist_pkg()

  def copy_info(self):
    # Basic setup #
    # Write VERSION
    with open(os.path.join(self.dest_dir, 'VERSION'), 'w') as f:
      f.write(self.version)
    # Write README
    for i in self.readme:
      shutil.copyfile(i, os.path.join(self.dest_dir, os.path.basename(i)))
    # Write LICENSE
    if os.path.isfile(self.license):
      shutil.copyfile(self.license, os.path.join(self.dest_dir, 'LICENSE'))
    # Actual Python installer script
    shutil.copyfile(self.install_script, os.path.join(self.dest_dir, 'bin', 'install.py'))
    # Write executable Bash script wrapping Python script
    with open(os.path.join(self.dest_dir, 'install'), 'w') as f:
      f.write(INSTALL_SH)
    st = os.stat(os.path.join(self.dest_dir, "install"))
    os.chmod(os.path.join(self.dest_dir, "install"), st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

  def copy_libtbx(self):
    # Copy over libtbx for setup.
    archive(
      os.path.join(libtbx_path),
      os.path.join(self.dest_dir, 'lib', 'libtbx')
    )

  def copy_dependencies(self):
    # Copy dependencies
    archive(
      os.path.join(self.root_dir, 'base'),
      os.path.join(self.dest_dir, 'base')
    )
    libtbx.auto_build.rpath.run(['--otherroot', os.path.join(self.root_dir, 'base'), os.path.join(self.dest_dir, 'base')])

  def copy_build(self):
    # Copy the entire build directory, minus .o files.
    archive(
      os.path.join(self.root_dir, 'build'),
      os.path.join(self.dest_dir, 'build')
    )
    libtbx.auto_build.rpath.run(['--otherroot', os.path.join(self.root_dir, 'base'), os.path.join(self.dest_dir, 'build')])

  def copy_modules(self):
    # Source modules #
    for module in set(self.installer.modules):
      archive(
        os.path.join(self.root_dir, 'modules', module),
        os.path.join(self.dest_dir, 'modules', module)
      )

  def copy_doc(self):
    # Copy doc
    archive(
      os.path.join(self.root_dir, 'doc'),
      os.path.join(self.dest_dir, 'doc')
    )

  def fix_permissions(self):
    subprocess.check_call([
      'chmod',
      '-R',
      'u+rw,a+rX',
      self.dest_dir
      ])

  def write_environment_files(self):
    """Generate shell scripts in the top-level installation directory."""
    print "Generating %s environment setup scripts..."%self.installer.product_name
    fmt = {'env_prefix':self.installer.product_name.upper(), 'version':self.version}
    # bash
    with open(os.path.join(self.dest_dir, '%s_env.sh'%self.installer.product_name.lower()), 'w') as f:
      f.write(BASHRC%fmt)
    # tcsh
    with open(os.path.join(self.dest_dir, '%s_env.csh'%self.installer.product_name.lower()), 'w') as f:
      f.write(CSHRC%fmt)

  def make_dist(self):
    makedirs(self.dist_dir)
    tar(
      os.path.basename(self.dest_dir),
      os.path.join(self.dist_dir, '%s.tar.gz'%os.path.basename(self.dest_dir)),
      cwd=os.path.join(self.dest_dir, '..')
    )

  def make_dist_pkg(self):
    if (not os.access("/Applications", os.W_OK|os.X_OK)) :
      print "Can't access /Applications - skipping .pkg build"
      return

    # This has to match other convetions...
    pkg_prefix = "/Applications"
    app_root_dir = os.path.join(pkg_prefix,
                                '%s-%s'%(self.installer.dest_dir_prefix,
                                         self.version))
    print self.dest_dir
    print os.path.exists(os.path.join(self.dest_dir, "install"))
    subprocess.check_call([
      os.path.join(self.dest_dir, 'install'),
      '--prefix', pkg_prefix,
    ], cwd=self.dest_dir)

    tmp = os.path.join(self.root_dir, 'tmp')
    makedirs(tmp)
    makedirs(self.dist_dir)
    os.chdir(tmp) # UGH X 1000.
    from libtbx.auto_build import create_mac_pkg
    create_mac_pkg.run(args=[
        "--package_name", self.installer.product_name,
        "--organization", self.installer.organization,
        "--version", self.version,
        "--license", self.license,
        "--dist-dir", self.dist_dir,
        #"--no_compress",
        app_root_dir
    ])

def run(args):
  parser = OptionParser()
  parser.add_option("--version", dest="version", action="store",
    help="Package version", default=time.strftime("%Y_%m_%d",time.localtime()))
  parser.add_option("--binary", dest="binary", action="store_true",
    help="Include base and build directories", default=False)
  parser.add_option("--root_dir", dest="root_dir", action="store",
    help="Environment root")
  parser.add_option("--dist_dir", dest="dist_dir", action="store",
    help="Archive output directory")
  parser.add_option("--readme", dest="readme", action="append",
    help="Readme file", default=[])
  parser.add_option("--license", dest="license", action="store",
    help="License file", default=os.path.join(libtbx_path, "LICENSE_2_0.txt"))
  parser.add_option("--install_script", dest="install_script",
    help="Final installation script", default=None, metavar="FILE")
  options, args_ = parser.parse_args(args=args)
  assert len(args_) == 1, "Destination directory required argument."
  setup = SetupInstaller(
    dest_dir=args_[0],
    root_dir=options.root_dir,
    dist_dir=options.dist_dir,
    version=options.version,
    readme=options.readme,
    license=options.license,
    install_script=options.install_script,
    binary=options.binary,
  )
  setup.run()

if (__name__ == "__main__") :
  sys.exit(run(sys.argv[1:]))

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

def archive(source, destination, tarfile=None):
  assert not os.path.exists(destination), "File exists: %s"%destination
  print "Copying: %s -> %s"%(source, destination)
  if not os.path.exists(source):
    print "Warning: source does not exist! Skipping: %s"%source
    return
  shutil.copytree(
    source,
    destination,
    ignore=shutil.ignore_patterns('*.pyc', '*.pyo', '.svn', '.git', '.swp', '.sconsign'),
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
    self.modules = set(kwargs.get('modules') or [])
    self.base_modules = set(kwargs.get('base_modules') or [])
    # 
    self.root = os.path.abspath(kwargs.get('root'))
    self.dest = os.path.abspath(kwargs.get('dest'))
    dest_tar = kwargs.get('dest_tar') or '%s.tar.gz'%os.path.basename(self.dest)
    self.dest_tar = os.path.abspath(dest_tar)

    self.license = kwargs.get('license')
    self.readme = kwargs.get('readme') or [os.path.join(libtbx_path, 'COPYRIGHT_2_0.txt')]
    # Load the installer class, get the list of modules.
    assert os.path.isfile(self.install_script)
    installer_module = imp.load_source('install_script', self.install_script)
    installer = installer_module.installer()
    self.modules |= set(installer.modules)
    self.base_modules |= set(installer.base_modules)

  def run(self):
    # Setup directory structure
    print "Installer will be %s"%self.dest
    assert not os.path.exists(self.dest), "Installer dir exists: %s"%self.dest
    os.makedirs(self.dest)
    for i in ['bin', 'lib']:
      os.makedirs(os.path.join(self.dest, i))
    self.copy_info()
    self.copy_libtbx()
    self.copy_dependencies()
    self.copy_build()
    self.copy_modules()
    self.copy_doc()
    self.make_dist()
    
  def make_dist(self):
    try:
      os.makedirs(os.path.dirname(self.dest_tar))
    except Exception, e:
      print "Warning: %s"%e
    tar(
      os.path.basename(self.dest),
      self.dest_tar,
      cwd=os.path.join(self.dest, '..')
    )

  def copy_info(self):
    # Basic setup #
    # Write VERSION
    with open(os.path.join(self.dest, 'VERSION'), 'w') as f:
      f.write(self.version)
    # Write README
    for i in self.readme:
      shutil.copyfile(i, os.path.join(self.dest, os.path.basename(i)))
    # Write LICENSE
    if os.path.isfile(self.license):
      shutil.copyfile(self.license, os.path.join(self.dest, 'LICENSE'))
    # Actual Python installer script
    shutil.copyfile(self.install_script, os.path.join(self.dest, 'bin', 'install.py'))
    # Write executable Bash script wrapping Python script
    with open(os.path.join(self.dest, 'install'), 'w') as f:
      f.write(INSTALL_SH)
    st = os.stat(os.path.join(self.dest, "install"))
    os.chmod(os.path.join(self.dest, "install"), st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
  
  def copy_libtbx(self):
    # Copy over libtbx for setup.
    archive(
      os.path.join(libtbx_path), 
      os.path.join(self.dest, 'lib', 'libtbx')
    )

  def copy_dependencies(self):
    # Copy dependencies
    archive(
      os.path.join(self.root, 'base'),
      os.path.join(self.dest, 'base')
    )
    libtbx.auto_build.rpath.run(['--otherroot', os.path.join(self.root, 'base'), os.path.join(self.dest, 'base')])

  def copy_build(self):
    # Compiled modules
    archive(
      os.path.join(self.root, 'build', 'lib'),
      os.path.join(self.dest, 'build', 'lib')
    )
    # executables
    build_dir = os.path.join(self.root, 'build')
    for j in [i for i in os.listdir(build_dir) if os.path.isdir(os.path.join(build_dir, i, "exe"))]:
      archive(
        os.path.join(self.root, 'build', j, 'exe'),
        os.path.join(self.dest, 'build', j, 'exe')
      )
    libtbx.auto_build.rpath.run(['--otherroot', os.path.join(self.root, 'base'), os.path.join(self.dest, 'build')])

  def copy_modules(self):
    # Source modules #
    for module in self.modules:
      archive(
        os.path.join(self.root, 'modules', module),
        os.path.join(self.dest, 'modules', module)
      )
      
  def copy_doc(self):
    # Copy doc
    archive(
      os.path.join(self.root, 'doc'),
      os.path.join(self.dest, 'doc')
    )
    
def run (args) :
  parser = OptionParser()
  parser.add_option("--version", dest="version", action="store",
    help="Package version", default=time.strftime("%Y_%m_%d",time.localtime()))
  parser.add_option("--binary", dest="binary", action="store_true",
    help="Setup for binary installer only (no source packages)", default=False)
  parser.add_option("--root", dest="root", action="store",
    help="Environment root", default=os.getcwd())
  parser.add_option("--dest", dest="dest", action="store",
    help="Destination folder", default=os.getcwd())
  parser.add_option("--dest_tar", dest="dest_tar", action="store",
    help="Tar archive filename")
  parser.add_option("--readme", dest="readme", action="append",
    help="Readme file", default=[])
  parser.add_option("--license", dest="license", action="store",
    help="License file", default=os.path.join(libtbx_path, "LICENSE_2_0.txt"))
  parser.add_option("--install_script", dest="install_script",
    help="Final installation script", default=None, metavar="FILE")
  parser.add_option("--module", dest="modules", action="append",
    help="Local modules to include")
  parser.add_option("--base-module", dest="base_modules", action="append",
    help="Additional local modules placed in base/ directory")
  options, args_ = parser.parse_args(args=args)
  setup = SetupInstaller(
    version=options.version,
    root=options.root,
    dest=options.dest,
    readme=options.readme,
    license=options.license,
    install_script=options.install_script,
    modules=options.modules,
    base_modules=options.base_modules,
    dest_tar=options.dest_tar
  )
  setup.run()

if (__name__ == "__main__") :
  sys.exit(run(sys.argv[1:]))

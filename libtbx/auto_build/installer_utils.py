
from __future__ import division
import warnings
import shutil
import time
import stat
import re
import os.path as op
import os
import sys

# XXX CCTBX itself requires at least Python 2.5 (and some packages such as
# Phenix require 2.7+), but this script is intended to bootstrap an
# installation on older systems as well
def check_python_version () :
  if (not sys.version_info >= (2,3)) :
    raise Exception("Python version 2.3 or greater required to run this "+
      "script.")
  elif (sys.version_info < (2,4)) : # subprocess module not available
    warnings.warn("You are running an obsolete version of Python; this script "+
      "should still run, but not all functionality is available.",
      DeprecationWarning)
    time.sleep(2)

def import_subprocess () :
  if (sys.version_info[1] >= 7) :
    import subprocess
  else :
    # XXX HACK
    libtbx_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    if (not libtbx_path in sys.path) :
      sys.path.append(libtbx_path)
    import subprocess_with_fixes as subprocess
  return subprocess

def call (args, log=sys.stdout, shell=True, cwd=None) :
  subprocess = import_subprocess()
  rc = None
  # shell=True requires string as args.
  if shell and isinstance(args, list) :
    args = " ".join(args)
  p = subprocess.Popen(
    args=args,
    shell=shell,
    cwd=cwd,
    bufsize=-1,
    stdin=None,
    stdout=log, #subprocess.PIPE,
    stderr=subprocess.STDOUT,
    universal_newlines=True,
    close_fds=False)
  #o, e = p.communicate()
  #log.write(o)
  log.flush()
  p.wait()
  log.flush()
  rc = p.returncode
  if (rc != 0) :
    raise RuntimeError("Call to '%s' failed with exit code %d" % (args, rc))

def check_output (*popenargs, **kwargs) :
  # Back-port of Python 2.7 subprocess.check_output.
  subprocess = import_subprocess()
  process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
  output, unused_err = process.communicate()
  retcode = process.poll()
  if retcode:
    raise RuntimeError("Call to '%s' failed with exit code %d" % (popenargs, retcode))
  return output

def untar (pkg_name, log=sys.stdout, verbose=False, change_ownership=False,
    check_output_path=True) :
  assert os.path.isfile(pkg_name), pkg_name
  verbose_flag = owner_flag = ""
  if (verbose) :
    verbose_flag = "v"
  if (change_ownership) :
    owner_flag = "o"
  cmd = [ "tar", "x%s%sf" % (owner_flag, verbose_flag) ]
  if (pkg_name.endswith("gz")) :
    cmd = ["tar", "zx%s%sf" % (owner_flag, verbose_flag) ]
  elif (pkg_name.endswith("bz2")) :
    cmd = ["tar", "jx%s%sf" % (owner_flag, verbose_flag) ]
  args = cmd + [pkg_name]
  call(" ".join(args), log)
  dir_name = re.sub(".tgz", "",
               re.sub(".tar.gz", "",
                 re.sub(".tar.bz2", "",
                   os.path.basename(pkg_name))))
  if (check_output_path) :
    if (not os.path.isdir(dir_name)) :
      time.sleep(1)
      if (not os.path.isdir(dir_name)) :
        raise RuntimeError("Expected directory '%s' not found!" % dir_name)
    return os.path.abspath(dir_name)
  return None

def detect_osx_version () :
  uname = os.uname()
  version = uname[2]
  major, minor, rev = version.split(".")
  return int(major)

def copy_file (src_path, dest_path, executable=None) :
  assert os.path.isfile(src_path)
  open(dest_path, "wb").write(open(src_path, "rb").read())
  if os.access(src_path, os.X_OK) or executable :
    mode = os.stat(dest_path).st_mode
    os.chmod(dest_path, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

# shutil.copytree replacement - handles symlinks intelligently and also
# circumvents this bug:
# http://bugs.python.org/issue14662
# this is not a general solution, but it is good enough for the installer
# process (at least on Unix).
def copy_tree (src_path, dest_path, verbose=False, log=sys.stdout) :
  if (sys.version_info > (2,7,4)) :
    shutil.copytree(src_path, dest_path, symlinks=True)
  else :
    assert os.path.isdir(src_path), src_path
    assert not os.path.exists(dest_path), dest_path
    if (verbose) :
      print >> log, "creating %s" % dest_path
    os.makedirs(dest_path)
    for path_name in os.listdir(src_path) :
      node_src_path = os.path.join(src_path, path_name)
      node_dest_path = os.path.join(dest_path, path_name)
      if os.path.islink(node_src_path) :
        target_path = os.readlink(node_src_path)
        os.symlink(target_path, node_dest_path)
      elif os.path.isfile(node_src_path) :
        if (verbose) :
          print >> log, "  copy %s -> %s" % (node_src_path, node_dest_path)
        copy_file(node_src_path, node_dest_path)
      elif os.path.isdir(node_src_path) :
        copy_tree(node_src_path, node_dest_path)
      else :
        if (verbose) :
          print >> log, "  skipping %s" % node_src_path

def get_os_version () :
  uname = os.uname()
  kernel_version = uname[2]
  if (uname[0] == "Darwin") :
    os_versions = {
      "14" : "10.10",
      "13" : "10.9",
      "12" : "10.8",
      "11" : "10.7",
      "10" : "10.6",
      "9"  : "10.5",
      "8"  : "10.4",
    }
    return os_versions.get(kernel_version.split(".")[0], "unknown")
  return kernel_version

def machine_type () :
  import os
  uname = os.uname()
  if (uname[0] == "Linux") :
    platform = "intel-linux-2.6"
  elif (uname[0] == "Darwin") :
    platform = "mac-intel-osx"
  else :
    platform = uname[0]
  if (uname[-1] == "x86_64") :
    platform += "-x86_64"
  elif (platform == "mac-intel-osx") :
    version_fields = uname[2].split(".")
    major_version = int(version_fields[0])
    if (major_version >= 10) :
      platform += "-x86_64"
  return platform

def regenerate_relative_symlinks (dir_name, log=sys.stdout) :
  """
  Rewrite all symlinks in a directory with relative paths (allows later
  relocation).  This is mostly just for Mac OS where there are a lot of links
  into the Python.framework bundle.
  """
  old_cwd = os.getcwd()
  os.chdir(dir_name)
  for file_name in os.listdir(dir_name) :
    # e.g. /usr/bin/cmd
    full_path = op.join(dir_name, file_name)
    if op.islink(full_path) :
      real_path = op.realpath(full_path) # e.g. /usr/cctbx/bin/cmd
      prefix = op.commonprefix([full_path, real_path]) # /usr
      rel_path_link = op.relpath(full_path, prefix) # bin/cmd
      rel_path_target = op.relpath(real_path, prefix) # cctbx/bin/cmd
      parent_dirs = [ ".." for d in op.split(rel_path_link) ][1:]
      if (len(parent_dirs) > 0) :
        parent_dir = op.join(*parent_dirs)
      else :
        parent_dir = "."
      new_path = op.join(parent_dir, rel_path_target) # ../../cctbx/bin/cmd
      print >> log, "  creating symlink to %s" % new_path
      os.remove(file_name)
      os.symlink(new_path, file_name)

def find_and_delete_files (dir_name, file_name=None, file_ext=None) :
  """
  Recursively walk through a directory and delete any files (or directories)
  with the specified file name or extension.  Effectively equivalent to
  running 'find . -name "PATTERN" | xargs rm -rf'.
  """
  deleted = []
  for dirname, dirnames, filenames in os.walk(dir_name) :
    for dn in dirnames :
      if (dn == file_name) :
        full_path = op.join(dirname, dn)
        shutil.rmtree(full_path)
        deleted.append(full_path)
    for fn in filenames :
      full_path = op.join(dirname, fn)
      if (fn == file_name) :
        os.remove(full_path)
        deleted.append(full_path)
      elif (file_ext is not None) and fn.endswith(file_ext) :
        os.remove(full_path)
        deleted.append(full_path)
  return deleted

def archive_dist (dir_name, create_tarfile=True, use_shutil=True) :
  """
  Create a clean copy of a source repository, optionally bundling it as a
  gzipped tar file.
  """
  dir_name = op.abspath(dir_name)
  assert op.isdir(dir_name) and (dir_name != os.getcwd())
  module_name = op.basename(dir_name)
  local_path = op.join(os.getcwd(), module_name)
  if (use_shutil) :
    shutil.copytree(dir_name, local_path)
  else :
    copy_tree(dir_name, local_path)
  if op.exists(op.join(local_path, ".svn")) :
    try :
      call("svnversion %s > %s/.svnversion" % (module_name, module_name),
        log=sys.stdout)
    except RuntimeError, e :
      print e
  find_and_delete_files(local_path, file_ext=".pyc")
  find_and_delete_files(local_path, file_ext=".pyo")
  find_and_delete_files(local_path, file_ext=".swp")
  find_and_delete_files(local_path, file_name=".svn")
  find_and_delete_files(local_path, file_name=".git")
  find_and_delete_files(local_path, file_name=".sconsign")
  if (create_tarfile) :
    call("tar -czf %s.tar.gz %s" % (module_name, module_name), log=sys.stdout)
    tar_file = op.join(os.getcwd(), module_name + ".tar.gz")
    assert op.isfile(tar_file)
    shutil.rmtree(local_path)
    return tar_file
  return op.join(os.getcwd(), module_name)

def strip_libs (dir_name, log) :
  for dirname, dirnames, filenames in os.walk(dir_name) :
    for fn in filenames :
      full_path = op.join(dirname, fn)
      if fn.endswith(".so") and op.isfile(full_path) :
        try :
          call("strip %s" % full_path, log=log)
        except RuntimeError :
          pass

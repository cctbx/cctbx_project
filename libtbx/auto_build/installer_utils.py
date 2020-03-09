from __future__ import absolute_import, division, print_function, with_statement

import errno
import os
import os.path as op
import platform
import shutil
import stat
import subprocess
import sys
import tarfile
import time

# CCTBX itself requires Python 2.7, but this script is intended to bootstrap an
# installation on older systems as well
def check_python_version():
  if sys.hexversion < 0x2060000:
    sys.exit("Python version 2.6 or greater required to run this script")

def call(args, log=None, shell=True, cwd=None, verbose=False, env=None):
  if log is None: log = sys.stdout
  # shell=True requires string as args.
  if shell and isinstance(args, list):
    args = " ".join(args)
  if verbose:
    stdout = subprocess.PIPE
  else:
    stdout = log
  p = subprocess.Popen(
    args=args,
    shell=shell,
    cwd=cwd,
    bufsize=-1,
    stdin=None,
    stdout=stdout,
    stderr=subprocess.STDOUT,
    universal_newlines=True,
    close_fds=False,
    env=env)
  if verbose:
    while p.poll() is None:
      line = p.stdout.readline()
      # this will cause a deadlock if process is writing to stderr
      # stderr is redirected to stdout, so this cannot happen
      if line:
        print(": " + line.strip())
        log.write(line)
  #o, e = p.communicate()
  #log.write(o)
  log.flush()
  p.wait()
  log.flush()
  rc = p.returncode
  if rc != 0:
    raise RuntimeError("Call to '%s' failed with exit code %d" % (args, rc))

def check_output(*popenargs, **kwargs):
  # Back-port of Python 2.7 subprocess.check_output.
  try:
    process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
  except OSError as exc:
    if exc.errno == errno.ENOENT:
      raise OSError("No such file or directory (%s)" % popenargs[0])
    raise
  output, unused_err = process.communicate()
  retcode = process.poll()
  if retcode:
    raise RuntimeError("Call to '%s' failed with exit code %d" % (popenargs, retcode))
  return output

def untar(pkg_name, log=sys.stdout, verbose=False, change_ownership=False,
    check_output_path=True):
  assert os.path.isfile(pkg_name), pkg_name
  verbose_flag = owner_flag = ""
  if verbose:
    verbose_flag = "v"
  if change_ownership:
    owner_flag = "o"
  cmd = [ "tar", "x%s%sf" % (owner_flag, verbose_flag) ]
  if pkg_name.endswith("gz"):
    cmd = ["tar", "zx%s%sf" % (owner_flag, verbose_flag) ]
  elif pkg_name.endswith("bz2"):
    cmd = ["tar", "jx%s%sf" % (owner_flag, verbose_flag) ]
  args = cmd + [pkg_name]
  if os.name != 'nt':
    call(" ".join(args), log)
  else:
    # Note: This code breaks
    #   - extracting compressed files
    #   - logging
    #   - function parameters verbose and change_ownership
    #   - insufficient trapping of errors
    # Should import tar_extract from  bootstrap.py to avoid code duplication
    tar = tarfile.open(pkg_name)
    tar.extractall()
    tar.close()

  dir_name = os.path.basename(pkg_name) \
                 .replace(".tar.bz2", "") \
                 .replace(".tar.gz", "") \
                 .replace(".tgz", "") \
                 .replace(".tar", "")
  if check_output_path:
    if not os.path.isdir(dir_name):
      time.sleep(1)
      if not os.path.isdir(dir_name):
        if os.path.isdir(dir_name.capitalize()):
          dir_name = dir_name.capitalize()
        else:
          raise RuntimeError("Expected directory '%s' not found!" % dir_name)
    return os.path.abspath(dir_name)
  return None

def detect_osx_version():
  uname = os.uname()
  version = uname[2]
  major, minor, rev = version.split(".")
  return int(major)

def copy_file(src_path, dest_path, executable=None):
  assert os.path.isfile(src_path)
  with open(src_path, "rb") as fi:
    with open(dest_path, "wb") as fo:
      fo.write(fi.read())
  if os.access(src_path, os.X_OK) or executable:
    mode = os.stat(dest_path).st_mode
    os.chmod(dest_path, mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

# shutil.copytree replacement - handles symlinks intelligently and also
# circumvents this bug:
# http://bugs.python.org/issue14662
# this is not a general solution, but it is good enough for the installer
# process (at least on Unix).
def copy_tree(src_path, dest_path, verbose=False, log=sys.stdout):
  if sys.hexversion >= 0x2070400:
    return shutil.copytree(src_path, dest_path, symlinks=True)
  assert os.path.isdir(src_path), src_path
  assert not os.path.exists(dest_path), dest_path
  if verbose:
    print("creating %s" % dest_path, file=log)
  os.makedirs(dest_path)
  for path_name in os.listdir(src_path):
    node_src_path = os.path.join(src_path, path_name)
    node_dest_path = os.path.join(dest_path, path_name)
    if os.path.islink(node_src_path):
      target_path = os.readlink(node_src_path)
      os.symlink(target_path, node_dest_path)
    elif os.path.isfile(node_src_path):
      if verbose:
        print("  copy %s -> %s" % (node_src_path, node_dest_path), file=log)
      copy_file(node_src_path, node_dest_path)
    elif os.path.isdir(node_src_path):
      copy_tree(node_src_path, node_dest_path)
    elif verbose:
      print("  skipping %s" % node_src_path, file=log)

def get_os_version():
  uname = os.uname()
  kernel_version = uname[2]
  if uname[0] == "Darwin":
    # Use Python's platform module to determine OSX version
    release, _, _ = platform.mac_ver()
    # Convert X.Y.Z to X.Y - discard minor version
    release_major = ".".join(release.split(".")[:2])
    assert release, "Could not determine OSX version"
    return release_major
  return kernel_version

def machine_type():
  if sys.platform == "win32":
    mtype = "intel-windows"
    arch, os_type = platform.architecture()
    if arch == "64bit":
      mtype += "-x86_64"
    return mtype
  uname = os.uname()
  if uname[0] == "Linux":
    mtype = "intel-linux-2.6"
  elif uname[0] == "Darwin":
    mtype = "mac-intel-osx"
  else:
    mtype = uname[0]
  if uname[-1] == "x86_64":
    mtype += "-x86_64"
  elif mtype == "mac-intel-osx":
    version_fields = uname[2].split(".")
    major_version = int(version_fields[0])
    if major_version >= 10:
      mtype += "-x86_64"
    if major_version >= 13:
      mtype += "-10.9"
  return mtype

def regenerate_relative_symlinks(dir_name, log=sys.stdout):
  """
  Rewrite all symlinks in a directory with relative paths (allows later
  relocation).  This is mostly just for Mac OS where there are a lot of links
  into the Python.framework bundle.
  """
  old_cwd = os.getcwd()
  os.chdir(dir_name)
  for file_name in os.listdir(dir_name):
    # e.g. /usr/bin/cmd
    full_path = op.join(dir_name, file_name)
    if op.islink(full_path):
      real_path = op.realpath(full_path) # e.g. /usr/cctbx/bin/cmd
      prefix = op.commonprefix([full_path, real_path]) # /usr
      rel_path_link = op.relpath(full_path, prefix) # bin/cmd
      rel_path_target = op.relpath(real_path, prefix) # cctbx/bin/cmd
      parent_dirs = [ ".." for d in op.split(rel_path_link) ][1:]
      if parent_dirs:
        parent_dir = op.join(*parent_dirs)
      else:
        parent_dir = "."
      new_path = op.join(parent_dir, rel_path_target) # ../../cctbx/bin/cmd
      print("  creating symlink to %s" % new_path, file=log)
      os.remove(file_name)
      os.symlink(new_path, file_name)

def find_and_delete_files(dir_name, file_name=None, file_ext=None):
  """
  Recursively walk through a directory and delete any files (or directories)
  with the specified file name or extension.  Effectively equivalent to
  running 'find . -name "PATTERN" | xargs rm -rf'.
  """
  deleted = []
  for dirname, dirnames, filenames in os.walk(dir_name):
    for dn in dirnames:
      if dn == file_name:
        full_path = op.join(dirname, dn)
        shutil.rmtree(full_path)
        deleted.append(full_path)
    for fn in filenames:
      full_path = op.join(dirname, fn)
      if fn == file_name:
        os.remove(full_path)
        deleted.append(full_path)
      elif file_ext is not None and fn.endswith(file_ext):
        os.remove(full_path)
        deleted.append(full_path)
  return deleted

def archive_dist(dir_name, create_tarfile=True, use_shutil=True):
  """
  Create a clean copy of a source repository, optionally bundling it as a
  gzipped tar file.
  """
  dir_name = op.abspath(dir_name)
  assert op.isdir(dir_name) and dir_name != os.getcwd()
  module_name = op.basename(dir_name)
  local_path = op.join(os.getcwd(), module_name)
  if use_shutil:
    shutil.copytree(dir_name, local_path)
  else:
    copy_tree(dir_name, local_path)
  if op.exists(op.join(local_path, ".svn")):
    try:
      call("svnversion %s > %s/.svnversion" % (module_name, module_name),
        log=sys.stdout)
    except RuntimeError as e:
      print(e)
  find_and_delete_files(local_path, file_ext=".pyc")
  find_and_delete_files(local_path, file_ext=".pyo")
  find_and_delete_files(local_path, file_ext=".swp")
  find_and_delete_files(local_path, file_name=".svn")
  find_and_delete_files(local_path, file_name=".git")
  find_and_delete_files(local_path, file_name=".sconsign")
  if create_tarfile:
    #call("tar -czf %s.tar.gz %s" % (module_name, module_name), log=sys.stdout)
    tar = tarfile.open("%s.tar.gz" %module_name, "w:gz")
    tar.add(module_name)
    tar.close()

    tar_file = op.join(os.getcwd(), module_name + ".tar.gz")
    assert op.isfile(tar_file)
    shutil.rmtree(local_path)
    return tar_file
  return op.join(os.getcwd(), module_name)

def strip_libs(dir_name, log):
  for dirname, dirnames, filenames in os.walk(dir_name):
    for fn in filenames:
      full_path = op.join(dirname, fn)
      if fn.endswith(".so") and op.isfile(full_path):
        try:
          call("strip %s" % full_path, log=log)
        except RuntimeError:
          pass

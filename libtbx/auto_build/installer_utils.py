
from __future__ import division
import warnings
import time
import stat
import re
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

def call (args, log) :
  rc = None
  if (sys.version_info[1] >= 7) :
    import subprocess
  else :
    # XXX HACK
    libtbx_path = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
    if (not libtbx_path in sys.path) :
      sys.path.append(libtbx_path)
    import subprocess_with_fixes as subprocess
  if isinstance(args, list) :
    args = " ".join(args)
  p = subprocess.Popen(
    args=args,
    shell=True,
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
    os.chmod(dest_path, mode | stat.S_IXUSR)

# shutil.copytree replacement - circumvents this bug:
# http://bugs.python.org/issue14662
# this is not a general solution, but it is good enough for the installer
# process (at least on Unix)
def copy_tree (src_path, dest_path, verbose=False, log=sys.stdout) :
  assert os.path.isdir(src_path) and not os.path.exists(dest_path)
  if (verbose) :
    print >> log, "creating %s" % dest_path
  os.makedirs(dest_path)
  for path_name in os.listdir(src_path) :
    node_src_path = os.path.join(src_path, path_name)
    node_dest_path = os.path.join(dest_path, path_name)
    if os.path.isfile(node_src_path) :
      if (verbose) :
        print >> log, "  copy %s -> %s" % (node_src_path, node_dest_path)
      copy_file(node_src_path, node_dest_path)
    elif os.path.isdir(node_src_path) :
      copy_tree(node_src_path, node_dest_path)
    else :
      if (verbose) :
        print >> log, "  skipping %s" % node_src_path

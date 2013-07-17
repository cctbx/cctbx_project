
from __future__ import division
import warnings
import time
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

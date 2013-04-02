
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

def call (args, log, join_stdout_stderr=True) :
  rc = None
  if (sys.version_info[1] >= 7) :
    import subprocess
  else :
    import subprocess_with_fixes as subprocess
  stderr = subprocess.PIPE
  if (join_stdout_stderr) :
    stderr = subprocess.STDOUT
  p = subprocess.Popen(
    args=args,
    shell=True,
    stdout=subprocess.PIPE,
    stderr=stderr)
  o, e = p.communicate()
  log.write(o)
  rc = p.returncode
  if (rc != 0) :
    raise RuntimeError("Call to '%s' failed with exit code %d" % (args, rc))

def untar (pkg_name, log=sys.stdout, verbose=False) :
  assert os.path.isfile(pkg_name), pkg_name
  verbose_flag = ""
  if (verbose) :
    verbose_flag = "v"
  cmd = "tar x%sf" % verbose_flag
  if (pkg_name.endswith("gz")) :
    cmd = "tar zx%sf" % verbose_flag
  elif (pkg_name.endswith("bz2")) :
    cmd = "tar jx%sf" % verbose_flag
  args = "%s %s" % (cmd, pkg_name)
  call(args, log)
  dir_name = re.sub(".tgz", "",
               re.sub(".tar.gz", "",
                 re.sub(".tar.bz2", "",
                   os.path.basename(pkg_name))))
  if (not os.path.isdir(dir_name)) :
    time.sleep(1)
    if (not os.path.isdir(dir_name)) :
      raise RuntimeError("Expected directory '%s' not found!" % dir_name)
  return os.path.abspath(dir_name)

def detect_osx_version () :
  uname = os.uname()
  version = uname[2]
  major, minor, rev = version.split(".")
  return int(major)

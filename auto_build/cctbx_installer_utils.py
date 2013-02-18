
from __future__ import division
import re
import os
import sys

def call (args, log, join_stdout_stderr=True) :
  rc = None
  try :
    import subprocess
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
  except ImportError : # TODO needs testing
    child_stdin, child_stdout, child_stderr = os.popen3(args, "t")
    child_stdin.close()
    log.write(child_stdout)
    rc = 0
  if (rc != 0) :
    raise RuntimeError("Call to '%s' failed with exit code %d" % (args, rc))

def untar (pkg_name, log=sys.stdout, verbose=False) :
  assert os.path.isfile(pkg_name)
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
  assert os.path.isdir(dir_name)
  return os.path.abspath(dir_name)

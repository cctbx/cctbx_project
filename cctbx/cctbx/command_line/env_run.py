#! /usr/bin/env python

import sys, os, os.path
assert len(sys.argv) >= 3, "Insufficient command line arguments."
try:
  cmd_root = os.environ[sys.argv[1]]
except:
  print >> sys.stderr, 'FATAL ERROR: Environment variable "%s" not defined.' % (
    sys.argv[1])
  sys.exit(1)
cmd = "%s %s %s" % (
  sys.executable,
  os.path.join(cmd_root, sys.argv[2]),
  " ".join(sys.argv[3:]))
os.system(cmd)

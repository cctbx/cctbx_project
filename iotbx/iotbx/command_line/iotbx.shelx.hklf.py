#! /bin/env python
from iotbx.shelx import hklf
from scitbx.python_utils import easy_pickle
import sys

def run():
  if (len(sys.argv) != 3):
    print >> sys.stderr, "usage: hklf [read|write] file"
    sys.exit(1)
  assert sys.argv[1] in ("write", "read")
  if (sys.argv[1] == "read"):
    easy_pickle.dump(sys.argv[2]+".pickle", hklf.read(open(sys.argv[2], "r")))
  else:
    hklf.write(miller_array=easy_pickle.load(sys.argv[2]))

if (__name__ == "__main__"):
  run()

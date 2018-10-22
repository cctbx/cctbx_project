from __future__ import division
import compileall
from six.moves import cStringIO as StringIO
import sys, os

def run():
  dirs = sys.argv[1:]
  if (len(dirs) == 0):
    dirs = [os.getcwd()]
  sys.stdout = StringIO()
  sys.stderr = StringIO()
  for dir in dirs:
    compileall.compile_dir(dir, 100)

if (__name__ == "__main__"):
  run()

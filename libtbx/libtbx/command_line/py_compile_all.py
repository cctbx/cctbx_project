import compileall
import cStringIO
import sys

def run():
  dirs = sys.argv[1:]
  if (len(dirs) == 0):
    dirs = ["."]
  sys.stdout = cStringIO.StringIO()
  sys.stderr = cStringIO.StringIO()
  for dir in dirs:
    compileall.compile_dir(dir, 100)

if (__name__ == "__main__"):
  run()

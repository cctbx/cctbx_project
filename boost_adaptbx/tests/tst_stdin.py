import boost.python_file
ext = boost.python.import_ext("python_file_test_ext")
import libtbx.load_env
from libtbx import easy_run
try: import subprocess_with_fixes as subprocess
except ImportError: import subprocess
import sys, os

def exercise():
  proc = subprocess.Popen(args='libtbx.python %s --core' % __file__,
                          shell=True,
                          stdin=subprocess.PIPE,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
  output, error = proc.communicate("Veni Vidi Vici")
  assert not error, error
  assert not output, output

def read_from_stdin():
  written = ext.test_read(sys.stdin, "read")
  assert written == "Veni, Vidi, Vici, [ fail, eof ]", written

def run(core):
  if not core:
    exercise()
    print 'OK'
  else:
    read_from_stdin()

if __name__ == '__main__':
  run(core='--core' in sys.argv[1:])

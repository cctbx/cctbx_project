import boost.python
from boost.python import streambuf
ext = boost.python.import_ext("boost_adaptbx_python_streambuf_test_ext")
try: from libtbx import subprocess_with_fixes as subprocess
except ImportError: import subprocess
import sys

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
  written = ext.test_read(streambuf(sys.stdin), "read")
  assert written == "Veni, Vidi, Vici, [ fail, eof ]", written

def run(core):
  if not core:
    exercise()
    print 'OK'
  else:
    read_from_stdin()

if __name__ == '__main__':
  run(core='--core' in sys.argv[1:])

from __future__ import absolute_import, division, print_function

import boost_adaptbx.boost.python as bp
from boost_adaptbx.boost.python import ostream
ext = bp.import_ext("boost_adaptbx_python_streambuf_test_ext")
import subprocess
import sys

def exercise():
  assert __file__.find('"') < 0
  proc = subprocess.Popen(args='libtbx.python "%s" --core' % __file__,
                          shell=True,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
  output, error = proc.communicate()
  assert not error, error
  assert output == b"2 times 1.6 equals 3.2", output

def write_to_stdout():
  ext.test_write(ostream(sys.stdout), "write")

def run(core):
  if not core:
    exercise()
    print('OK')
  else:
    write_to_stdout()

if __name__ == '__main__':
  run(core='--core' in sys.argv[1:])

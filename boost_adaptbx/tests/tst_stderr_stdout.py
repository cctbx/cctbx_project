import boost.python_file
import sys
import gc

def run():
  ext = boost.python.import_ext("python_file_test_ext")
  ext.call_with_stderr_stdout_do_nothing(sys.stderr, sys.stdout)
  gc.collect()
  print "OK"

if (__name__ == "__main__"):
  run()

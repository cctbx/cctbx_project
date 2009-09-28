import libtbx.object_oriented_patterns as oop
import boost.python_file
import sys
import gc

class without_tell(oop.proxy):
  """ sys.stdout and sys.stderr don't have a functional method 'tell',
  except on MacOS X. This proxy brings the latter back into the ranks """

  def tell(self):
    raise NotImplementedError("Test of stdout/stderr / C++ stream bridge")


def run():
  ext = boost.python.import_ext("python_file_test_ext")
  ext.call_with_stderr_stdout_do_nothing(
    without_tell(sys.stderr), # bug trigger on MacOS X
    sys.stdout)
  gc.collect()
  print "OK"

if (__name__ == "__main__"):
  run()

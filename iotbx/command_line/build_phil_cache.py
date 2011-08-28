
import libtbx.load_env
from libtbx import phil
from libtbx.utils import import_python_object
import os
import sys

def run (args) :
  phil_path = args[-1]
  phil_inp = import_python_object(
    import_path=phil_path,
    error_prefix="",
    target_must_be="",
    where_str="").object
  if (isinstance(phil_inp, str)) :
    phil_object = phil.parse(phil_inp)
  elif (hasattr(phil_inp, "__call__")) :
    phil_object = phil_inp()
  else :
    print type(phil_inp)
    assert isinstance(phil_inp, phil.scope)
    phil_object = phil_inp
  cache_dir = os.path.join(libtbx.env.build_path, "phil_cache")
  if (not os.path.isdir(cache_dir)) :
    os.mkdir(cache_dir)
  full_path = os.path.join(cache_dir, phil_path)
  f = open("%s.phil" % full_path, "w")
  phil_object.show(out=f, attributes_level=3)
  f.close()
  print "Wrote master parameters to %s.phil" % full_path

if (__name__ == "__main__") :
  run(sys.argv[1:])

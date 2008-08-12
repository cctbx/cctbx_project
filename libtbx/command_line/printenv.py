import libtbx.load_env
import os

def run():
  var_names = os.environ.keys()
  var_names.sort()
  for var_name in var_names:
    print "%s=%s" % (var_name, os.environ[var_name])

if (__name__ == "__main__"):
  run()

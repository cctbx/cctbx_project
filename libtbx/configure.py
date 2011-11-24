import sys, os

def run():
  if (not hasattr(sys, "version_info")
      or sys.version_info[0] < 2
      or (sys.version_info[0] == 2 and sys.version_info[1] < 3)):
    print
    print "*" * 78
    print "FATAL: Python 2.3 or higher is required."
    print "Version currently in use:", sys.version
    print "*" * 78
    print
    return
  sys.path.insert(1, os.path.join(sys.path[0], "pythonpath"))
  sys.path[0] = os.path.dirname(sys.path[0])
  import libtbx.env_config
  libtbx.env_config.cold_start(sys.argv)
  print "Done."

if (__name__ == "__main__"):
  run()

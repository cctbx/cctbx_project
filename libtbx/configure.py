import sys, os

def run():
  if (not hasattr(sys, "version_info")
      or sys.version_info[0] < 2
      or (sys.version_info[0] == 2 and sys.version_info[1] < 2)):
    print
    print "*" * 78
    print "FATAL: Python 2.2 or higher is required."
    print "Version currently in use:", sys.version
    print "*" * 78
    print
    return
  if (os.name == "nt"):
    open("shortpath.bat", "w").write("@echo off\necho %~s1\n")
  from os.path import normpath, join, dirname
  libtbx_libtbx = normpath(join(dirname(sys.argv[0]), "libtbx"))
  sys.path.insert(0, libtbx_libtbx)
  from libtbx.command_line import configure
  configure.cold_start(sys.argv)

if (__name__ == "__main__"):
  run()

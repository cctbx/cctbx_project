from os.path import normpath, join, dirname
import sys

def run():
  libtbx_libtbx = normpath(join(dirname(sys.argv[0]), "libtbx"))
  sys.path.insert(0, libtbx_libtbx)
  from libtbx.command_line import configure
  configure.cold_start(sys.argv)

if (__name__ == "__main__"):
  run()

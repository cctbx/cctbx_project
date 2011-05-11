import libtbx.env_config
import sys

if (__name__ == "__main__"):
  libtbx.env_config.warm_start(sys.argv)
  print "Done."

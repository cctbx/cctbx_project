from __future__ import absolute_import, division, print_function
import sys, os

def run():
  if sys.hexversion < 0x02070000:
    print()
    print("*" * 78)
    print("FATAL: Python 2.7 or higher is required.")
    print("Version currently in use:", sys.version)
    print("*" * 78)
    return False
  sys.path[0] = os.path.dirname(sys.path[0])
  import libtbx.env_config
  libtbx.env_config.cold_start(sys.argv)
  print("Done.")
  return True

if __name__ == "__main__":
  if not run():
    sys.exit(1)

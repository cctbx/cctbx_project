from __future__ import absolute_import, division, print_function
import sys, os

def run():
  if sys.hexversion < 0x02070000:
    print()
    print("*" * 78)
    print("FATAL: Python 2.7 or higher is required.")
    print("Version currently in use:", sys.version)
    print("*" * 78)
    print()
    return False
  # test for six and future
  try:
    import six
    import future
  except ImportError:
    print()
    print("*" * 78)
    print("FATAL: For Python 2/3 compatibility, CCTBX currently requires the six and\n  future modules.")
    print("To install, please use pip or conda (if available)")
    print("  pip install six future")
    print("  conda install six future")
    print("*" * 78)
    print()
    return False
  sys.path[0] = os.path.dirname(sys.path[0])
  import libtbx.env_config
  libtbx.env_config.cold_start(sys.argv)
  print("Done.")
  return True

if __name__ == "__main__":
  if not run():
    sys.exit(1)

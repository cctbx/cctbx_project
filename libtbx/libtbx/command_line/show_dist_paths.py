import libtbx.env
import sys
if (len(sys.argv) == 1):
  for path in libtbx.env.cache.dist_paths.values():
    print path
else:
  for arg in sys.argv[1:]:
    print libtbx.env.cache.dist_path(arg, None)

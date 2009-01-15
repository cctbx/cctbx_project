# LIBTBX_SET_DISPATCHER_NAME phenix.grow_density

from mmtbx import grow_density
import sys

if(__name__ == "__main__"):
  grow_density.run(sys.argv[1:], show_geometry_statistics = False)

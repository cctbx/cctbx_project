# LIBTBX_SET_DISPATCHER_NAME mmtbx.afitt

from __future__ import division
import sys
from mmtbx.geometry_restraints import afitt

if __name__ == "__main__":
  afitt.run(args=sys.argv[1:])

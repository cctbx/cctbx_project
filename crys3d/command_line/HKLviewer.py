# LIBTBX_SET_DISPATCHER_NAME phenix.HKLviewer
# LIBTBX_SET_DISPATCHER_NAME cctbx.HKLviewer
# LIBTBX_SET_DISPATCHER_NAME phasertng.HKLviewer
from __future__ import absolute_import, division, print_function

from crys3d.hklview import HKLviewer


if (__name__ == "__main__") :
  HKLviewer.run(isembedded=True)

# LIBTBX_SET_DISPATCHER_NAME cxi.monitor_detectors
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from xfel.cxi.gfx import wx_detectors
import sys

if (__name__ == "__main__") :
  wx_detectors.run(sys.argv[1:])

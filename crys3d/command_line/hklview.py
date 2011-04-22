# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from crys3d.wx_hklview import *
import sys

if (__name__ == "__main__") :
  from iotbx import file_reader
  f = file_reader.any_file(sys.argv[-1])
  ma = None
  for array in f.file_server.miller_arrays :
    if array.is_xray_amplitude_array() or array.is_xray_intensity_array() :
      ma = array
      break
  a = wx.App(0)
  f = HKLViewFrame(None, -1, "HKL viewer", size=(1024,768))
  f.set_miller_array(ma)
  f.Show()
  a.MainLoop()

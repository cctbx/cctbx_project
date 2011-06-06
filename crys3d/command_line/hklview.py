# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from crys3d.hklview.frames import *
import sys

def run (args) :
  from iotbx import file_reader
  ma = None
  if (len(args) == 0) :
    from cctbx import miller, crystal
    from cctbx.array_family import flex
    xs = crystal.symmetry((3,3,5,90,90,120), "P6")
    mi = flex.miller_index([
      (0,0,1),(0,0,2),(0,0,3),
      (0,1,0),(0,2,0),(0,3,0),
      (1,1,0),(1,2,0),(1,3,0)])
    d = flex.double([1.0, 2.0, 3.0, 5.0, 10.0, 15.0, 6.0, 9.0, 12.0])
    s = miller.set(xs, mi, anomalous_flag=False)
    ma = s.array(data=d).set_info("test")
  else :
    f = file_reader.any_file(args[-1])
    for array in f.file_server.miller_arrays :
      if array.is_xray_amplitude_array() or array.is_xray_intensity_array() :
        ma = array
        break
  a = wx.App(0)
  f = HKLViewFrame(None, -1, "HKL viewer", size=(1024,768))
  f.set_miller_array(ma)
  f.Show()
  a.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])

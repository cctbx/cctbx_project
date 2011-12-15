# LIBTBX_SET_DISPATCHER_NAME phenix.data_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import crys3d.hklview
from crys3d.hklview.frames import *
from cctbx.miller.display import master_phil
from wxtbx import icons
import iotbx.phil
import wx
import sys

def run (args) :
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
    settings = crys3d.hklview.settings()
  else :
    pcl = iotbx.phil.process_command_line_with_files(
      args=args,
      master_phil=master_phil,
      reflection_file_def="data")
    settings = pcl.work.extract()
  a = wx.App(0)
  app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_3d.GetBitmap())
  if (wx.VERSION >= (2,9)) :
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
  else :
    tb_icon = wx.TaskBarIcon()
  tb_icon.SetIcon(app_icon, "PHENIX data viewer")
  a.hklview_settings = settings
  f = HKLViewFrame(None, -1, "Reflection data viewer", size=(1024,768))
  f.Show()
  if (ma is not None) :
    f.set_miller_array(ma)
  elif (settings.data is not None) :
    f.load_reflections_file(settings.data)
  else :
    f.OnLoadFile(None)
  a.SetTopWindow(f)
  a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), f)
  a.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])

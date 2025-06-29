"""View pseudo-precession planes through a dataset"""
# LIBTBX_SET_DISPATCHER_NAME phenix.data_viewer
# LIBTBX_SET_DISPATCHER_NAME cctbx.data_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function
from crys3d.hklview.frames import *
from cctbx.miller.display import master_phil
from wxtbx import icons
import wxtbx.app
import iotbx.phil
import wx
import sys

def run (args) :
  ma = None
  show_2d = False
  if ("--2d" in args) :
    show_2d = True
    args.remove("--2d")
  pcl = iotbx.phil.process_command_line_with_files(
    args=args,
    master_phil=master_phil,
    reflection_file_def="data",
    pdb_file_def="symmetry_file",
    usage_string="phenix.data_viewer f_obs.mtz [options]")
  settings = pcl.work.extract()
  a = wxtbx.app.CCTBXApp(0)
  app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_3d.GetBitmap())
  if (wx.VERSION >= (2,9)) :
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
  else :
    tb_icon = wx.TaskBarIcon()
  tb_icon.SetIcon(app_icon, "PHENIX data viewer")
  a.hklview_settings = settings
  viewer_class = HKLViewFrame
  if (show_2d) :
    viewer_class = HKLViewFrame2D
  f = viewer_class(None, -1, "Reflection data viewer", size=(1024,768))
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


"""Similar to phenix.data_viewer, as side-by-side view"""
from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME phenix.compare_datasets
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from crys3d.hklview import comparison
from wxtbx import icons
import iotbx.phil
from libtbx.utils import Sorry, Usage
import wx
import os
import sys

# TODO
master_phil = iotbx.phil.parse("""
file_name_1 = None
  .type = path
label_1 = None
  .type = str
file_name_2 = None
  .type = path
label_2 = None
  .type = str
include scope cctbx.miller.display.master_phil
""", process_includes=True)

def run (args) :
  if ("--help" in args) or ("--options" in args) :
    raise Usage("""\
phenix.compare_datasets data1.mtz data2.mtz

Side-by-side visualization of pseudo-precession planes through reciprocal
space for a pair of datasets - essentially a duplex version of 2D view in
phenix.data_viewer.
""")
  a = wx.App(0)
  app_icon = wx.EmptyIcon()
  app_icon.CopyFromBitmap(icons.hklview_2d.GetBitmap())
  if (wx.VERSION >= (2,9)) :
    tb_icon = wx.TaskBarIcon(wx.TBI_DOCK)
  else :
    tb_icon = wx.TaskBarIcon()
  tb_icon.SetIcon(app_icon, "PHENIX data viewer")
  frame = comparison.ComparisonFrame(None, -1, "Dataset comparison")
  frame.Show()
  file_name_1 = file_name_2 = None
  for arg in args :
    if (os.path.isfile(arg)) :
      if (file_name_1 is None) :
        file_name_1 = arg
      elif (file_name_2 is None) :
        file_name_2 = arg
      else :
        raise Sorry("Only two files at a time are supported.")
  if (file_name_1 is not None) and (file_name_2 is not None) :
    frame.load_files(file_name_1, file_name_2)
  a.SetTopWindow(frame)
  a.Bind(wx.EVT_WINDOW_DESTROY, lambda evt: tb_icon.Destroy(), frame)
  a.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])


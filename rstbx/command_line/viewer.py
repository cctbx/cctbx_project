# LIBTBX_SET_DISPATCHER_NAME phenix.image_viewer
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from rstbx.viewer.frame import XrayFrame
import wx
import os
import sys

def run (args) :
  app = wx.App(0)
  frame = XrayFrame(None, -1, "X-ray image display", size=(800,720))
  if (len(args) == 1 and os.path.basename(args[0]) == "DISTL_pickle") :
    assert os.path.isfile(args[0])
    frame.load_distl_output(args[0])
  else :
    for file_name in args:
      assert os.path.isfile(file_name)
      frame.add_file_name(file_name)
    frame.load_image(args[0])
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])

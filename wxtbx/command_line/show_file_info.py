# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from wxtbx import info_panels
from libtbx.utils import Usage, Sorry
import wx
import os
import sys

def run (args) :
  if (len(args) == 0) or (not os.path.isfile(args[-1])):
    raise Usage("wxtbx.show_file_info [options] file_name")
  app = wx.App(0)
  from iotbx import file_reader
  f = file_reader.any_file(args[-1],
    raise_sorry_if_errors=True)
  f.show_summary()
  frame = None
  if (f.file_type == "pdb") :
    frame = info_panels.PDBFileInfo(None)
  elif (f.file_type == "hkl") :
    frame = info_panels.ReflectionFileInfo(None)
  elif (f.file_type == "img") :
    frame = info_panels.ImageFileInfo(None)
  if (frame is None) :
    raise Sorry("File type %s not supported for display." % f.file_type)
  frame.set_file(f.file_name)
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])

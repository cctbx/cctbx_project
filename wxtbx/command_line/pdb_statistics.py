from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME mmtbx.pdb_statistics
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from wxtbx.polygon_db_viewer import ConfigFrame
import wx

if (__name__ == "__main__"):
  app = wx.App(0)
  frame = ConfigFrame(None, -1, "PDB statistics from POLYGON database")
  frame.Show()
  app.MainLoop()

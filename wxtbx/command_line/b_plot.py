# LIBTBX_SET_DISPATCHER_NAME mmtbx.bfactor_plot
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT

from wxtbx import b_plot
import wx
from libtbx.utils import Usage
import sys

def run (args) :
  if (len(args) == 0) :
    raise Usage("mmtbx.bfactor_plot model.pdb")
  pdb_file = args[0]
  params = b_plot.master_phil.fetch().extract()
  params.b_plot.pdb_file = pdb_file
  result = b_plot.run(params=params)
  app = wx.App(0)
  b_plot.show_plot_frame(result)
  app.MainLoop()

if (__name__ == "__main__") :
  run(sys.argv[1:])

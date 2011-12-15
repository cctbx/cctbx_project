# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from wxtbx import plots
import wx
import sys

def run (args) :
  logfile = args[0]
  app = wx.App(0)
  frame = plots.loggraph(parent=None,
    title="Loggraph test",
    tables=None,
    processed_lines=open(logfile).readlines())
  frame.Show()
  app.MainLoop()

if __name__ == "__main__" :
  run(sys.argv[1:])

from __future__ import absolute_import, division, print_function
# LIBTBX_SET_DISPATCHER_NAME cctbx.plot_anomalous_scattering
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from wxtbx import anomalous_scattering
import wxtbx.app
import sys

def run(args):
  app = wxtbx.app.CCTBXApp(0)
  frame = anomalous_scattering.AnomPlotFrame(parent=None,
    title="Anomalous scattering plot")
  frame.Show()
  app.MainLoop()

if (__name__ == "__main__"):
  run(sys.argv[1:])

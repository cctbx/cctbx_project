# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_SET_DISPATCHER_NAME wxtbx.xtriage
# LIBTBX_SET_DISPATCHER_NAME mmtbx.xtriage_gui

from __future__ import absolute_import, division, print_function
import wxtbx.xtriage
import wxtbx.app
import sys

if (__name__ == "__main__"):
  from mmtbx.scaling import xtriage
  result = xtriage.run(args=sys.argv[1:],
    data_file_name="xtriage_data.pkl")
  app = wxtbx.app.CCTBXApp(0)
  frame = wxtbx.xtriage.XtriageFrame(parent=None,
    title="Xtriage",
    size=(900,600))
  frame.SetResult(result)
  frame.Show()
  app.MainLoop()

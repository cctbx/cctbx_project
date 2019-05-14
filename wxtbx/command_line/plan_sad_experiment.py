# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

from __future__ import absolute_import, division, print_function
import wxtbx.xtriage
import wxtbx.app
import sys

if (__name__ == "__main__"):
  from mmtbx.command_line import plan_sad_experiment
  result = plan_sad_experiment.run(args=sys.argv[1:])
  app = wxtbx.app.CCTBXApp(0)
  frame = wxtbx.xtriage.XtriageFrameSingleResult(
    parent=None,
    title="Plan SAD experiment",
    size=(900,600))
  frame.SetResult(result)
  frame.Show()
  app.MainLoop()

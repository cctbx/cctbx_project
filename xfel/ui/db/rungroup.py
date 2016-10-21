from __future__ import division
from xfel.ui.db import db_proxy

class Rungroup(db_proxy):
  def __init__(self, app, rungroup_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_rungroup" % app.params.experiment_tag, id = rungroup_id, **kwargs)
    self.rungroup_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "runs":
      startrun = self.app.get_run(run_number = self.startrun).run
      if self.endrun is None:
        return [r for r in self.app.get_all_runs() if r.run >= startrun]
      else:
        endrun = self.app.get_run(run_number = self.endrun).run
        return [r for r in self.app.get_all_runs() if r.run >= startrun and \
                                                      r.run <= endrun]
    else:
      return super(Rungroup, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name != "runs"
    super(Rungroup, self).__setattr__(name, value)

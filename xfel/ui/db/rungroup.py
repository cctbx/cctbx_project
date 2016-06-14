from __future__ import division
from xfel.ui.db import db_proxy

class Rungroup(db_proxy):
  def __init__(self, app, rungroup_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_rungroup" % app.params.experiment_tag, id = rungroup_id, **kwargs)
    self.rungroup_id = self.id

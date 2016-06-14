from __future__ import division
from xfel.ui.db import db_proxy

class Run(db_proxy):
  def __init__(self, app, run_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_run" % app.params.experiment_tag, id = run_id, **kwargs)
    self.run_id = self.id


from __future__ import division
from xfel.ui.db import db_proxy

class Rungroup(db_proxy):
  def __init__(self, app, rungroup_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_rungroup" % app.params.experiment_tag, id = rungroup_id, **kwargs)

    if rungroup_id is None:
      rungroup_id = 4 # create a new rungroup

    self.rungroup_id = rungroup_id

    # dummy values:
    startrun = 1
    endrun = 5
    detz_parameter = 580
    binning = None
    comment = ""

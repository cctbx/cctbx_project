from __future__ import division
from xfel.ui.db.tag import Tag
from xfel.ui.db import db_proxy

class Run(db_proxy):
  def __init__(self, dbobj, run_id = None, **kwargs):
    db_proxy.__init__(self, dbobj, id = run_id, **kwargs)

    if run_id is None:
      run_id = 4 # create a new run

    self.run_id = run_id

    # dummy values:
    self.run = 42

    self.tags = [Tag(self.dbobj) for i in xrange(2)]

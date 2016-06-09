from __future__ import division
from xfel.ui.db.tag import Tag

class Run(object):
  def __init__(self, run_id = None, **kwargs):
    if run_id is None:
      run_id = 4 # create a new run

    self.run_id = run_id

    # dummy values:
    self.run = 42

    self.tags = [Tag() for i in xrange(2)]

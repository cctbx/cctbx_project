from __future__ import division
from xfel.ui.db import db_proxy

class Tag(db_proxy):
  def __init__(self, dbobj, tag_id = None, **kwargs):
    db_proxy.__init__(self, dbobj, id = tag_id, **kwargs)

    if tag_id is None:
      tag_id = 4 # create a new tag

    self.tag_id = tag_id

    # dummy values:
    self.name = "default tag"
    self.comment = ""

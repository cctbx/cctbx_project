from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy

class Tag(db_proxy):
  def __init__(self, app, tag_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_tag" % app.params.experiment_tag, id = tag_id, **kwargs)
    self.tag_id = self.id


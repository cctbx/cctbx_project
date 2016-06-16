from __future__ import division
from xfel.ui.db import db_proxy

class Run(db_proxy):
  def __init__(self, app, run_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_run" % app.params.experiment_tag, id = run_id, **kwargs)
    self.run_id = self.id

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "tags":
      return self.app.get_run_tags(self.id)
    else:
      return super(Run, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name != "tags"
    super(Run, self).__setattr__(name, value)

  def add_tag(self, tag):
    query = "INSERT INTO `%s_run_tag` (run_id, tag_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

  def remove_tag(self, tag):
    query = "DELETE FROM `%s_run_tag` WHERE run_id = %d AND tag_id = %s" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

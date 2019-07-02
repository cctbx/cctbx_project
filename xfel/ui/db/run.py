from __future__ import absolute_import, division, print_function
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

  def __str__(self):
    try:
      return str(int(self.run))
    except ValueError:
      return "%d: %s"%(self.id, self.run)

  def add_tag(self, tag):
    query = "INSERT INTO `%s_run_tag` (run_id, tag_id) VALUES (%d, %d)" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

  def remove_tag(self, tag):
    query = "DELETE FROM `%s_run_tag` WHERE run_id = %d AND tag_id = %s" % (
      self.app.params.experiment_tag, self.id, tag.id)
    self.app.execute_query(query, commit=True)

  def get_rungroups(self):
   from xfel.ui.db.rungroup import Rungroup
   tag = self.app.params.experiment_tag
   query = """SELECT rg.id FROM `%s_rungroup` rg
              JOIN `%s_rungroup_run` rgr on rg.id = rgr.rungroup_id
              WHERE rgr.run_id = %d AND rg.active=True
              """ %(tag, tag, self.id)
   cursor = self.app.execute_query(query)
   rungroup_ids = ["%d"%i[0] for i in cursor.fetchall()]
   if len(rungroup_ids) == 0:
     return []
   return self.app.get_all_x(Rungroup, "rungroup", where = "WHERE rungroup.id IN (%s)"%",".join(rungroup_ids))

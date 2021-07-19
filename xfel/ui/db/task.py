from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy
from xfel.ui import load_phil_scope_from_dispatcher

task_types = ["indexing", "ensemble_refinement", "scaling", "merging", "phenix"]
task_dispatchers = [None, "cctbx.xfel.stripe_experiment", "cctbx.xfel.merge", "cctbx.xfel.merge", None]
task_scope = ["local", "local", "local", "global", "global"]

class Task(db_proxy):
  def __init__(self, app, task_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_task" % app.params.experiment_tag, id = task_id, **kwargs)
    self.task_id = self.id
    self._trial = None

  def __getattr__(self, name):
    # Called only if the property cannot be found
    if name == "trial":
      if self._trial is None and self.trial_id is not None:
        self._trial = self.app.get_trial(trial_id = self.trial_id)
      return self._trial
    elif name == "scope":
      return task_scope[task_types.index(self.type)]
    else:
      return super(Task, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name not in ['scope']
    if name == 'trial':
      self.trial_id = value.trial_id
      self._trial = value
    else:
      super(Task, self).__setattr__(name, value)

  @staticmethod
  def get_phil_scope(app, task_type):
    assert task_type in task_types
    if task_type == "indexing":
      dispatcher = app.params.dispatcher
      if dispatcher == 'cxi.xtc_process': #LABELIT
        from spotfinder.applications.xfel import cxi_phil
        return cxi_phil.cxi_versioned_extract().persist.phil_scope
    else:
      dispatcher = task_dispatchers[task_types.index(task_type)]

    if dispatcher is None:
      return dispatcher, None

    return dispatcher, load_phil_scope_from_dispatcher(dispatcher)

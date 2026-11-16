from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy
from xfel.ui import load_phil_scope_from_dispatcher

task_types = ["indexing", "ensemble_refinement", "scaling", "merging", "phenix"]
task_dispatchers = [None, "cctbx.xfel.time_varying_refinement", "cctbx.xfel.merge", "cctbx.xfel.merge", None]
task_scope = ["local", "local", "local", "global", "global"]

def canonical_parameters(task_type, parameters, phil_scope = None):
  """ The form of a task's parameters text used to decide whether two tasks
      are the same. PHIL-backed types are re-emitted as a diff against the
      dispatcher's scope so formatting differences (wrapping, choice rendering,
      ordering) do not count; phenix parameters are raw text. Text that does
      not parse compares as itself. """
  parameters = (parameters or "").strip()
  if phil_scope is None or task_type in ("indexing", "phenix"):
    return parameters
  from libtbx.phil import parse
  try:
    working = phil_scope.fetch(parse(parameters))
    # Round-trip through the extracted python values so equivalent literals
    # (2 and 2.0, quoted and bare strings) render the same way.
    formatted = phil_scope.format(python_object=working.extract())
    return phil_scope.fetch_diff(formatted).as_str(print_width=100000).strip()
  except Exception:
    return parameters

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
    elif name == "parent":
      if self.parent_task_id is None:
        return None
      return self.app.get_task(self.parent_task_id)
    elif name == "scope":
      return task_scope[task_types.index(self.type)]
    else:
      return super(Task, self).__getattr__(name)

  def __setattr__(self, name, value):
    assert name not in ['scope', 'parent']
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

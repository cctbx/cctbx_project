from __future__ import absolute_import, division, print_function
from xfel.ui.db import db_proxy
from xfel.ui import load_phil_scope_from_dispatcher

task_types = ["indexing", "ensemble_refinement", "scaling", "merging"]
task_dispatchers = [None, "cctbx.xfel.stripe_experiment", "cctbx.xfel.merge", "cctbx.xfel.merge"]

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
    else:
      return super(Task, self).__getattr__(name)

  def __setattr__(self, name, value):
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

    return dispatcher, load_phil_scope_from_dispatcher(dispatcher)


def submit_ensemble_refinement_job(app, job, identifier_string):
  from xfel.command_line.striping import Script
  from xfel.ui import settings_dir
  from xfel.command_line.cxi_mpi_submit import get_submission_id
  from xfel.ui.db import get_run_path
  from libtbx import easy_run
  import os
  configs_dir = os.path.join(settings_dir, "cfgs")
  target_phil_path = os.path.join(configs_dir, identifier_string + "_params.phil")
  with open(target_phil_path, 'w') as f:
    f.write(job.task.parameters)

  path = get_run_path(app.params.output_folder, job.trial, job.rungroup, job.run, job.task)
  os.mkdir(path)

  arguments = """
  mp.queue={}
  mp.nproc={}
  striping.results_dir={}
  striping.trial={}
  striping.run_group={}
  striping.run={}
  {}
  striping.chunk_size=3000
  striping.stripe=False
  striping.dry_run=True
  striping.output_folder={}
  """.format(app.params.mp.queue,
             app.params.mp.nproc,
             app.params.output_folder,
             job.trial.trial,
             job.rungroup.id,
             job.run.run,
             target_phil_path,
             path,
             ).split()

  commands = Script(arguments).run()
  submission_ids = []
  for command in commands:
    try:
      result = easy_run.fully_buffered(command=command)
      result.raise_if_errors()
    except Exception as e:
      if not "Warning: job being submitted without an AFS token." in str(e):
        raise e
    submission_ids.append(get_submission_id(result, app.params.mp.method))
  return ",".join(submission_ids)


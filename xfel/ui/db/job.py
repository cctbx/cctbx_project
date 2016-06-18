from __future__ import division
from xfel.ui.db import db_proxy

class Job(db_proxy):
  def __init__(self, app, job_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_job" % app.params.experiment_tag, id = job_id, **kwargs)
    self.job_id = self.id


# Support classes and functions for job submission

class _job(object):
  """Used to represent a job that may not have been submitted into the cluseter or database  yet"""
  def __init__(self, trial, rungroup, run):
    self.trial = trial
    self.rungroup = rungroup
    self.run = run

  def __eq__(self, other):
    return self.trial.id == other.trial_id and \
           self.rungroup.id == other.rungroup_id and \
           self.run.id == other.run_id

def submit_all_jobs(app):
  runs = app.get_all_runs()
  submitted_jobs = app.get_all_jobs()
  trials = app.get_all_trials(only_active = True)

  needed_jobs = []
  for trial in trials:
    for rungroup in trial.rungroups:
      if not rungroup.active: continue
      rg_start = app.get_run(run_id=rungroup.startrun)
      if rungroup.endrun is None:
        # open ended run group
        rg_runs = [r for r in runs if r.run >= rg_start.run]
      else:
        # closed run group
        rg_end = app.get_run(run_id=rungroup.endrun)
        rg_runs = [r for r in runs if r.run >= rg_start.run and r.run <= rg_end.run]
      for run in rg_runs:
        needed_jobs.append(_job(trial, rungroup, run))

  for job in needed_jobs:
    if job in submitted_jobs:
      continue

    print "Submitting job: trial %d, rungroup %d, run %d"%(job.trial.id, job.rungroup.id, job.run.id)
    app.create_job(
      trial_id = job.trial.id,
      rungroup_id = job.rungroup.id,
      run_id = job.run.id,
      status = "Submitted")

    submit_job(app, job)

def submit_job(app, job):
  import os, libtbx.load_env
  from xfel.ui import settings_dir
  configs_dir = os.path.join(settings_dir, "cfgs")
  if not os.path.exists(configs_dir):
    os.makedirs(configs_dir)

  target_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_params.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))
  phil = open(target_phil_path, "w")
  phil.write(job.trial.target_phil_str)
  phil.close()

  submit_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_submit.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))

  template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "submit.phil"))
  phil = open(submit_phil_path, "w")

  app.params.mp.nproc = 12
  app.params.mp.queue = "psanaq"

  d = dict(dry_run = app.params.dry_run,
    experiment = app.params.experiment,
    experiment_tag = app.params.experiment_tag,
    run_num = job.run.run,
    trial = job.trial.trial,
    rungroup = job.rungroup.rungroup_id,
    output_dir = app.params.output_folder,
    nproc = app.params.mp.nproc,
    queue = app.params.mp.queue,
    target = target_phil_path,
    host = app.params.db.host,
    dbname = app.params.db.name,
    user = app.params.db.user,
  )
  if app.params.db.password is not None and len(app.params.db.password) == 0:
    d['password'] = None
  else:
    d['password'] = app.params.db.password

  for line in template.readlines():
    phil.write(line.format(**d))

  template.close()
  phil.close()

  from xfel.command_line.cxi_mpi_submit import Script as submit_script
  submit_script().run([submit_phil_path])

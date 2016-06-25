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

  from xfel.command_line.xtc_process import phil_scope
  from iotbx.phil import parse
  trial_params = phil_scope.fetch(parse(job.trial.target_phil_str)).extract()
  if trial_params.format.file_format == "cbf":
    trial_params.format.cbf.detz_offset = job.rungroup.detz_parameter
    trial_params.format.cbf.override_energy = job.rungroup.energy
    trial_params.format.cbf.invalid_pixel_mask = job.rungroup.untrusted_pixel_mask_path
    trial_params.format.cbf.gain_mask_value = job.rungroup.gain_mask_level
  else:
    assert False
  trial_params.dispatch.process_percent = job.trial.process_percent

  working_phil = phil_scope.format(python_object=trial_params)
  diff_phil = phil_scope.fetch_diff(source=working_phil)

  target_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_params.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))
  phil = open(target_phil_path, "w")
  phil.write(diff_phil.as_str())
  phil.close()

  config_path = None
  if job.rungroup.calib_dir is not None or job.rungroup.config_str is not None:
    config_str = "[psana]\n"
    if job.rungroup.calib_dir is not None:
      config_str += "calib-dir=%s\n"%job.rungroup.calib_dir
    if job.rungroup.config_str is not None:
      modules = []
      for line in job.rungroup.config_str.split("\n"):
        assert isinstance(line, str)
        if line.startswith('['):
          modules.append(line.lstrip('[').rstrip(']'))
      assert len(modules) > 0
      config_str += "modules = %s\n"%(" ".join(modules))
      config_str += job.rungroup.config_str

    config_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d.cfg"%
      (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))
    cfg = open(config_path, 'w')
    cfg.write(config_str)
    cfg.close()

  submit_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_submit.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))

  template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "submit.phil"))
  phil = open(submit_phil_path, "w")

  d = dict(dry_run = app.params.dry_run,
    cfg = config_path,
    calib_dir = job.rungroup.calib_dir,
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

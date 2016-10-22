from __future__ import division
from xfel.ui.db import db_proxy

class Job(db_proxy):
  def __init__(self, app, job_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_job" % app.params.experiment_tag, id = job_id, **kwargs)
    self.job_id = self.id

  def get_log_path(self):
    from xfel.ui.db import get_run_path
    import os
    run_path = str(get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run))
    return os.path.join(run_path, "stdout", "log.out")

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
      assert rungroup.active
      rg_start = app.get_run(run_number=rungroup.startrun)
      if rungroup.endrun is None:
        # open ended run group
        rg_runs = [r for r in runs if r.run >= rg_start.run]
      else:
        # closed run group
        rg_end = app.get_run(run_number=rungroup.endrun)
        rg_runs = [r for r in runs if r.run >= rg_start.run and r.run <= rg_end.run]
      for run in rg_runs:
        needed_jobs.append(_job(trial, rungroup, run))

  all_jobs = [j for j in submitted_jobs] # shallow copy
  for job in needed_jobs:
    if job in submitted_jobs:
      continue

    print "Submitting job: trial %d, rungroup %d, run %d"%(job.trial.trial, job.rungroup.id, job.run.run)

    j = app.create_job(trial_id = job.trial.id,
                       rungroup_id = job.rungroup.id,
                       run_id = job.run.id,
                       status = "SUBMITTED")
    all_jobs.append(j)
    try:
      j.submission_id = submit_job(app, job)
    except Exception, e:
      print "Couldn't submit job:", str(e)
      j.status = "SUBMIT_FAIL"

def submit_job(app, job):
  import os, libtbx.load_env
  from xfel.ui import settings_dir
  configs_dir = os.path.join(settings_dir, "cfgs")
  if not os.path.exists(configs_dir):
    os.makedirs(configs_dir)
  target_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_params.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))
  backend = ['labelit', 'dials'][['cxi.xtc_process', 'cctbx.xfel.xtc_process'].index(app.params.dispatcher)]

  phil_str = job.trial.target_phil_str
  if job.rungroup.extra_phil_str is not None:
    phil_str += "\n" + job.rungroup.extra_phil_str

  if backend == 'dials':
    from xfel.command_line.xtc_process import phil_scope
    from iotbx.phil import parse
    trial_params = phil_scope.fetch(parse(phil_str)).extract()
    image_format = trial_params.format.file_format
    assert image_format in ['cbf', 'pickle']
  else:
    image_format = 'pickle'

  if job.rungroup.calib_dir is not None or job.rungroup.config_str is not None or backend == 'labelit' or image_format == 'pickle':
    config_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d.cfg"%
      (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))
  else:
    config_path = None

  # Dictionary for formating the submit phil and, if used, the labelit cfg file
  d = dict(
    # Generally for the LABELIT backend or image pickles
    address                   = job.rungroup.detector_address,
    default_calib_dir         = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0"),
    dark_avg_path             = job.rungroup.dark_avg_path,
    dark_stddev_path          = job.rungroup.dark_stddev_path,
    untrusted_pixel_mask_path = job.rungroup.untrusted_pixel_mask_path,
    detz_parameter            = job.rungroup.detz_parameter,
    gain_map_path             = job.rungroup.gain_map_path,
    gain_mask_level           = job.rungroup.gain_mask_level,
    beamx                     = job.rungroup.beamx,
    beamy                     = job.rungroup.beamy,
    energy                    = job.rungroup.energy,
    binning                   = job.rungroup.binning,
    two_theta_low             = 12.5, # FIXME
    two_theta_high            = 22.8, # FIXME
    # Generally for job submission
    dry_run                   = app.params.dry_run,
    dispatcher                = app.params.dispatcher,
    cfg                       = config_path,
    experiment                = app.params.experiment,
    run_num                   = job.run.run,
    output_dir                = app.params.output_folder,
    use_ffb                   = app.params.use_ffb,
    # Generally for both
    trial                     = job.trial.trial,
    rungroup                  = job.rungroup.rungroup_id,
    experiment_tag            = app.params.experiment_tag,
    calib_dir                 = job.rungroup.calib_dir,
    nproc                     = app.params.mp.nproc,
    queue                     = app.params.mp.queue,
    target                    = target_phil_path,
    host                      = app.params.db.host,
    dbname                    = app.params.db.name,
    user                      = app.params.db.user,
  )
  if app.params.db.password is not None and len(app.params.db.password) == 0:
    d['password'] = None
  else:
    d['password'] = app.params.db.password

  phil = open(target_phil_path, "w")

  if backend == 'dials':
    if trial_params.format.file_format == "cbf":
      trial_params.format.cbf.detz_offset = job.rungroup.detz_parameter
      trial_params.format.cbf.override_energy = job.rungroup.energy
      trial_params.format.cbf.invalid_pixel_mask = job.rungroup.untrusted_pixel_mask_path
      trial_params.format.cbf.gain_mask_value = job.rungroup.gain_mask_level
    trial_params.dispatch.process_percent = job.trial.process_percent

    working_phil = phil_scope.format(python_object=trial_params)
    diff_phil = phil_scope.fetch_diff(source=working_phil)

    phil.write(diff_phil.as_str())
  elif backend == 'labelit':
    phil.write(phil_str)
  else:
    assert False
  phil.close()

  if config_path is not None:
    if backend == 'dials':
      d['untrusted_pixel_mask_path'] = None # Don't pass a pixel mask to mod_image_dict as it will
                                            # will be used during dials processing directly

    config_str = "[psana]\n"
    if job.rungroup.calib_dir is not None:
      config_str += "calib-dir=%s\n"%job.rungroup.calib_dir
    modules = []
    if job.rungroup.config_str is not None:
      for line in job.rungroup.config_str.split("\n"):
        if line.startswith('['):
          modules.append(line.lstrip('[').rstrip(']'))
    if backend == 'labelit':
      modules.insert(0, 'my_ana_pkg.mod_radial_average')
      modules.extend(['my_ana_pkg.mod_hitfind:index','my_ana_pkg.mod_dump:index'])
    elif image_format == 'pickle':
      modules.extend(['my_ana_pkg.mod_image_dict'])
    if app.params.dump_shots:
      modules.insert(0, 'my_ana_pkg.mod_dump:shot')

    if len(modules) > 0:
      config_str += "modules = %s\n"%(" ".join(modules))

    if job.rungroup.config_str is not None:
      config_str += job.rungroup.config_str + "\n"

    if backend == 'labelit' or image_format == 'pickle':
      d['address'] = d['address'].replace('.','-').replace(':','|') # old style address
      if backend == 'labelit':
        template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "index_all.cfg"))
      elif image_format == 'pickle':
        template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "image_dict.cfg"))
      for line in template.readlines():
        config_str += line.format(**d)
      template.close()
      d['address'] = job.rungroup.detector_address

    cfg = open(config_path, 'w')
    cfg.write(config_str)
    cfg.close()

    if backend == 'dials':
      d['untrusted_pixel_mask_path'] = job.rungroup.untrusted_pixel_mask_path

  submit_phil_path = os.path.join(configs_dir, "%s_%s_r%04d_t%03d_rg%03d_submit.phil"%
    (app.params.experiment, app.params.experiment_tag, job.run.run, job.trial.trial, job.rungroup.id))

  template = open(os.path.join(libtbx.env.find_in_repositories("xfel/ui/db/cfgs"), "submit.phil"))
  phil = open(submit_phil_path, "w")

  if backend == 'labelit':
    d['target'] = None # any target phil will be in mod_hitfind

  for line in template.readlines():
    phil.write(line.format(**d))

  d['target'] = target_phil_path

  template.close()
  phil.close()

  from xfel.command_line.cxi_mpi_submit import Script as submit_script
  return submit_script().run([submit_phil_path])

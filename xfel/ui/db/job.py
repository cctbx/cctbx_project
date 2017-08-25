from __future__ import division
from xfel.ui.db import db_proxy

known_job_statuses = ["DONE", "ERR", "PEND", "RUN", "SUSP", "PSUSP", "SSUSP", "UNKWN", "EXIT", "DONE", "ZOMBI", "DELETED", "SUBMIT_FAIL", "SUBMITTED"]
finished_job_statuses = ["DONE", "EXIT", "DELETED", "UNKWN", "ERR", "SUBMIT_FAIL"]

class Job(db_proxy):
  def __init__(self, app, job_id = None, **kwargs):
    db_proxy.__init__(self, app, "%s_job" % app.params.experiment_tag, id = job_id, **kwargs)
    self.job_id = self.id

  def get_log_path(self):
    from xfel.ui.db import get_run_path
    import os
    run_path = str(get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run))
    return os.path.join(run_path, "stdout", "log.out")

  def delete(self, output_only=False):
    from xfel.ui.db import get_run_path
    import os, shutil

    if self.status not in finished_job_statuses:
      print "Job is not finished (status = %s)"%self.status
      return

    if self.status == "DELETED":
      return

    job_folder = get_run_path(self.app.params.output_folder, self.trial, self.rungroup, self.run)
    if os.path.exists(job_folder):
      print "Deleting job folder for job", self.id
      shutil.rmtree(job_folder)
    else:
      print "Cannot find job folder (%s)"%job_folder

    # Have to be careful to delete from the tables in the right order
    tag = self.app.params.experiment_tag

    def delete_and_commit(query):
      cursor = self.app.execute_query(query, commit=True)
      print "(%d)"%cursor.rowcount

    print "Deleting cell_bin entries",
    query = """DELETE cell_bin FROM `%s_cell_bin` cell_bin
               JOIN `%s_crystal` crystal ON crystal.id = cell_bin.crystal_id
               JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    ids = {}
    for item in "crystal", "beam", "detector":
      print "Listing %s ids"%item,
      query = """SELECT %s.id FROM `%s_%s` %s
                 JOIN `%s_experiment` expr ON expr.%s_id = %s.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
                 item, tag, item, item, tag, item, item, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print "(%d)"%len(item_ids)
      ids[item] = ",".join(item_ids)

    if len(self.trial.isoforms) == 0:
      print "Listing bin entries",
      query = """SELECT bin.id FROM `%s_bin` bin
                 JOIN `%s_cell` cell ON bin.cell_id = cell.id
                 JOIN `%s_crystal` crystal ON crystal.cell_id = cell.id
                 JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d
                 AND cell.trial_id is NULL""" % (
                 tag, tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print "(%d)"%len(item_ids)
      bin_ids = ",".join(item_ids)

      print "Listing cell entries",
      query = """SELECT cell.id FROM `%s_cell` cell
                 JOIN `%s_crystal` crystal ON crystal.cell_id = cell.id
                 JOIN `%s_experiment` expr ON expr.crystal_id = crystal.id
                 JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
                 JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
                 JOIN `%s_event` evt ON evt.id = ie_e.event_id
                 WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d
                 AND cell.trial_id IS NULL""" % (
                 tag, tag, tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
      cursor = self.app.execute_query(query)
      item_ids = ["%d"%i[0] for i in cursor.fetchall()]
      print "(%d)"%len(item_ids)
      cell_ids = ",".join(item_ids)

    print "Deleting experiment entries",
    query = """DELETE expr FROM `%s_experiment` expr
               JOIN `%s_imageset` imgset ON imgset.id = expr.imageset_id
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    for item in "crystal", "beam", "detector":
      if len(ids[item]) > 0:
        print "Deleting %s entries"%item,
        query = """DELETE %s FROM `%s_%s` %s
                   WHERE %s.id IN (%s)""" % (
                   item, tag, item, item, item, ids[item])
        delete_and_commit(query)

    if len(self.trial.isoforms) == 0 and len(bin_ids) > 0:
      print "Deleting bin entries",
      query = """DELETE bin FROM `%s_bin` bin
                 WHERE bin.id IN (%s)""" % (
                 tag, bin_ids)
      delete_and_commit(query)

    if len(self.trial.isoforms) == 0 and len(cell_ids) > 0:
      print "Deleting cell entries",
      query = """DELETE cell FROM `%s_cell` cell
                 WHERE cell.id IN (%s)""" % (
                 tag, cell_ids)
      delete_and_commit(query)

    print "Listing imageset entries",
    query = """SELECT imgset.id FROM `%s_imageset` imgset
               JOIN `%s_imageset_event` ie_e ON ie_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = ie_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    cursor = self.app.execute_query(query)
    item_ids = ["%d"%i[0] for i in cursor.fetchall()]
    print "(%d)"%len(item_ids)
    imageset_ids = ",".join(item_ids)

    print "Deleting imageset_event entries",
    query = """DELETE is_e FROM `%s_imageset_event` is_e
               JOIN `%s_event` evt ON evt.id = is_e.event_id
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    if len(imageset_ids) > 0:
      print "Deleting imageset entries",
      query = """DELETE imgset FROM `%s_imageset` imgset
                 WHERE imgset.id IN (%s)""" % (
                 tag, imageset_ids)
      delete_and_commit(query)

    print "Deleting event entries",
    query = """DELETE evt FROM `%s_event` evt
               WHERE evt.run_id = %d AND evt.trial_id = %d AND evt.rungroup_id = %d""" % (
               tag, self.run.id, self.trial.id, self.rungroup.id)
    delete_and_commit(query)

    self.status = "DELETED"

  def remove_from_db(self):
    assert self.status == "DELETED"

    print "Removing job %d from the db"%self.id,
    tag = self.app.params.experiment_tag
    query = """DELETE job FROM `%s_job` job
               WHERE job.id = %d""" % (
               tag, self.id)
    cursor = self.app.execute_query(query, commit=True)
    print "(%d)"%cursor.rowcount

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
    from xfel.command_line.xtc_process import phil_scope as orig_phil_scope
    from iotbx.phil import parse
    if job.rungroup.two_theta_low is not None or job.rungroup.two_theta_high is not None:
      override_str = """
      radial_average {
        enable = True
        show_plots = False
        verbose = False
        output_bins = False
      }
      """
      phil_scope = orig_phil_scope.fetch(parse(override_str))
    else:
      phil_scope = orig_phil_scope

    trial_params = phil_scope.fetch(parse(phil_str)).extract()
    image_format = trial_params.format.file_format = job.rungroup.format
    assert image_format in ['cbf', 'pickle']
    if image_format == 'cbf':
      if "rayonix" in job.rungroup.detector_address.lower():
        mode = trial_params.format.cbf.mode = "rayonix"
      elif "cspad" in job.rungroup.detector_address.lower():
        mode = trial_params.format.cbf.mode = "cspad"
      else:
        assert False, "Couldn't figure out what kind of detector is specified by address %s"%job.rungroup.detector_address
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
    two_theta_low             = job.rungroup.two_theta_low,
    two_theta_high            = job.rungroup.two_theta_high,
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
      trial_params.input.address = job.rungroup.detector_address
      trial_params.format.cbf.detz_offset = job.rungroup.detz_parameter
      trial_params.format.cbf.override_energy = job.rungroup.energy
      trial_params.format.cbf.invalid_pixel_mask = job.rungroup.untrusted_pixel_mask_path
      if mode == 'cspad':
        trial_params.format.cbf.cspad.gain_mask_value = job.rungroup.gain_mask_level
      elif mode == 'rayonix':
        trial_params.format.cbf.rayonix.bin_size = job.rungroup.binning
        trial_params.format.cbf.rayonix.override_beam_x = job.rungroup.beamx
        trial_params.format.cbf.rayonix.override_beam_y = job.rungroup.beamy
    trial_params.dispatch.process_percent = job.trial.process_percent

    if job.rungroup.two_theta_low is not None or job.rungroup.two_theta_high is not None:
      trial_params.radial_average.two_theta_low = job.rungroup.two_theta_low
      trial_params.radial_average.two_theta_high = job.rungroup.two_theta_high

    working_phil = phil_scope.format(python_object=trial_params)
    diff_phil = orig_phil_scope.fetch_diff(source=working_phil)

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
      modules.insert(0, 'my_ana_pkg.mod_radial_average')
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

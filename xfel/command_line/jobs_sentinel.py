# LIBTBX_SET_DISPATCHER_NAME xpp.jobs_sentinel

"""
Program to monitor for new runs and submit them for processing. Default is to check for
and submit new runs at 0.1 Hz.

Example usage:
xpp.jobs_sentinel db.name=xppi6115 db.user=xppi6115 experiment=xppi6115 experiment_tag=debug
"""
from __future__ import absolute_import, division, print_function

import iotbx.phil
import libtbx.load_env
from libtbx.utils import Usage, Sorry
import sys, os, time
from six.moves import zip

master_phil = """
  experiment = None
    .type = str
  experiment_tag = None
    .type = str
  config_dir = "."
    .type = str
  output_dir = "."
    .type = str
  dry_run = False
    .type = bool
  mp {
    nproc = 1
      .type = int
      .help = Number of processes
    queue = "psanacsq"
      .type = str
      .help = Queue to submit MPI job to
  }
  db {
    host = psdb.slac.stanford.edu
      .type = str
    name = None
      .type = str
    user = None
      .type = str
    password = None
      .type = str
  }
"""

def parse_entry(entry):
  result = []
  for item in entry:
    if item == "NULL":
      result.append(None)
    else:
      result.append(item)
  return result

def write_hitfind_cfg(params, dbobj, trial_id, trial, rungroup_id):
  assert os.path.exists(params.config_dir)

  config_path = os.path.join(params.config_dir, "%s_%s_t%03d_rg%03d.cfg"%(params.experiment, params.experiment_tag, trial, rungroup_id))

  cmd = "SELECT target_phil_path FROM %s_trials where %s_trials.trial_id = %d"%(params.experiment_tag, params.experiment_tag, trial_id)
  cursor = dbobj.cursor()
  cursor.execute(cmd)
  assert cursor.rowcount == 1
  target_phil_path = parse_entry(cursor.fetchall()[0])[0]

  cmd = "SELECT detz_parameter, beamx, beamy, untrusted_pixel_mask_path, dark_avg_path, dark_stddev_path, gain_map_path, binning FROM %s_rungroups where %s_rungroups.rungroup_id = %d"%(params.experiment_tag, params.experiment_tag, rungroup_id)
  cursor = dbobj.cursor()
  cursor.execute(cmd)
  assert cursor.rowcount == 1
  detz_parameter, beamx, beamy, untrusted_pixel_mask_path, dark_avg_path, dark_stddev_path, gain_map_path, binning = parse_entry(cursor.fetchall()[0])

  template = open(os.path.join(libtbx.env.find_in_repositories("xfel/xpp/cfgs"),"index_all.cfg"))
  cfg = open(config_path, 'w')

  d = dict(
    default_calib_dir         = libtbx.env.find_in_repositories("xfel/metrology/CSPad/run4/CxiDs1.0_Cspad.0"),
    trial_id                  = trial_id,
    rungroup_id               = rungroup_id,
    dark_avg_path             = dark_avg_path,
    dark_stddev_path          = dark_stddev_path,
    untrusted_pixel_mask_path = untrusted_pixel_mask_path,
    detz_parameter            = detz_parameter,
    gain_map_path             = gain_map_path,
    beamx                     = beamx,
    beamy                     = beamy,
    binning                   = binning,
    db_name                   = params.db.name,
    db_experiment_tag         = params.experiment_tag,
    db_user                   = params.db.user,
    db_password               = params.db.password,
    target_phil_path          = target_phil_path)

  for line in template.readlines():
    cfg.write(line.format(**d))

  template.close()
  cfg.close()
  return config_path

def submit_job(params, dbobj, trial_id, trial, rungroup_id, run, config_path):
  assert os.path.exists(config_path)
  submit_phil_path = os.path.join(params.config_dir, "%s_%s_r%04d_t%03d_rg%03d_submit.phil"%(params.experiment, params.experiment_tag, run, trial, rungroup_id))


  template = open(os.path.join(libtbx.env.find_in_repositories("xfel/xpp/cfgs"), "submit.phil"))
  phil = open(submit_phil_path, "w")

  d = dict(dry_run = params.dry_run,
    cfg = config_path,
    experiment = params.experiment,
    run_num = run,
    trial = trial,
    rungroup = rungroup_id,
    output_dir = params.output_dir,
    nproc = params.mp.nproc,
    queue = params.mp.queue
  )

  for line in template.readlines():
    phil.write(line.format(**d))

  template.close()
  phil.close()

  from xfel.command_line.submit_job import Script as submit_script
  submit_script().run([submit_phil_path])

def run(args):
  try:
    from cxi_xdr_xes.cftbx.cspad_ana import db as db
  except ImportError:
    raise Sorry("Trial logging not supported for this installation. Contact the developers for access.")

  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil)
  params = phil.work.extract()

  if params.db.host is None:
    raise Usage("Please provide a host name")
  if params.db.name is None:
    raise Usage("Please provide a database name")
  if params.db.user is None:
    raise Usage("Please provide a user name")
  if params.db.password is None:
    import getpass
    password = getpass.getpass()
  else:
    password = params.db.password

  while True:
    print("Checking for new jobs to submit")
    # need to get a new connection each iteration to prevent caching
    try:
      dbobj = db.dbconnect(host=params.db.host, db=params.db.name, username=params.db.user, password=password)
    except Exception as e:
      raise Sorry(e)

    # Get the set of known runs in the experiment database
    cmd = "SELECT run from %s_runs"%params.experiment_tag
    cursor = dbobj.cursor()
    cursor.execute(cmd)
    known_runs = [int(parse_entry(entry)[0]) for entry in cursor.fetchall()]

    # Get the list of active trials
    cmd = "SELECT trial_id, trial from %s_trials WHERE active = True"%params.experiment_tag
    cursor = dbobj.cursor()
    cursor.execute(cmd)
    entries = [parse_entry(entry) for entry in cursor.fetchall()]
    active_trial_ids = [int(entry[0]) for entry in entries]
    active_trials = [int(entry[1]) for entry in entries]

    for trial_id, trial in zip(active_trial_ids, active_trials):
      # Get the active rungroups for this trial
      cmd = "SELECT rungroups_id from %s_trial_rungroups WHERE trials_id = %d AND active = True"%(params.experiment_tag, trial_id)
      cursor = dbobj.cursor()
      cursor.execute(cmd)
      active_rungroup_ids = [int(parse_entry(entry)[0]) for entry in cursor.fetchall()]

      for rungroup_id in active_rungroup_ids:
        # Get the list of runs for this rungroup
        cmd = "SELECT startrun, endrun from %s_rungroups WHERE %s_rungroups.rungroup_id = %d"%(params.experiment_tag, params.experiment_tag, rungroup_id)
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        assert cursor.rowcount == 1
        startrun, endrun = parse_entry(cursor.fetchall()[0])

        cmd = "SELECT run_id, run from %s_runs WHERE run >= %d"%(params.experiment_tag, startrun)
        if endrun is not None:
          cmd += " AND run <= %d"%endrun
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        entries = [parse_entry(entry) for entry in cursor.fetchall()]
        run_ids = [int(entry[0]) for entry in entries]
        runs = [int(entry[1]) for entry in entries]

        # Find the submitted runs for this trial/rungroup combination
        cmd = "SELECT runs_id from %s_jobs WHERE trials_id = %d and rungroups_id = %d"%(params.experiment_tag, trial_id, rungroup_id)
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        submitted_run_ids = [int(parse_entry(entry)[0]) for entry in cursor.fetchall()]

        # Submit new runs
        for run_id, run in zip(run_ids, runs):
          if run_id in submitted_run_ids:
            continue

          print("Submitting run %d into trial %d using rungroup %d"%(run_id, trial, rungroup_id))

          config_path = write_hitfind_cfg(params, dbobj, trial_id, trial, rungroup_id)
          if submit_job(params, dbobj, trial_id, trial, rungroup_id, run, config_path):
            pass

          cmd = "INSERT INTO %s_jobs (runs_id, trials_id, rungroups_id, status) VALUES (%d, %d, %d, '%s')"%(params.experiment_tag, run_id, trial_id, rungroup_id, "submitted")
          cursor = dbobj.cursor()
          cursor.execute(cmd)
          dbobj.commit()

    time.sleep(10)

if __name__ == "__main__":
  run(sys.argv[1:])

from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME xpp.jobs_sentinel

"""
Program to monitor for new runs and submit them for processing. Default is to check for
and submit new runs at 0.1 Hz.

Example usage:
xpp.jobs_sentinel db.name=xppi6115 db.user=xppi6115 experiment=xppi6115 experiment_tag=debug
"""

import iotbx.phil
from libtbx.utils import Usage, Sorry
import sys, time

master_phil = """
  experiment = None
    .type = str
  experiment_tag = None
    .type = str
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

def write_hitfind_cfg(trial, rungroup):
  pass

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
    print "Checking for new jobs to submit"
    # need to get a new connection each iteration to prevent caching
    try:
      dbobj = db.dbconnect(host=params.db.host, db=params.db.name, username=params.db.user, password=password)
    except Exception, e:
      raise Sorry(e)

    # Get the set of known runs in the experiment database
    cmd = "SELECT run from %s_runs"%params.experiment_tag
    cursor = dbobj.cursor()
    cursor.execute(cmd)
    known_runs = [int(entry[0]) for entry in cursor.fetchall()]

    # Get the list of active trials
    cmd = "SELECT trial_id from %s_trials WHERE active = True"%params.experiment_tag
    cursor = dbobj.cursor()
    cursor.execute(cmd)
    active_trial_ids = [int(entry[0]) for entry in cursor.fetchall()]

    for trial_id in active_trial_ids:
      # Get the active rungroups for this trial
      cmd = "SELECT rungroups_id from %s_trial_rungroups WHERE trials_id = %d AND active = True"%(params.experiment_tag, trial_id)
      cursor = dbobj.cursor()
      cursor.execute(cmd)
      active_rungroup_ids = [int(entry[0]) for entry in cursor.fetchall()]

      for rungroup_id in active_rungroup_ids:
        # Get the list of runs for this rungroup
        cmd = "SELECT startrun, endrun from %s_rungroups WHERE %s_rungroups.rungroup_id = %d"%(params.experiment_tag, params.experiment_tag, rungroup_id)
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        assert cursor.rowcount == 1
        startrun, endrun = cursor.fetchall()[0]

        cmd = "SELECT run_id, run from %s_runs WHERE run >= %d"%(params.experiment_tag, startrun)
        if endrun is not None:
          cmd += " AND run <= %d"%endrun
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        entries = cursor.fetchall()
        run_ids = [int(entry[0]) for entry in entries]
        runs = [int(entry[1]) for entry in entries]

        # Find the submitted runs for this trial/rungroup combination
        cmd = "SELECT runs_id from %s_jobs WHERE trials_id = %d and rungroups_id = %d"%(params.experiment_tag, trial_id, rungroup_id)
        cursor = dbobj.cursor()
        cursor.execute(cmd)

        submitted_run_ids = [int(entry[0]) for entry in cursor.fetchall()]

        # Submit new runs
        for run_id, run in zip(run_ids, runs):
          if run_id in submitted_run_ids:
            continue

          print "Submitting run %d into trial %d using rungroup %d"%(run_id, trial_id, rungroup_id)

          cmd = "INSERT INTO %s_jobs (runs_id, trials_id, rungroups_id, status) VALUES (%d, %d, %d, '%s')"%(params.experiment_tag, run_id, trial_id, rungroup_id, "submitted")
          cursor = dbobj.cursor()
          cursor.execute(cmd)
          dbobj.commit()

    time.sleep(10)

if __name__ == "__main__":
  run(sys.argv[1:])

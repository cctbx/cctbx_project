# LIBTBX_SET_DISPATCHER_NAME xpp.runs_sentinel

"""
Program to monitor runs using SLACs dataexport webservice and to log them in the
experiment managment database when they arrive. Default is to check for and log
new runs at 0.1 Hz.

Example usage:
xpp.runs_sentinel db.name=xppi6115 db.user=xppi6115 experiment=xppi6115 experiment_tag=debug web.user=<username> web.password=<password>
"""
from __future__ import absolute_import, division, print_function

import iotbx.phil
from libtbx.utils import Usage, Sorry
import sys, time

master_phil = """
  experiment = None
    .type = str
  experiment_tag = None
    .type = str
  run_tags = None
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
  web {
    user = None
      .type = str
    password = None
      .type = str
  }
"""

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

  try:
    dbobj = db.dbconnect(host=params.db.host, db=params.db.name, username=params.db.user, password=password)
  except Exception as e:
    raise Sorry(e)

  from xfel.xpp.simulate import file_table
  query = "https://pswww.slac.stanford.edu/ws-auth/dataexport/placed?exp_name=%s"%(params.experiment)

  # set up extra run tags, if provided
  if params.run_tags is not None:
    extra1 = ", tags"
    extra2 = ",'%s'"%params.run_tags
  else:
    extra1 = ""
    extra2 = ""

  while True:
    # Get the set of known runs in the experiment database
    cmd = "SELECT run from %s_runs"%params.experiment_tag
    cursor = dbobj.cursor()
    cursor.execute(cmd)

    # Get the set of runs from SLAC's database
    FT = file_table(params,query)

    # Find the delta
    known_runs = [int(entry[0]) for entry in cursor.fetchall()]
    unknown_runs = [run for run in FT.rundict if run not in known_runs]

    print("%d new runs"%len(unknown_runs))

    # Enter any new runs into the experiment database
    if len(unknown_runs) > 0:
      cmd = "INSERT INTO %s_runs (run%s) VALUES "%(params.experiment_tag, extra1)
      comma = ""
      for run in unknown_runs:
        cmd += comma + "(%d%s)"%(run, extra2)
        comma = ", "

      cursor = dbobj.cursor()
      cursor.execute(cmd)
      dbobj.commit()

    time.sleep(10)

if __name__ == "__main__":
  run(sys.argv[1:])

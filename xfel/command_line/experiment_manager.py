from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME xpp.experiment_manager

import iotbx.phil
from libtbx.utils import Usage, Sorry
import sys, os
import libtbx.load_env

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

def get_bool_from_user(prompt, default=True):
  yes = ['yes','y', 'ye']
  no = ['no','n']
  if default:
    yes.append('')
    prompt += " y/n [y] "
  else:
    no.append('')
    prompt += " y/n [n] "

  while True:
    response = raw_input(prompt).lower()
    if response in yes:
      return True
    if response in no:
      return False

class initialize(object):
  expected_tables = ["runs", "jobs", "rungroups", "trials", "trial_rungroups", "isoforms", "frames", "observations", "hkls"]

  def __init__(self, params, dbobj):
    self.dbobj = dbobj
    self.params = params

  def __call__(self):
    if self.params.experiment is None:
      self.params.experiment = raw_input("Administrate which experiment? ")

    if self.params.experiment_tag is None:
      if get_bool_from_user("Use experiment name as experiment tag?"):
        self.params.experiment_tag = self.params.experiment
      else:
        self.params.experiment_tag = raw_input("Input an experiment tag: ")

    print "Administering experiment", self.params.experiment, "using tag", self.params.experiment_tag
    if get_bool_from_user("Drop exisiting tables for %s?"%self.params.experiment_tag, default=False):
      self.drop_tables()

    if not self.verify_tables():
      self.create_tables()
      if not self.verify_tables():
        raise Sorry("Couldn't create experiment tables")

  def verify_tables(self):
    print "Checking tables...",

    bools = []
    for table in self.expected_tables:
      cursor = self.dbobj.cursor()
      cmd = "SHOW TABLES LIKE '%s'"%(self.params.experiment_tag + "_" + table)
      cursor.execute(cmd)
      bools.append(cursor.rowcount > 0)

    if bools.count(True) == len(self.expected_tables):
      print "good to go"
      return True
    elif bools.count(False) == len(self.expected_tables):
      print "experiment tag not found"
      return False
    else:
      print "some tables are missing"
      return False

  def drop_tables(self):
    print "Dropping tables..."
    cmd = "SET FOREIGN_KEY_CHECKS=0; "
    for table in self.expected_tables:
      cmd += "DROP TABLE IF EXISTS %s; "%(self.params.experiment_tag + "_" + table);
    cmd += "SET FOREIGN_KEY_CHECKS=1;"
    cursor = self.dbobj.cursor()
    cursor.execute(cmd);

  def create_tables(self):
    print "Creating tables..."
    sql_path = os.path.join(libtbx.env.find_in_repositories("xfel/xpp"), "experiment_schema.sql")
    assert os.path.exists(sql_path)

    reading_create = False
    cmd = []
    f = open(sql_path)
    for line in f:
      line = line.strip()
      if not reading_create and "CREATE TABLE" in line:
        line = line.replace("`mydb`.`","`%s`.`%s_")%(self.params.db.name, self.params.experiment_tag)
        reading_create = True

      if "REFERENCES" in line:
        line = line.replace("`mydb`.`","`%s`.`%s_")%(self.params.db.name, self.params.experiment_tag)

      if reading_create:
        cmd.append(line)
        if ";" in line:
          cmd = " ".join(cmd)
          cursor = self.dbobj.cursor()
          try:
            cursor.execute(cmd)
          except Exception, e:
            print "Failed to create table. SQL command:"
            print cmd
            print e

          cmd = []
          reading_create = False

    f.close()

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
  except Exception, e:
    raise Sorry(e)

  initialize(params, dbobj)()

  print "Done"

if __name__ == "__main__":
  run(sys.argv[1:])

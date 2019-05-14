from __future__ import absolute_import, division, print_function
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

def get_optional_input(prompt):
  response = raw_input(prompt)
  if response == "":
    return "NULL"
  else:
    return response

class initialize(object):
  expected_tables = ["runs", "jobs", "rungroups", "trials", "trial_rungroups", "isoforms", "frames", "observations", "hkls"]

  def __init__(self, params, dbobj, interactive = True, drop_tables = None):
    self.dbobj = dbobj
    self.params = params
    self.interactive = interactive
    self.do_drop_tables = drop_tables

  def __call__(self):
    if self.params.experiment is None:
      self.params.experiment = raw_input("Administrate which experiment? ")

    if self.params.experiment_tag is None:
      if get_bool_from_user("Use experiment name as experiment tag?"):
        self.params.experiment_tag = self.params.experiment
      else:
        self.params.experiment_tag = raw_input("Input an experiment tag: ")

    assert self.params.experiment is not None and self.params.experiment_tag is not None and len(self.params.experiment_tag) > 0

    print("Administering experiment", self.params.experiment, "using tag", self.params.experiment_tag)
    if self.interactive and self.do_drop_tables is None:
      if get_bool_from_user("Drop existing tables for %s?"%self.params.experiment_tag, default=False):
        self.drop_tables()
    elif self.do_drop_tables == True:
      self.drop_tables()

    if not self.verify_tables():
      self.create_tables()
      if not self.verify_tables():
        raise Sorry("Couldn't create experiment tables")

  def verify_tables(self):
    print("Checking tables...", end=' ')

    bools = []
    for table in self.expected_tables:
      cursor = self.dbobj.cursor()
      cmd = "SHOW TABLES LIKE '%s'"%(self.params.experiment_tag + "_" + table)
      cursor.execute(cmd)
      bools.append(cursor.rowcount > 0)

    if bools.count(True) == len(self.expected_tables):
      print("good to go")
      return True
    elif bools.count(False) == len(self.expected_tables):
      print("experiment tag not found")
      return False
    else:
      print("some tables are missing")
      return False

  def drop_tables(self):
    if self.interactive and raw_input("Are you sure? Type drop: ").lower() != "drop":
      return

    print("Dropping tables...")
    for table in self.expected_tables:
      cmd = "SHOW TABLES LIKE '%s'"%(self.params.experiment_tag + "_" + table)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      if cursor.rowcount > 0:
        cmd = "SET FOREIGN_KEY_CHECKS=0; "
        cmd += "DROP TABLE IF EXISTS %s; "%(self.params.experiment_tag + "_" + table);
        cmd += "SET FOREIGN_KEY_CHECKS=1;"
        cursor = self.dbobj.cursor()
        cursor.execute(cmd);

  def create_tables(self, sql_path = None):
    print("Creating tables...")
    if sql_path is None:
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

      if "CONSTRAINT" in line:
        line = line.replace("`fk_","`fk_%s_")%self.params.experiment_tag

      if reading_create:
        cmd.append(line)
        if ";" in line:
          cmd = " ".join(cmd)
          cursor = self.dbobj.cursor()
          try:
            cursor.execute(cmd)
          except Exception as e:
            print("Failed to create table. SQL command:")
            print(cmd)
            print(e)

          cmd = []
          reading_create = False

    f.close()

class option_chooser(object):
  menu_string = ""
  options = None

  def __init__(self, params, dbobj):
    self.params = params
    self.dbobj = dbobj

  def __call__(self):
    while True:
      print("############################")
      print(self.menu_string)
      for o in self.options:
        print(o, ":", self.options[o].selection_string)
      print("h : get help on a selection")
      print("q : exit this menu")

      i = raw_input("Selection: ").lower()
      print("############################")
      if i == 'h':
        i = raw_input("Help on which item? ")
        if i in self.options:
          print(self.options[i].help_string)
        else:
          print("Option not found:", i)
      elif i == 'q':
        return
      elif i not in self.options:
        print("Option not found:", i)
      else:
        self.options[i](self.params, self.dbobj)()

class db_action(object):
  def __init__(self, params, dbobj):
    self.params = params
    self.dbobj = dbobj

class isoforms_menu(option_chooser):
  selection_string = "Manage isoforms"
  help_string = "Use this menu to manage known crystal isoforms refined during indexing and integration"

  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)
    self.menu_string = "Managing isoforms"
    self.options = {}

class runs_menu(option_chooser):
  selection_string = "Manage runs"
  help_string = "Use this menu to manage runs for this experiment, including adding and removing tags"

  class list_runs(db_action):
    selection_string = "List all runs"
    help_string = "List all runs being logged so far, including any tags on them"

    def __call__(self):
      cmd = "SELECT * from %s_runs"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      if cursor.rowcount == 0:
        print("No runs logged yet")
      else:
        for entry in cursor.fetchall():
          id, run, tags = entry
          print("Run %s. Tags: %s"%(run, tags))

  class update_tags(db_action):
    selection_string = "Update run tags"
    help_string = "Add or remove tags from runs"

    def get_runs_from_user(self, prompt):
      if get_bool_from_user("%s All runs? "%prompt):
        return ""
      else:
        start = raw_input("%s Starting with run: "%prompt)
        end = raw_input("%s Ending with run: "%prompt)
        return "WHERE %s_runs.run >= %s and %s_runs.run <= %s"%(
          self.params.experiment_tag, start, self.params.experiment_tag, end)

    def __call__(self):
      choice = raw_input("Add or remove tag? [a/r]: ").lower()
      if choice == 'a':
        prompt = "Add tag."
      elif choice == 'r':
        prompt = "Remove tag."
      else:
        print("Choice not recognized")
        return
      tag = raw_input("Enter tag: ")
      cursor = self.dbobj.cursor()

      cmd = "SELECT * from %s_runs %s"%(self.params.experiment_tag, self.get_runs_from_user(prompt))
      cursor.execute(cmd)
      for entry in cursor.fetchall():
        run_id, run, tags = entry
        if choice == 'a':
          if tags is None or tags == "":
            tags = ""
            comma = ""
          else:
            comma = ","
          if tag not in tags.split(','):
            tags += comma + tag
        else:
          if tags is None or tag not in tags:
            continue
          tags = tags.split(',')
          while tag in tags:
            tags.remove(tag)
          tags = ','.join(tags)

        cmd = "UPDATE %s_runs SET tags='%s' WHERE run_id=%s"%(self.params.experiment_tag, tags, run_id)
        cursor.execute(cmd)

      self.dbobj.commit()
      print("Run tags updated")

  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)
    self.menu_string = "Managing runs"
    self.options = {
      'l': runs_menu.list_runs,
      'u': runs_menu.update_tags
    }

class trials_menu(option_chooser):
  selection_string = "Manage trials"
  help_string = "Use this menu to manage trials for this experiment, including add new trials or inactivating old ones"

  class list_trials(db_action):
    selection_string = "List all trials"
    help_string = "List all trials being used with this experiment tag"

    def __call__(self):
      cmd = "SELECT * from %s_trials"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      if cursor.rowcount == 0:
        print("No trials set up yet")
      else:
        for entry in cursor.fetchall():
          id, trial, active, target, comment = entry
          active = bool(active)
          print("Trial %s. Active: %s, target: %s, comment: %s"%(trial, active, target, comment))

  class add_trial(db_action):
    selection_string = "Add trial"
    help_string = "Add a trial for use with this experiment tag"

    def __call__(self):
      trial = raw_input("Trial number: ")
      active = get_bool_from_user("Make trial active? ")
      target = raw_input("Path to target phil file: ")
      comment = raw_input("Add a comment: ")

      cmd = "INSERT INTO %s_trials (trial,active,target_phil_path,comment) VALUES (%s,%s,'%s','%s')"%(self.params.experiment_tag, trial, active, target, comment)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Trial added")

  class update_trial(db_action):
    selection_string = "Update trial"
    help_string = "Activate/inactivate a trial or change its comment"

    def __call__(self):
      trial = raw_input("Trial number: ")
      active = get_bool_from_user("Make trial active? ")
      comment = raw_input("New comment: ")

      cmd = "UPDATE %s_trials SET active=%s, comment='%s' WHERE trial=%s"%(self.params.experiment_tag, active, comment, trial)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Trial updated")

  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)
    self.menu_string = "Managing trials"
    self.options = {
      'l': trials_menu.list_trials,
      'a': trials_menu.add_trial,
      'u': trials_menu.update_trial
    }

class link_trials_to_rungroup_menu(option_chooser):
  selection_string = "Link trials to run groups"
  help_string = "Use this menu to link trials to specific run groups for processing"

  class list_links(db_action):
    selection_string = "List all trial/rungroup links"
    help_string = "List all links between trials and rungroups"

    def __call__(self):
      cmd = "SELECT %s_trial_rungroups.trial_rungroup_id, %s_trials.trial_id, %s_rungroups.rungroup_id, %s_trials.trial, %s_rungroups.startrun, %s_rungroups.endrun, %s_trial_rungroups.active FROM %s_trials JOIN %s_trial_rungroups ON %s_trials.trial_id = %s_trial_rungroups.trials_id JOIN %s_rungroups ON %s_rungroups.rungroup_id = %s_trial_rungroups.rungroups_id"%tuple([self.params.experiment_tag]*14)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      if cursor.rowcount == 0:
        print("No links set up yet")
      else:
        print("Link Trial  RG Start run End run   Active")
        for entry in cursor.fetchall():
          link_id, trial_id, rungroup_id, trial, startrun, endrun, active = entry
          active = bool(active)
          if endrun is None:
            print("% 4d % 5d % 3d % 9d       +   %s"%(link_id, trial, rungroup_id, startrun, active))
          else:
            print("% 4d % 5d % 3d % 9d % 7d   %s"%(link_id, trial, rungroup_id, startrun, endrun, active))

  class add_link(db_action):
    selection_string = "Add trial/rungroup link"
    help_string = "Link a trial to a rungroup for processing"

    def __call__(self):
      trial = raw_input("Trial: ")
      rungroup_id = raw_input("Run group id: ")
      active = get_bool_from_user("Make link active? ")

      cmd = "SELECT trial_id from %s_trials where %s_trials.trial = %s"%(self.params.experiment_tag, self.params.experiment_tag, trial)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      assert cursor.rowcount <= 1
      if cursor.rowcount == 0:
        print("Trial %s not found."%trial)
        return
      trial_id = cursor.fetchall()[0][0]

      cmd = "SELECT rungroup_id from %s_rungroups where %s_rungroups.rungroup_id = %s"%(self.params.experiment_tag, self.params.experiment_tag, rungroup_id)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      assert cursor.rowcount <= 1
      if cursor.rowcount == 0:
        print("Rungroup %s not found."%rungroup_id)
        return

      cmd = "INSERT INTO %s_trial_rungroups (trials_id, rungroups_id, active) VALUES (%s,%s,%s)"%(self.params.experiment_tag, trial_id, rungroup_id, active)

      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Link added")

  class update_link(db_action):
    selection_string = "Update trial/rungroup links"
    help_string = "Activate or inactivate trial/rungroup links"

    def __call__(self):
      link_id = raw_input("Trial/rungroup link to change: ")
      active = get_bool_from_user("Make link active? ")

      cmd = "UPDATE %s_trial_rungroups SET active=%s WHERE trial_rungroup_id=%s"%(self.params.experiment_tag, active, link_id)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Link updated")

  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)
    self.menu_string = "Managing trial/rungroup links"
    self.options = {
      'l': link_trials_to_rungroup_menu.list_links,
      'a': link_trials_to_rungroup_menu.add_link,
      'u': link_trials_to_rungroup_menu.update_link
    }

class rungroups_menu(option_chooser):
  selection_string = "Manage run groups"
  help_string = "Use this menu to manage settings on groups of runs"

  class list_rungroups(db_action):
    selection_string = "List all run groups"
    help_string = "List all run groups being used with this experiment tag"

    def __call__(self):
      cmd = "SELECT rungroup_id,startrun,endrun,comment from %s_rungroups"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      if cursor.rowcount == 0:
        print("No run groups set up yet")
      else:
        print("RG  Start run End run   Comment")
        for entry in cursor.fetchall():
          id, startrun, endrun, comment = entry
          if endrun is None:
            print("% 3d % 9d       +   %s"%(id, startrun, comment))
          else:
            print("% 3d % 9d % 7d   %s"%(id, startrun, endrun, comment))

  class add_rungroup(db_action):
    selection_string = "Add run group"
    help_string = "Add a run group defining a set of runs with the same parameters"

    def __call__(self):
      startrun = raw_input("Start run: ")
      endrun = get_optional_input("End run (leave blank if the last run in this group hasn't been collected yet): ")
      detz_parameter = raw_input("Detz parameter (CXI: detz_offset, XPP: distance): ")
      beamx = get_optional_input("Beam center x (leave blank to not override): ")
      beamy = get_optional_input("Beam center y (leave blank to not override): ")
      pixelmask = get_optional_input("Path to untrusted pixel mask (if available): ")
      darkavg = get_optional_input("Path to dark average image (if available): ")
      darkstddev = get_optional_input("Path to dark standard deviation image (if available): ")
      gainmap = get_optional_input("Path to gain map image (if available): ")
      binning = get_optional_input("Binning (if applicable): ")
      #usecase = get_optional_input("Use case (list available cases here):")
      comment = raw_input("Add a comment: ")

      cmd = "INSERT INTO %s_rungroups (startrun,endrun,detz_parameter,beamx,beamy,untrusted_pixel_mask_path,dark_avg_path,dark_stddev_path,gain_map_path,binning,comment) VALUES (%s,%s,%s,%s,%s,'%s','%s','%s','%s',%s,'%s')"%(self.params.experiment_tag, startrun, endrun, detz_parameter, beamx, beamy, pixelmask, darkavg, darkstddev, gainmap, binning, comment)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Run group added")

  class update_rungroup(db_action):
    selection_string = "Update run group"
    help_string = "Change run group start and end values and comments"

    def __call__(self):
      rungroup = raw_input("Run group number: ")
      startrun = raw_input("New start run: ")
      endrun = get_optional_input("New end run (leave blank if the last run in this group hasn't been collected yet): ")
      comment = raw_input("New comment: ")

      cmd = "UPDATE %s_rungroups SET startrun=%s, endrun=%s, comment='%s' WHERE rungroup_id=%s"%(self.params.experiment_tag, startrun, endrun, comment, rungroup)
      cursor = self.dbobj.cursor()
      cursor.execute(cmd)
      self.dbobj.commit()
      print("Run group updated")

  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)
    self.menu_string = "Managing run groups"
    self.options = {
      'l': rungroups_menu.list_rungroups,
      'a': rungroups_menu.add_rungroup,
      'u': rungroups_menu.update_rungroup
    }

class top_menu(option_chooser):
  def __init__(self, params, dbobj):
    option_chooser.__init__(self, params, dbobj)

    self.menu_string = "Top level menu, administering %s using tag %s"%(self.params.experiment, self.params.experiment_tag)

    self.options = {
      'r': runs_menu,
      'i': isoforms_menu,
      't': trials_menu,
      'g': rungroups_menu,
      'l': link_trials_to_rungroup_menu
    }

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

  initialize(params, dbobj)()
  top_menu(params, dbobj)()

  print("Done")

if __name__ == "__main__":
  run(sys.argv[1:])

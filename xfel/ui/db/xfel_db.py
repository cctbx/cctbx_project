from __future__ import absolute_import, division, print_function

import os, time
import libtbx.load_env
from libtbx.utils import Sorry

from xfel.ui.db.trial import Trial
from xfel.ui.db.run import Run
from xfel.ui.db.rungroup import Rungroup
from xfel.ui.db.tag import Tag
from xfel.ui.db.job import Job
from xfel.ui.db.stats import Stats
from xfel.ui.db.experiment import Cell, Bin, Isoform, Event

from xfel.ui.db import get_db_connection
from six.moves import range
import six
from six.moves import zip

try:
  from MySQLdb import OperationalError
except ImportError:
  raise Sorry('Mysql not available')

from xfel.command_line.experiment_manager import initialize as initialize_base

class initialize(initialize_base):
  expected_tables = ["run", "job", "rungroup", "trial", "tag", "run_tag", "event", "trial_rungroup",
                     "imageset", "imageset_event", "beam", "detector", "experiment",
                     "crystal", "cell", "cell_bin", "bin", "isoform", "rungroup_run"]

  def __getattr__(self, prop):
    if prop == "dbobj":
      return get_db_connection(self.params)
    raise AttributeError("%s not found"%prop)

  def __setattr__(self, prop, val):
    if prop == "dbobj": pass
    return super(initialize, self).__setattr__(prop, val)

  def __init__(self, params, dbobj):
    initialize_base.__init__(self, params, dbobj, interactive = False, drop_tables = None)

  def create_tables(self, sql_path = None):
    if sql_path is None:
      sql_path = os.path.join(libtbx.env.find_in_repositories("xfel/ui/db"), "schema.sql")

    return initialize_base.create_tables(self, sql_path)

  def verify_tables(self):
    self.expected_tables.pop(self.expected_tables.index('rungroup_run'))
    tables_ok = super(initialize, self).verify_tables()
    self.expected_tables.append('rungroup_run')

    if tables_ok:
      # Maintain backwards compatibility with SQL tables v2: 09/24/16
      query = "SHOW columns FROM `%s_event`"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(query)
      columns = cursor.fetchall()
      column_names = zip(*columns)[0]
      if 'two_theta_low' not in column_names and 'two_theta_high' not in column_names:
        query = """
          ALTER TABLE %s_event
          ADD COLUMN two_theta_low DOUBLE NULL,
          ADD COLUMN two_theta_high DOUBLE NULL
        """%self.params.experiment_tag
        cursor.execute(query)
      elif 'two_theta_low' in column_names and 'two_theta_high' in column_names:
        pass
      else:
        assert False

      # Maintain backwards compatibility with SQL tables v2: 09/28/16
      query = "SHOW columns FROM `%s_job`"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(query)
      columns = cursor.fetchall()
      column_names = zip(*columns)[0]
      if 'submission_id' not in column_names:
        query = """
          ALTER TABLE %s_job
          ADD COLUMN submission_id VARCHAR(45) NULL
        """%self.params.experiment_tag
        cursor.execute(query)

      # Maintain backwards compatibility with SQL tables v2: 12/12/16
      query = "SHOW columns FROM `%s_rungroup`"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(query)
      columns = cursor.fetchall()
      column_names = zip(*columns)[0]
      for needed_column, column_format in zip(['format', 'two_theta_low', 'two_theta_high'],
                                              ["VARCHAR(45) NOT NULL DEFAULT 'pickle'",
                                               "DOUBLE NULL", "DOUBLE NULL"]):
        if needed_column not in column_names:
          query = """
            ALTER TABLE %s_rungroup
            ADD COLUMN %s %s
          """%(self.params.experiment_tag, needed_column, column_format)
          cursor.execute(query)

      # Maintain backwards compatibility with SQL tables v2: 06/23/17
      query = "SHOW columns FROM `%s_trial`"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(query)
      columns = cursor.fetchall()
      column_names = zip(*columns)[0]
      if 'd_min' not in column_names:
        query = """
          ALTER TABLE `%s_trial`
          ADD COLUMN d_min FLOAT NULL
        """%self.params.experiment_tag
        cursor.execute(query)
      query = "SHOW columns FROM `%s_cell`"%self.params.experiment_tag
      cursor = self.dbobj.cursor()
      cursor.execute(query)
      columns = cursor.fetchall()
      column_names = zip(*columns)[0]
      if 'trial_id' not in column_names:
        query = """
          ALTER TABLE `%s_cell`
          ADD COLUMN trial_id INT NULL,
          ADD CONSTRAINT fk_cell_trial1 FOREIGN KEY (trial_id) REFERENCES `%s_trial` (id) ON DELETE NO ACTION,
          ADD INDEX fk_cell_trial1_idx (trial_id ASC)
        """%(self.params.experiment_tag, self.params.experiment_tag)
        cursor.execute(query)

      # Maintain backwards compatibility with SQL tables v3: 02/1/19
      query = "SHOW TABLES LIKE '%s_rungroup_run'"%(self.params.experiment_tag)
      cursor.execute(query)
      if cursor.rowcount == 0:
        print("Upgrading to version 4 of mysql database schema")
        query = """
        CREATE TABLE IF NOT EXISTS `%s`.`%s_rungroup_run` (
          `rungroup_id` INT NOT NULL,
          `run_id` INT NOT NULL,
          PRIMARY KEY (`rungroup_id`, `run_id`),
          INDEX `fk_rungroup_has_run_run1_idx` (`run_id` ASC),
          INDEX `fk_rungroup_has_run_rungroup1_idx` (`rungroup_id` ASC),
          CONSTRAINT `fk_%s_rungroup_has_run_rungroup1`
            FOREIGN KEY (`rungroup_id`)
            REFERENCES `%s`.`%s_rungroup` (`id`)
            ON DELETE NO ACTION
            ON UPDATE NO ACTION,
          CONSTRAINT `fk_%s_rungroup_has_run_run1`
            FOREIGN KEY (`run_id`)
            REFERENCES `%s`.`%s_run` (`id`)
            ON DELETE NO ACTION
            ON UPDATE NO ACTION)
        ENGINE = InnoDB;
        """%(self.params.db.name, self.params.experiment_tag, self.params.experiment_tag,
             self.params.db.name, self.params.experiment_tag, self.params.experiment_tag,
             self.params.db.name, self.params.experiment_tag)
        cursor.execute(query)
        # Convert from startrun/endrun to linked table connecting run and rungroup and an 'open' flag in rungroup
        query = "ALTER TABLE `%s_rungroup` ADD COLUMN open TINYINT(1) NOT NULL DEFAULT 0"%self.params.experiment_tag
        cursor.execute(query)
        query = "SELECT id, startrun, endrun FROM `%s_rungroup`"%self.params.experiment_tag
        cursor.execute(query)
        for rungroup_id, startrun, endrun in cursor.fetchall():
          if endrun is None:
            # This run is 'open', so get the last run available and update the open bit for the rungroup
            query = 'SELECT run FROM `%s_run`'%self.params.experiment_tag
            cursor.execute(query)
            endrun = max(zip(*cursor.fetchall())[0])
            query = 'UPDATE `%s_rungroup` set open = 1 where id = %d'%(self.params.experiment_tag, rungroup_id)
            cursor.execute(query)
          # Add all thr runs to the rungroup
          for run in range(startrun, endrun+1):
            query = 'SELECT id FROM `%s_run` run WHERE run.run = %d'%(self.params.experiment_tag, run)
            cursor.execute(query)
            rows = cursor.fetchall(); assert len(rows) <= 1
            if len(rows) > 0:
              run_id = rows[0][0]
              query = "INSERT INTO `%s_rungroup_run` (rungroup_id, run_id) VALUES (%d, %d)" % ( \
                self.params.experiment_tag, rungroup_id, run_id)
              cursor.execute(query)
        self.dbobj.commit()
        query = "ALTER TABLE `%s_rungroup` DROP COLUMN startrun, DROP COLUMN endrun"%self.params.experiment_tag
        cursor.execute(query)
        # remove NOT NULL (and default for format column)
        query = "ALTER TABLE `%s_rungroup` MODIFY COLUMN format varchar(45)"%self.params.experiment_tag
        cursor.execute(query)
        query = "ALTER TABLE `%s_rungroup` MODIFY COLUMN detector_address varchar(100)"%self.params.experiment_tag
        cursor.execute(query)
        query = "ALTER TABLE `%s_rungroup` MODIFY COLUMN detz_parameter double"%self.params.experiment_tag
        cursor.execute(query)
        # Retype and add new columns
        query = "ALTER TABLE `%s_run` MODIFY COLUMN run varchar(45) NOT NULL"%self.params.experiment_tag
        cursor.execute(query)
        query = "ALTER TABLE `%s_run` ADD COLUMN path varchar(4097)"%self.params.experiment_tag
        cursor.execute(query)
        column_names = ['number', 'd_min', 'd_max', 'total_hkl']
        column_types = ['int', 'double', 'double', 'int']
        for column_name, column_type in zip(column_names, column_types):
          query = "ALTER TABLE `%s_bin` MODIFY COLUMN %s %s"%(self.params.experiment_tag, column_name, column_type)
          cursor.execute(query)
    return tables_ok

  def set_up_columns_dict(self, app):
    columns_dict = {}
    for table in self.expected_tables:
      table_name = "%s_%s" % (self.params.experiment_tag, table)
      query = "SHOW COLUMNS FROM `%s`" % (table_name)
      cursor = app.execute_query(query)
      columns_dict[table_name] = [c[0] for c in cursor.fetchall() if c[0] != 'id']
    return columns_dict

class db_application(object):
  def __init__(self, params):
    self.params = params

  def execute_query(self, query, commit = False):
    if self.params.db.verbose:
      from time import time
      st = time()
      self.query_count += 1

    retry_count = 0
    retry_max = 10
    sleep_time = 0.1
    while retry_count < retry_max:
      try:
        dbobj = get_db_connection(self.params)
        cursor = dbobj.cursor()
        cursor.execute(query)
        if commit:
          dbobj.commit()

        if self.params.db.verbose:
          et = time() - st
          if et > 1:
            print('Query % 6d SQLTime Taken = % 10.6f seconds' % (self.query_count, et), query[:min(len(query),145)])
        return cursor
      except OperationalError as e:
        if "Can't connect to MySQL server" not in str(e):
          print(query)
          raise e
        retry_count += 1
        print("Couldn't connect to MYSQL, retry", retry_count)
        time.sleep(sleep_time)
        sleep_time *= 2
      except Exception as e:
        print("Couldn't execute MYSQL query.  Query:")
        print(query)
        print("Exception:")
        print(str(e))
        raise e
    raise Sorry("Couldn't execute MYSQL query. Too many reconnects. Query: %s"%query)

class xfel_db_application(db_application):
  def __init__(self, params, drop_tables = False, verify_tables = False):
    super(xfel_db_application, self).__init__(params)
    self.query_count = 0
    dbobj = get_db_connection(params)
    self.init_tables = initialize(params, dbobj) # only place where a connection is held

    if drop_tables:
      self.drop_tables()

    if verify_tables and not self.verify_tables():
      self.create_tables()
      print('Creating experiment tables...')
      if not self.verify_tables():
        raise Sorry("Couldn't create experiment tables")

    self.columns_dict = self.init_tables.set_up_columns_dict(self)

  def list_lcls_runs(self):
    from xfel.xpp.simulate import file_table
    query = "https://pswww.slac.stanford.edu/ws-auth/dataexport/placed?exp_name=%s" % (self.params.facility.lcls.experiment)
    FT = file_table(self.params.facility.lcls, query, enforce80=self.params.facility.lcls.web.enforce80, enforce81=self.params.facility.lcls.web.enforce81)
    runs = FT.get_runs()
    for r in runs: r['run'] = str(r['run'])
    return runs

  def verify_tables(self):
    return self.init_tables.verify_tables()

  def create_tables(self):
    return self.init_tables.create_tables()

  def drop_tables(self):
    return self.init_tables.drop_tables()

  def create_trial(self, d_min = 1.5, n_bins = 10, **kwargs):
    # d_min and n_bins only used if isoforms are in this trial

    trial = Trial(self, d_min = d_min, **kwargs)
    if trial.target_phil_str is not None:
      from iotbx.phil import parse
      dispatcher = self.params.dispatcher
      if dispatcher == 'cxi.xtc_process':
        from spotfinder.applications.xfel import cxi_phil
        trial_params = cxi_phil.cxi_versioned_extract().persist.phil_scope.fetch(parse(trial.target_phil_str)).extract()
        isoforms = trial_params.isoforms
      else:
        from xfel.ui import load_phil_scope_from_dispatcher
        phil_scope = load_phil_scope_from_dispatcher(self.params.dispatcher)
        trial_params = phil_scope.fetch(parse(trial.target_phil_str)).extract()
        isoforms = trial_params.indexing.stills.isoforms
      if len(isoforms) > 0:
        for isoform in isoforms:
          print("Creating isoform", isoform.name)
          db_isoform = Isoform(self,
                               name = isoform.name,
                               trial_id = trial.id)
          a, b, c, alpha, beta, gamma = isoform.cell.parameters()
          cell = self.create_cell(cell_a = a, cell_b = b, cell_c = c,
                                  cell_alpha = alpha, cell_beta = beta, cell_gamma = gamma,
                                  lookup_symbol = isoform.lookup_symbol,
                                  isoform_id = db_isoform.id)
          from cctbx.crystal import symmetry

          cs = symmetry(unit_cell = isoform.cell,space_group_symbol=str(isoform.lookup_symbol))
          mset = cs.build_miller_set(anomalous_flag=False, d_min=d_min)
          binner = mset.setup_binner(n_bins=n_bins)
          for i in binner.range_used():
            d_max, d_min = binner.bin_d_range(i)
            Bin(self, number = i, d_min = d_min, d_max = d_max,
                total_hkl = binner.counts_complete()[i], cell_id = cell.id)
      elif dispatcher == 'cxi.xtc_process':
        pass # TODO: labelit target
      else:
        if trial_params.indexing.known_symmetry.unit_cell is not None and \
            trial_params.indexing.known_symmetry.space_group is not None:
          print("Creating target cell")
          unit_cell = trial_params.indexing.known_symmetry.unit_cell
          symbol = str(trial_params.indexing.known_symmetry.space_group)
          a, b, c, alpha, beta, gamma = unit_cell.parameters()
          cell = self.create_cell(cell_a = a, cell_b = b, cell_c = c,
                                  cell_alpha = alpha, cell_beta = beta, cell_gamma = gamma,
                                  lookup_symbol = symbol,
                                  trial_id = trial.id)
          from cctbx.crystal import symmetry

          cs = symmetry(unit_cell = unit_cell, space_group_symbol = symbol)
          mset = cs.build_miller_set(anomalous_flag=False, d_min=d_min)
          binner = mset.setup_binner(n_bins=n_bins)
          for i in binner.range_used():
            d_max, d_min = binner.bin_d_range(i)
            Bin(self, number = i, d_min = d_min, d_max = d_max,
                total_hkl = binner.counts_complete()[i], cell_id = cell.id)
    return trial

  def get_trial_isoforms(self, trial_id):
    where = "WHERE trial_id = %d"%trial_id
    return self.get_all_x(Isoform, "isoform", where)

  def get_trial_cell(self, trial_id):
    where = "WHERE trial_id = %d"%trial_id
    cells = self.get_all_x(Cell, "cell", where)
    assert len(cells) <= 1
    if len(cells) == 0:
      return None
    else:
      return cells[0]

  def get_trial_cells(self, trial_id, rungroup_id = None, run_id = None):
    # Use big queries to assist listing lots of cells. Start with list of cells for this trial
    tag = self.params.experiment_tag
    if rungroup_id is not None or run_id is not None:
      assert rungroup_id is not None and run_id is not None
      extra_where = "AND evt.run_id = %d AND evt.rungroup_id = %d"%(run_id, rungroup_id)
    else:
      extra_where = ""

    where  = """JOIN `%s_crystal` crystal ON crystal.cell_id = cell.id
               JOIN `%s_experiment` expt ON expt.crystal_id = crystal.id
               JOIN `%s_imageset` imgset ON imgset.id = expt.imageset_id
               JOIN `%s_imageset_event` is_e ON is_e.imageset_id = imgset.id
               JOIN `%s_event` evt ON evt.id = is_e.event_id
               JOIN `%s_trial` trial ON evt.trial_id = trial.id
               JOIN `%s_rungroup` rg ON evt.rungroup_id = rg.id
               WHERE trial.id = %d AND rg.active = True %s""" % (
               tag, tag, tag, tag, tag, tag, tag, trial_id, extra_where)
    cells = self.get_all_x(Cell, 'cell', where)
    where = " JOIN `%s_cell` cell ON bin.cell_id = cell.id "%(tag) + where
    return self.link_cell_bins(cells, where = where)

  def link_cell_bins(self, cells, where = None):
    tag = self.params.experiment_tag
    cells_d = {cell.id:cell for cell in cells}

    # Get all the bin ids for bins associated with these cells and assemble the bin objects
    if where is None:
      cell_ids = ", ".join([str(key) for key in cells_d.keys()])
      where = """ WHERE bin.cell_id IN (%s)""" % (cell_ids)
    bins = self.get_all_x(Bin, 'bin', where=where)

    # Link in all the bins
    for bin in bins:
      if bin.cell_id in cells_d: # might not be there if new data has arrived
        cells_d[bin.cell_id]._bins.append(bin)

    for cell in cells:
      cell._bins_set = True

    return cells

  def create_cell(self, **kwargs):
    return Cell(self, **kwargs)

  def get_cell(self, cell_id = None, name = None):
    assert [cell_id, name].count(None) == 1
    if name is not None:
      query = "SELECT id FROM `%s_cell` WHERE name = '%s'"%(self.params.experiment_tag, name)
      cursor = self.execute_query(query)
      results = cursor.fetchall()
      assert len(results) in [0,1]
      if len(results) == 0:
        return None
      cell_id = int(results[0][0])

    return Cell(self, cell_id=cell_id)

  def get_cell_bins(self, cell_id):
    query = "SELECT id FROM `%s_bin` WHERE cell_id = %d" % \
            (self.params.experiment_tag, cell_id)
    cursor = self.execute_query(query)
    ids = [str(i[0]) for i in cursor.fetchall()]
    if len(ids) == 0:
      return []
    where = "WHERE id IN (%s)" % ", ".join(ids)
    return self.get_all_x(Bin, 'bin', where)

  def get_all_x(self, cls, name, where = None):
    table_name = "%s_%s" % (self.params.experiment_tag, name)
    columns = ["%s.%s"%(name, c) for c in self.columns_dict[table_name]]
    query = "SELECT %s.id, %s FROM `%s` %s" % (name, ", ".join(columns), table_name, name)
    if where is not None:
      query += " " + where
    cursor = self.execute_query(query)
    results = []
    for row in cursor.fetchall():
      d = {key:value for key, value in zip(self.columns_dict[table_name], row[1:])}
      d["%s_id"%name] = row[0]
      results.append(cls(self, **d))
    return results

  def get_all_x_with_subitems(self, cls, name, where = None, sub_items = None):
    """ Assemble a list of db_proxy objects, where each one references a sub object
        @param cls principal class
        @param name table name not including experiment tag
        @param where optional constraints on what to get
        @param sub_items: array of tuples (cls,name) of items this table references
        @return list of assembled db_proxy objects"""

    # Assemble a list of all columns to be extracted via a big join. Then entries in the
    # columns array will be of the form "name.column", IE: tag.comment
    table_name = "%s_%s" % (self.params.experiment_tag, name)
    columns = ["%s.%s"%(name, c) for c in self.columns_dict[table_name]]
    columns.append("%s.id"%name)

    if sub_items is None:
      sub_items = []
    sub_table_names = ["%s_%s"%(self.params.experiment_tag, i[1]) for i in sub_items]
    for i, sub_item in enumerate(sub_items):
      scls, sname = sub_item
      columns.extend(["%s.%s"%(sname, c) for c in self.columns_dict[sub_table_names[i]]])
      columns.append("%s.id"%sname)

    # the main item being extracted is in the FROM statement and is given a nickname which
    # is the table name without the experiment tag
    query = "SELECT %s FROM `%s` %s" % (", ".join(columns), table_name, name)

    # Join statements to bring in the sub tables
    for i, sub_item in enumerate(sub_items):
      scls, sname = sub_item
      query += " JOIN `%s` %s ON %s.id = %s.%s_id"% (
        sub_table_names[i], sname, sname, name, sname)

    if where is not None:
      query += " " + where
    cursor = self.execute_query(query)

    results = []
    for row in cursor.fetchall():
      # Each row will be a complete item and sub items in column form. Assemble one
      # dictionary (d) for the main item and a dictionary of dictionaries (sub_ds)
      # for each of the sub items
      d = {}
      sub_ds = {sub_item[1]:(sub_item[0], {}) for sub_item in sub_items}
      for key, value in zip(columns, row):
        n, c = key.split('.') # nickname n, column name c
        if n == name:
          d[c] = value # this column came from the main table
        else:
          sub_ds[n][1][c] = value # this column came from a sub table

      # pop the id column as it is passed as name_id to the db_proxy class (ie Job(job_id = 2))
      _id = d.pop("id")
      d["%s_id"%name] = _id
      results.append(cls(self, **d)) # instantiate the main class
      for sub_d_n, sub_d in six.iteritems(sub_ds):
        _id = sub_d[1].pop("id")
        sub_d[1]["%s_id"%sub_d_n] = _id
        setattr(results[-1], sub_d_n, sub_d[0](self, **sub_d[1])) # instantiate the sub items
    return results

  def get_trial(self, trial_id = None, trial_number = None):
    assert [trial_id, trial_number].count(None) == 1
    if trial_id is None:
      trials = [t for t in self.get_all_trials() if t.trial == trial_number]
      assert len(trials) == 1
      return trials[0]
    else:
      return Trial(self, trial_id)

  def get_trial_rungroups(self, trial_id, only_active = False):
    tag = self.params.experiment_tag
    if only_active:
      query = """SELECT t_rg.rungroup_id FROM `%s_trial_rungroup` t_rg
                 JOIN `%s_rungroup` rg ON rg.id = t_rg.rungroup_id
                 WHERE t_rg.trial_id = %d AND rg.active = True""" % (tag, tag, trial_id)
    else:
      query = """SELECT t_rg.rungroup_id FROM `%s_trial_rungroup` t_rg
                 WHERE t_rg.trial_id = %d""" % (tag, trial_id)
    cursor = self.execute_query(query)
    rungroup_ids = ["%d"%i[0] for i in cursor.fetchall()]
    if len(rungroup_ids) == 0:
      return []
    return self.get_all_x(Rungroup, "rungroup", where = "WHERE rungroup.id IN (%s)"%",".join(rungroup_ids))

  def get_trial_runs(self, trial_id):
    tag = self.params.experiment_tag
    query = """SELECT run.id FROM `%s_run` run
               JOIN `%s_rungroup_run` rgr ON run.id = rgr.run_id
               JOIN `%s_rungroup` rg ON rg.id = rgr.rungroup_id
               JOIN `%s_trial_rungroup` t_rg on t_rg.rungroup_id = rg.id
               JOIN `%s_trial` trial ON trial.id = t_rg.trial_id
               WHERE trial.id=%d AND rg.active=True
               """ %(tag, tag, tag, tag, tag, trial_id)
    cursor = self.execute_query(query)
    run_ids = ["%d"%i[0] for i in cursor.fetchall()]
    if len(run_ids) == 0:
      return []
    return self.get_all_x(Run, "run", where = "WHERE run.id IN (%s)"%",".join(run_ids))

  def get_trial_tags(self, trial_id):
    tag = self.params.experiment_tag
    query = """SELECT tag.id FROM `%s_tag` tag
               JOIN `%s_run_tag` rt ON rt.tag_id = tag.id
               JOIN `%s_run` run ON run.id = rt.run_id
               JOIN `%s_rungroup_run` rgr ON run.id = rgr.run_id
               JOIN `%s_rungroup` rg ON rg.id = rgr.rungroup_id
               JOIN `%s_trial_rungroup` t_rg ON t_rg.rungroup_id = rg.id
               JOIN `%s_trial` trial ON trial.id = t_rg.trial_id
               WHERE trial.id=%d AND rg.active=True
               """ % (tag, tag, tag, tag, tag, tag, tag, trial_id)
    cursor = self.execute_query(query)
    tag_ids = ["%d"%i[0] for i in cursor.fetchall()]
    if len(tag_ids) == 0:
      return []
    return self.get_all_x(Tag, "tag", where = "WHERE tag.id IN (%s)"%",".join(tag_ids))

  def get_all_trials(self, only_active = False):
    if only_active:
      return [t for t in self.get_all_x(Trial, "trial") if t.active]
    else:
      return self.get_all_x(Trial, "trial")

  def create_run(self, **kwargs):
    return Run(self, **kwargs)

  def get_run(self, run_id = None, run_number = None):
    assert [run_id, run_number].count(None) == 1
    if run_id is not None:
      return Run(self, run_id)
    runs = [r for r in self.get_all_runs() if r.run == str(run_number)]
    if len(runs) == 0:
      raise Sorry("Couldn't find run %s"%run_number)
    assert len(runs) == 1
    return runs[0]

  def get_all_runs(self):
    return self.get_all_x(Run, "run")

  def create_rungroup(self, **kwargs):
    return Rungroup(self, **kwargs)

  def get_rungroup(self, rungroup_id):
    return Rungroup(self, rungroup_id)

  def get_rungroup_runs(self, rungroup_id):
    query = """SELECT run.id FROM `%s_run` run
               JOIN `%s_rungroup_run` rgr ON run.id = rgr.run_id
               WHERE rgr.rungroup_id = %d"""%(self.params.experiment_tag,
                 self.params.experiment_tag, rungroup_id)
    cursor = self.execute_query(query)
    run_ids = ["%d"%i[0] for i in cursor.fetchall()]
    if len(run_ids) == 0:
      return []
    return self.get_all_x(Run, "run", where = "WHERE run.id IN (%s)"%",".join(run_ids))

  def get_all_rungroups(self, only_active = False):
    if only_active:
      return [rg for rg in self.get_all_x(Rungroup, "rungroup") if rg.active]
    else:
      return self.get_all_x(Rungroup, "rungroup")

  def create_tag(self, **kwargs):
    return Tag(self, **kwargs)

  def get_tag(self, tag_id):
    return Tag(self, tag_id)

  def get_run_tags(self, run_id):
    query = "SELECT tag_id from `%s_run_tag` WHERE `%s_run_tag`.run_id = %d" % \
            (self.params.experiment_tag, self.params.experiment_tag, run_id)
    cursor = self.execute_query(query)
    tag_ids = [str(i[0]) for i in cursor.fetchall()]
    if len(tag_ids) == 0:
      return []
    where = "WHERE id IN (%s)" % ", ".join(tag_ids)
    return self.get_all_x(Tag, 'tag', where)

  def get_all_tags(self):
    return self.get_all_x(Tag, "tag")

  def delete_x(self, item, item_id):
    query = "DELETE FROM `%s` WHERE id = %d"%(item.table_name, item_id)
    self.execute_query(query, commit = True)

  def delete_tag(self, tag = None, tag_id = None):
    assert [tag, tag_id].count(None) == 1

    if tag_id is None:
      tag_id = tag.tag_id
    else:
      tag = self.get_tag(tag_id)

    query = "DELETE FROM `%s_run_tag` WHERE tag_id = %d" % (self.params.experiment_tag, tag_id)
    self.execute_query(query, commit=True)

    self.delete_x(tag, tag_id)

  def create_job(self, **kwargs):
    return Job(self, **kwargs)

  def get_job(self, job_id):
    return Job(self, job_id)

  def get_all_jobs(self, active = False):
    if active:
      where = "WHERE trial.active = True AND rungroup.active = True"
    else:
      where = None
    return self.get_all_x_with_subitems(Job, "job", sub_items = [(Trial, 'trial'), (Run, 'run'), (Rungroup, 'rungroup')], where = where)

  def delete_job(self, job = None, job_id = None):
    assert [job, job_id].count(None) == 1
    if job_id is None:
      job_id = job.job_id
    else:
      job = self.get_job(job_id)

    self.delete_x(job, job_id)

  def get_all_events(self, trial = None, runs = None, only_indexed = True, isoform = None, where = None):
    tag = self.params.experiment_tag
    if where is None:
      final_where = ""
    else:
      final_where = where.strip()
    where = ""
    if only_indexed or isoform is not None:
      where = " JOIN `%s_imageset_event` is_e ON event.id = is_e.event_id"%tag
    if trial is not None:
      if runs is None:
        runs = trial.runs
      if len(runs) == 0:
        return []
      if isoform is None:
        where += " WHERE"
      else:
        where += """ JOIN `%s_imageset` imgset ON imgset.id = is_e.imageset_id
                     JOIN `%s_experiment` exp ON exp.imageset_id = imgset.id
                     JOIN `%s_crystal` crystal ON crystal.id = exp.crystal_id
                     JOIN `%s_cell` cell ON cell.id = crystal.cell_id
                     JOIN `%s_isoform` isoform ON isoform.id = cell.isoform_id
                     WHERE isoform.id = %s AND"""%(tag, tag, tag, tag, tag, isoform.id)
      where += " event.trial_id = %d AND event.run_id in (%s)" % (
        trial.id, ", ".join([str(r.id) for r in runs]))

      if 'rungroup_id' in self.columns_dict["%s_%s" % (tag, 'event')]: # some backwards compatibility, as event.rungroup_id was added late to the schema
        rungroups = ", ".join([str(rg.id) for rg in trial.rungroups])
        if len(rungroups) > 0:
          where += " AND event.rungroup_id in (%s)"%rungroups

    if len(where) > 0 and len(final_where) > 0:
      final_where = "AND " + final_where.lstrip("WHERE")
    return self.get_all_x(Event, "event", where + " " + final_where)

  def get_stats(self, **kwargs):
    return Stats(self, **kwargs)

class standalone_run_finder(object):
  def __init__(self, params):
    self.params = params

  def list_standalone_runs(self):
    import glob

    runs = []
    if self.params.facility.standalone.monitor_for == 'folders':
      for foldername in sorted(os.listdir(self.params.facility.standalone.data_dir)):
        path = os.path.join(self.params.facility.standalone.data_dir, foldername)
        if not os.path.isdir(path): continue
        if self.params.facility.standalone.composite_files:
          for filepath in sorted(glob.glob(os.path.join(path, self.params.facility.standalone.template))):
            filename = os.path.basename(filepath)
            runs.append((foldername + "_" + os.path.splitext(filename)[0], filepath))
        else:
          files = sorted(glob.glob(os.path.join(path, self.params.facility.standalone.template)))
          if len(files) > 0:
            runs.append((foldername, os.path.join(path, self.params.facility.standalone.template)))
    elif self.params.facility.standalone.monitor_for == 'files':
      if not self.params.facility.standalone.composite_files:
        print("Warning, monitoring a single folder for single image files is inefficient")
      path = self.params.facility.standalone.data_dir
      for filepath in sorted(glob.glob(os.path.join(path, self.params.facility.standalone.template))):
        filename = os.path.basename(filepath)
        runs.append((os.path.splitext(filename)[0], filepath))

    return runs

class cheetah_run_finder(standalone_run_finder):
  def __init__(self, params):
    super(cheetah_run_finder, self).__init__(params)
    self.known_runs = []
    self.known_bad_runs = []

  def list_standalone_runs(self):
    runs = super(cheetah_run_finder, self).list_standalone_runs()

    good_runs = []
    for name, path in runs:
      if name in self.known_runs:
        good_runs.append((name, path))
      elif name not in self.known_bad_runs:
        status_file = os.path.join(os.path.dirname(path), 'status.txt')
        if not os.path.exists(status_file): continue
        status_lines = [line for line in open(status_file).readlines() if 'Status' in line]
        if len(status_lines) != 1: continue
        l = status_lines[0].strip()
        status = l.split(',')[-1].split('=')[-1]
        print(name, status)
        if status == 'Finished':
          good_runs.append((name, path))
          self.known_runs.append(name)
        elif 'error' in status.lower():
          self.known_bad_runs.append(name)

    return good_runs

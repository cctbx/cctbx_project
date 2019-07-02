from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex

from cctbx.miller import set as mset

from cctbx.uctbx import unit_cell
from cctbx.crystal import symmetry
import time

from libtbx.development.timers import Timer

def phil_validation(params):
  return True

def application(params, loop = True):
  from cxi_xdr_xes.cftbx.cspad_ana import db as cxidb
  dbobj = cxidb.dbconnect(params.db.host, params.db.name, params.db.user, params.db.password)
  cursor = dbobj.cursor()
  PM = progress_manager(params, cursor)
  PM.setup_isoforms(cursor)
  PM.setup_runtags(cursor)
  isoforms = PM.isoforms
  del dbobj

  while 1:
    dbobj = cxidb.dbconnect(params.db.host, params.db.name, params.db.user, params.db.password)
    cursor = dbobj.cursor()

    results = {}
    print("Looking for data...")

    for tag in params.run_tags.split(','):
      for isoform in isoforms:
        M = PM.get_HKL(cursor,isoform=isoform,run_tags=tag)
        cell = isoforms[isoform]['cell']
        miller_set = mset(anomalous_flag = False, crystal_symmetry=symmetry(unit_cell=cell, space_group_symbol=isoforms[isoform]['lookup_symbol']), indices=M)
        miller_set.show_comprehensive_summary()

        miller_set.setup_binner(d_min=params.resolution, n_bins=params.n_bins)
        given = miller_set.binner().counts_given()
        ccomplete = miller_set.binner().counts_complete()
        for i_bin in miller_set.binner().range_used():
            sel         = miller_set.binner().selection(i_bin)
            self_sel    = miller_set.select(sel)
            d_max,d_min = self_sel.d_max_min()
            compl       = self_sel.completeness(d_max = d_max)

            n_ref       = sel.count(True)
            if ccomplete[i_bin] == 0.:
              multiplicity = 0.
            else:
              res_highest   = d_min
              multiplicity  = given[i_bin]/ccomplete[i_bin]
            d_range     = miller_set.binner().bin_legend(
                   i_bin = i_bin, show_bin_number = False, show_counts = True)
            fmt = "%3d: %-24s %4.2f %6d mult=%4.2f"
            print(fmt % (i_bin,d_range,compl,n_ref,
                          multiplicity))
        print()
        if len(tag) > 0:
          key = "%s %s"%(tag, isoform)
        else:
          key = isoform
        given_used = flex.int(given).select(flex.size_t(miller_set.binner().range_used()))
        ccomplete_used = flex.int(ccomplete).select(flex.size_t(miller_set.binner().range_used()))
        results[key] = dict(
          multiplicity = flex.sum(given_used)/flex.sum(ccomplete_used),
          completeness = miller_set.completeness(),
          multiplicity_highest = multiplicity,
          completeness_highest = compl,
          resolution_highest = res_highest
        )
    del dbobj
    if not loop:
      return results
    time.sleep(10)

from xfel.merging.database.merging_database import manager
class progress_manager(manager):
  def __init__(self,params, cursor):
    self.params = params
    if params.experiment_tag is None:
      self.db_experiment_tag = params.experiment
    else:
      self.db_experiment_tag = params.experiment_tag

    query = "SELECT trial_id FROM %s_trials WHERE %s_trials.trial = %d"%(self.db_experiment_tag, self.db_experiment_tag, params.trial)
    cursor.execute(query)
    assert cursor.rowcount == 1
    self.trial_id = int(cursor.fetchall()[0][0])

  def setup_isoforms(self, cursor):
    isoform_ids = []
    query = "SELECT DISTINCT(isoforms_isoform_id) from %s_frames WHERE %s_frames.trials_id = %d"%(self.db_experiment_tag, self.db_experiment_tag, self.trial_id)
    cursor.execute(query)
    for entry in cursor.fetchall():
      isoform_ids.append(int(entry[0]))

    d = {}
    for isoform_id in isoform_ids:
      query = "SELECT * FROM %s_isoforms WHERE %s_isoforms.isoform_id = %d"%(self.db_experiment_tag, self.db_experiment_tag, isoform_id)
      cursor.execute(query)
      assert cursor.rowcount == 1
      if_id, name, cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma, lookup_symbol = cursor.fetchall()[0]
      assert isoform_id == int(if_id)
      cell = unit_cell((cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma))
      d[name] = dict(
        isoform_id = isoform_id,
        cell = cell,
        lookup_symbol = lookup_symbol)
    self.isoforms = d

  def setup_runtags(self, cursor):
    if self.params.run_tags is None:
      self.params.run_tags = ""

    if len(self.params.run_tags) > 0:
      return

    all_tags = []
    name = self.db_experiment_tag
    query = """ SELECT rg.startrun, rg.endrun
                  FROM %s_rungroups rg
                  JOIN %s_trial_rungroups trg
                    ON rg.rungroup_id = trg.rungroups_id
                    AND trg.trials_id = %d
                    AND trg.active
               """%(name, name, self.trial_id)
    cursor.execute(query)

    for entry in cursor.fetchall():
      startrun, endrun = entry
      query = """SELECT DISTINCT(runs.tags)
                   FROM %s_runs runs
                   WHERE runs.run_id >= %d
              """%(name, startrun)
      if endrun is not None:
        query += " AND runs.run_id <= %d"%endrun
      cursor.execute(query)

      for row in cursor.fetchall():
        tags = row[0]
        if tags is None: continue
        if len(tags) == 0: continue
        all_tags.extend(tags.split(','))
    self.params.run_tags = ','.join(set(all_tags))

    print("Found these run tags:", self.params.run_tags)

  def get_HKL(self,cursor,isoform,run_tags):
    name = self.db_experiment_tag
    if run_tags is not None:
      extrajoin = "JOIN %s_runs runs ON frames.runs_run_id = runs.run_id"%name
      for tag in run_tags.split():
        tag = tag.strip()
        extrajoin += " AND runs.tags LIKE '%%%s%%'"%tag
    else:
      extrajoin = ""
    if self.params.include_negatives:
      extrawhere = ""
    else:
      extrawhere = "WHERE obs.i >= 0"

    query = """SELECT hkls.h, hkls.k, hkls.l
               FROM %s_observations obs
               JOIN %s_hkls hkls ON obs.hkls_id = hkls.hkl_id
               JOIN %s_isoforms isos ON hkls.isoforms_isoform_id = isos.isoform_id
                 AND isos.isoform_id = %d
               JOIN %s_frames frames ON obs.frames_id = frames.frame_id
                 AND frames.trials_id = %d
               JOIN %s_trial_rungroups trg
                 ON trg.trials_id = frames.trials_id
                 AND trg.rungroups_id = frames.rungroups_id
                 AND trg.active
               %s
               %s"""%(
               name, name, name, self.isoforms[isoform]['isoform_id'], name, self.trial_id, name, extrajoin, extrawhere)

    #print query
    if len(run_tags) > 0:
      print("%s, isoform %s"%(run_tags, isoform))
    else:
      print("isoform %s"%isoform)
    T = Timer("Reading db...")
    cursor.execute(query)
    del T
    T = Timer("Getting results...")
    ALL = cursor.fetchall()
    del T
    T = Timer("Parsing results...")
    indices = flex.miller_index([(a[0],a[1],a[2]) for a in ALL])
    del T

    print("ISOFORM %s:"%(isoform),len(indices),"result")
    return indices

  def connection(self):
    pass

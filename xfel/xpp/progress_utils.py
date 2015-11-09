from __future__ import division
from cctbx.array_family import flex

from cctbx.miller import set as mset

from cctbx.uctbx import unit_cell
from cctbx.crystal import symmetry
import time


def phil_validation(params):
  return True

def application(params):
  PM = progress_manager(params)
  from cxi_xdr_xes.cftbx.cspad_ana import db as cxidb

  while 1:
    dbobj = cxidb.dbconnect(params.db.host, params.data, params.db.user, params.db.password)
    cursor = dbobj.cursor()
    for isoform in ["A","B"]:
      M = PM.get_HKL(cursor,isoform=isoform)
      cell = dict(B=unit_cell((118.5,225.1,313.3)), A=unit_cell((118.4,225.6,333.2,90,90,90)))[isoform]
      miller_set = mset(anomalous_flag = False, crystal_symmetry=symmetry(unit_cell=cell, space_group_symbol="P m m m"), indices=M)
      miller_set.show_comprehensive_summary()

      miller_set.setup_binner(d_min=2.5, n_bins=15)
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
            multiplicity= given[i_bin]/ccomplete[i_bin]
          d_range     = miller_set.binner().bin_legend(
                 i_bin = i_bin, show_bin_number = False, show_counts = True)
          fmt = "%3d: %-24s %4.2f %6d mult=%4.2f"
          print fmt % (i_bin,d_range,compl,n_ref,
                        multiplicity)
      print
    del dbobj
    time.sleep(10)

from xfel.cxi.merging_database import manager
class progress_manager(manager):
  def __init__(self,params):
    self.params = params
    self.db_name = params.data

  def get_HKL(self,cursor,isoform):
    name = self.db_name
    query = '''select %s_hkls.h,%s_hkls.k,%s_hkls.l
               from %s_observations,%s_isoforms,%s_hkls,%s_trials
               WHERE %s_trials.trial = %s
               AND %s_trials.trial_id = %s_observations.frames_trials_id
               AND %s_hkls.isoforms_id = %s_isoforms.isoform_id
               AND %s_isoforms.name = '%s'
               AND %s_observations.hkls_id = %s_hkls.hkl_id
               ORDER BY %s_hkls.h,%s_hkls.k,%s_hkls.l'''%(
            name, name, name, name, name, name, name, name,
            self.params.trial,
            name, name, name, name, name,
            isoform,
            name, name, name, name, name)
    #print query
    cursor.execute(query)
    ALL = cursor.fetchall()
    indices = flex.miller_index([(a[0],a[1],a[2]) for a in ALL])

    print "ISOFORM %s:"%(isoform),len(indices),"result"
    return indices

  def connection(self):
    pass

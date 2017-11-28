# LIBTBX_SET_DISPATCHER_NAME mpi.merge
from __future__ import division
import sys,time
from xfel.command_line import cxi_merge

from xfel.command_line.single_node_merge import get_observations
cxi_merge.get_observations = get_observations

from xfel.command_line.single_node_merge import scaling_manager as scaling_manager_base
class scaling_manager_mpi(scaling_manager_base):

  def scale_all (self, file_names) :
    tar_list,integration_pickle_names = file_names
    t1 = time.time()
    if self.params.backend == 'MySQL':
      from xfel.cxi.merging_database import manager
    elif self.params.backend == 'SQLite':
      from xfel.cxi.merging_database_sqlite3 import manager
    else:
      from xfel.cxi.merging_database_fs import manager

    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())

    # MPI
    nproc = self.params.nproc
    assert nproc >= 1
    if nproc > 1:
      self._scale_all_parallel(tar_list, db_mgr)
    else:
      self._scale_all_serial(tar_list, db_mgr)
    # done
    db_mgr.join()

    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(integration_pickle_names))
    print >> self.log, "  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais
    print >> self.log, "  %d rejected for unit cell outliers" % \
      self.n_wrong_cell
    print >> self.log, "  %d rejected for low signal" % \
      self.n_low_signal
    print >> self.log, "  %d rejected due to up-front poor correlation under min_corr parameter" % \
      self.n_low_corr
    print >> self.log, "  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error
    for key in self.failure_modes.keys():
      print >>self.log, "  %d rejected due to %s"%(self.failure_modes[key], key)

    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_corr + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell \
               + sum([val for val in self.failure_modes.itervalues()])
    assert checksum == len(integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print >> self.log, "Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min)

    if self.params.raw_data.sdfac_refine:
      self.scale_errors()

    if self.params.raw_data.errors_from_sample_residuals:
      self.errors_from_residuals()

  def _scale_all_parallel (self, file_names, db_mgr) :
    from mpi4py import MPI
    comm = MPI.COMM_SELF.Spawn("mpi.worker2", args=[], #libtbx.python",args=['worker2.py'],
                           maxprocs=self.params.nproc)
    transmitted_info = dict(file_names=file_names,
                            miller_set=self.miller_set,
                            model = self.i_model,
                            params = self.params )
    comm.bcast(transmitted_info, root = MPI.ROOT)
    reports = comm.gather("no data",root=MPI.ROOT)
    print "reports",reports
    comm.Disconnect()
    for worker_sm in reports:
      self._add_all_frames(worker_sm)

  def _scale_all_serial (self, tar_list, db_mgr) :
    """
    Scale frames sequentially (single-process).  The return value is
    picked up by the callback.
    """
    self.tar_to_scale_frame_adapter(tar_list, db_mgr)
    return (self)

cxi_merge.scaling_manager = scaling_manager_mpi

if (__name__ == "__main__"):
  result = cxi_merge.run(args=sys.argv[1:])
  if result is None:
    sys.exit(1)

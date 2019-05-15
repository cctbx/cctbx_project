# LIBTBX_SET_DISPATCHER_NAME dev.mpi.cluster_two_merge
from __future__ import absolute_import, division, print_function
from six.moves import range
import sys,time

from xfel.command_line.cxi_merge import scaling_manager as scaling_manager_base
import six
class scaling_manager_mpi(scaling_manager_base):

  def mpi_initialize (self, file_names) :
    self.integration_pickle_names = file_names
    self.t1 = time.time()

    assert self.params.backend == 'FS' # for prototyping rank=0 marshalling
    from xfel.merging.database.merging_database_fs import manager
    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())
    self.master_db_mgr = db_mgr

  def mpi_finalize (self) :
    t2 = time.time()
    print("", file=self.log)
    print("#" * 80, file=self.log)
    print("FINISHED MERGING", file=self.log)
    print("  Elapsed time: %.1fs" % (t2 - self.t1), file=self.log)
    print("  %d of %d integration files were accepted" % (
      self.n_accepted, len(self.integration_pickle_names)), file=self.log)
    print("  %d rejected due to wrong Bravais group" % \
      self.n_wrong_bravais, file=self.log)
    print("  %d rejected for unit cell outliers" % \
      self.n_wrong_cell, file=self.log)
    print("  %d rejected for low signal" % \
      self.n_low_signal, file=self.log)
    print("  %d rejected due to up-front poor correlation under min_corr parameter" % \
      self.n_low_corr, file=self.log)
    print("  %d rejected for file errors or no reindex matrix" % \
      self.n_file_error, file=self.log)
    for key in self.failure_modes.keys():
      print("  %d rejected due to %s"%(self.failure_modes[key], key), file=self.log)

    checksum = self.n_accepted  + self.n_file_error \
               + self.n_low_corr + self.n_low_signal \
               + self.n_wrong_bravais + self.n_wrong_cell \
               + sum([val for val in six.itervalues(self.failure_modes)])
    assert checksum == len(self.integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print("Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min), file=self.log)

    if self.params.raw_data.sdfac_refine or self.params.raw_data.errors_from_sample_residuals:
      if self.params.raw_data.sdfac_refine:
        if self.params.raw_data.error_models.sdfac_refine.minimizer == 'simplex':
          from xfel.merging.algorithms.error_model.sdfac_refine import sdfac_refine as error_modeler
        elif self.params.raw_data.error_models.sdfac_refine.minimizer == 'lbfgs':
          from xfel.merging.algorithms.error_model.sdfac_refine_lbfgs import sdfac_refine_refltable_lbfgs as error_modeler
        elif self.params.raw_data.error_models.sdfac_refine.minimizer == 'LevMar':
          from xfel.merging.algorithms.error_model.sdfac_refine_levmar import sdfac_refine_refltable_levmar as error_modeler

      if self.params.raw_data.errors_from_sample_residuals:
        from xfel.merging.algorithms.error_model.errors_from_residuals import errors_from_residuals as error_modeler

      error_modeler(self).adjust_errors()

from xfel.merging.command_line.dev_cxi_merge import Script as base_Script

class Script(base_Script):

  def validate(self):
    base_Script.validate(self)
    if (self.params.rescale_with_average_cell):
      raise Usage("""Rescaling_with_average_cell not supported with MPI
      (Would require a second round of scaling, inefficient).""")

  def run(self,comm,timing=False):
    rank = comm.Get_rank()
    size = comm.Get_size()
    from time import time as tt

    # set things up
    if rank == 0:
      if timing: print("SETUP START RANK=%d TIME=%f"%(rank,tt()))
      self.initialize()
      self.validate()
      self.read_models()
      scaler_master = self.scaler_class(
        miller_set=self.miller_set,
        i_model=self.i_model,
        params=self.params,
        log=self.out)
      scaler_master.mpi_initialize(self.frame_files)

      transmitted_info = dict(file_names=self.frame_files,
                              miller_set=self.miller_set,
                              model = self.i_model,
                              params = self.params )
      if timing: print("SETUP END RANK=%d TIME=%f"%(rank,tt()))

    else:
      if timing: print("SETUP START RANK=%d TIME=%f"%(rank,tt()))
      transmitted_info = None
      if timing: print("SETUP END RANK=%d TIME=%f"%(rank,tt()))

    if timing: print("BROADCAST START RANK=%d TIME=%f"%(rank,tt()))
    transmitted_info = comm.bcast(transmitted_info, root = 0)
    if timing: print("BROADCAST END RANK=%d TIME=%f"%(rank,tt()))

    # now actually do the work
    if timing: print("SCALER_WORKER_SETUP START RANK=%d TIME=%f"%(rank,tt()))
    scaler_worker = self.scaler_class(transmitted_info["miller_set"],
                                       transmitted_info["model"],
                                       transmitted_info["params"],
                                       log = sys.stdout)
    if timing: print("SCALER_WORKER_SETUP END RANK=%d TIME=%f"%(rank,tt()))
    assert scaler_worker.params.backend == 'FS' # only option that makes sense
    from xfel.merging.database.merging_database_fs import manager2 as manager
    db_mgr = manager(scaler_worker.params)
    file_names = [transmitted_info["file_names"][i] for i in range(len(transmitted_info["file_names"])) if i%size == rank]
    if timing: print("SCALER_WORKERS START RANK=%d TIME=%f"%(rank, tt()))
    scaler_worker._scale_all_serial(file_names, db_mgr)
    if timing: print("SCALER_WORKERS END RANK=%d TIME=%f"%(rank, tt()))
    scaler_worker.finished_db_mgr = db_mgr
    # might want to clean up a bit before returning
    del scaler_worker.log
    del scaler_worker.params
    del scaler_worker.miller_set
    del scaler_worker.i_model
    del scaler_worker.reverse_lookup

    # gather reports and all add together
    if timing: print("GATHER START RANK=%d TIME=%f"%(rank, tt()))
    reports = comm.gather(scaler_worker,root=0)
    if timing: print("GATHER END RANK=%d TIME=%f"%(rank, tt()))
    if rank == 0:
      print("Processing reports from %d ranks"%(len(reports)))
      ireport = 0
      for item in reports:
        if timing: print("SCALER_MASTER_ADD START RANK=%d TIME=%f"%(rank, tt()))
        scaler_master._add_all_frames(item)

        if timing: print("SCALER_MASTER_ADD END RANK=%d TIME=%f"%(rank, tt()))
        print("processing %d calls from report %d"%(len(item.finished_db_mgr.sequencer),ireport)); ireport += 1

        for call_instance in item.finished_db_mgr.sequencer:
          if call_instance["call"] == "insert_frame":
            if timing: print("SCALER_MASTER_INSERT_FRAME START RANK=%d TIME=%f"%(rank, tt()))
            frame_id_zero_base = scaler_master.master_db_mgr.insert_frame(**call_instance["data"])
            if timing: print("SCALER_MASTER_INSERT_FRAME END RANK=%d TIME=%f"%(rank, tt()))
          elif call_instance["call"] == "insert_observation":
            if timing: print("SCALER_MASTER_INSERT_OBS START RANK=%d TIME=%f"%(rank, tt()))
            call_instance["data"]['frame_id_0_base'] = [frame_id_zero_base] * len(call_instance["data"]['frame_id_0_base'])
            scaler_master.master_db_mgr.insert_observation(**call_instance["data"])
            if timing: print("SCALER_MASTER_INSERT_OBS END RANK=%d TIME=%f"%(rank, tt()))

      if timing: print("SCALER_MASTER_FINALISE START RANK=%d TIME=%f"%(rank, tt()))
      scaler_master.master_db_mgr.join() # database written, finalize the manager
      scaler_master.mpi_finalize()
      if timing: print("SCALER_MASTER_FINALISE END RANK=%d TIME=%f"%(rank, tt()))

      return self.finalize(scaler_master)

if (__name__ == "__main__"):
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()

  script = Script(scaling_manager_mpi)
  result = script.run(comm=comm,timing=False)
  if rank == 0:
    script.show_plot(result)
  print("DONE")

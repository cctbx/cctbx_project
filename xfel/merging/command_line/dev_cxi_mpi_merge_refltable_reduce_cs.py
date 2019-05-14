from __future__ import absolute_import, division, print_function
from six.moves import range
import sys

from time import time as tt
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# LIBTBX_SET_DISPATCHER_NAME dev.cxi.mpi_merge_refltable_reduce_cs
#
# $Id$

from xfel.merging.command_line.dev_mpi_cluster_two_merge import scaling_manager_mpi, Script
from xfel.merging.command_line.dev_cxi_merge_refltable import refltable_scaling_manager

#from xfel.merging.command_line.dev_cxi_merge_refltable import merging_reflection_table as mrt
#from xfel.merging.command_line.dev_cxi_merge_refltable import merging_crystal_table as mct
#from xfel.cxi.merging_utils import null_data

class refltable_scaling_manager_mpi(scaling_manager_mpi, refltable_scaling_manager):
  pass

from xfel.merging.command_line.dev_cxi_merge import Script as base_Script

class Script(base_Script):

  @staticmethod
  def mpi_merge_op(data0, data1, datatype):
    data0.n_accepted += data1.n_accepted
    data0.n_file_error += data1.n_file_error
    data0.n_low_corr += data1.n_low_corr
    data0.n_low_signal += data1.n_low_signal
    data0.n_processed += data1.n_processed
    data0.n_wrong_bravais += data1.n_wrong_bravais
    data0.n_wrong_cell += data1.n_wrong_cell

    data0.completeness += data1.completeness
    data0.completeness_predictions += data1.completeness_predictions
    data0.summed_N += data1.summed_N
    data0.summed_weight += data1.summed_weight
    data0.summed_wt_I += data1.summed_wt_I

    print("CORR_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.corr_values.extend(data1.corr_values)
    print("CORR_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    print("DMIN_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.d_min_values.extend(data1.d_min_values)
    print("DMIN_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    print("REJFRAC_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.rejected_fractions.extend(data1.rejected_fractions)
    print("REJFRAC_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    print("WAVELENGTH_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.wavelength.extend(data1.wavelength)
    print("WAVELENGTH_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    print("UCVAL_ADDCELLS START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.uc_values.add_cells(data1.uc_values)
    print("UCVAL_ADDCELLS END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    print("DICT_MERGE START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.failure_modes = {k : data0.failure_modes.get(k, 0) + data1.failure_modes.get(k,0) for k in set(data0.failure_modes.keys()) | set(data1.failure_modes.keys())}
    print("DICT_MERGE END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    print("ISIGI_CID START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    next_crystal_id = len(data0.crystal_table)
    data1.ISIGI['crystal_id'] += next_crystal_id
    print("ISIGI_CID END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    print("ISIGI_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.ISIGI.extend(data1.ISIGI)
    print("ISIGI_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    print("CTABLE_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.crystal_table.extend(data1.crystal_table)
    print("CTABLE_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    print("SEQ_ADD START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
    data0.finished_db_mgr.sequencer += data1.finished_db_mgr.sequencer
    print("SEQ_ADD END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    if not data0.params.short_circuit:
      print("OBS_EXTEND START RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))
      data0.observations.extend(data1.observations)
      print("OBS_EXTEND END RANK=%d:%d TIME=%f;"%(data0.myRank,data1.myRank,tt()))

    return data0

  @staticmethod
  def combine_frames (data0, data1, datatype) :
    """
    Combine the scaled data from a frame with the current overall dataset.
    Also accepts None or null_data objects, when data are unusable but we
    want to record the file as processed.
    """
    data0.n_processed += 1
    if (data1 is None) :
      return data0
    if (isinstance(data1, null_data)) :
      if (data1.file_error) :
        data0.n_file_error += 1
      elif (data1.low_signal) :
        data0.n_low_signal += 1
      elif (data1.wrong_bravais) :
        data0.n_wrong_bravais += 1
      elif (data1.wrong_cell) :
        data0.n_wrong_cell += 1
      elif (getattr(data1,"reason",None) is not None):
        if str(data1.reason)!="":
          data0.failure_modes[str(data1.reason)] = data0.failure_modes.get(str(data1.reason),0) + 1
        elif repr(type(data1.reason))!="":
          data0.failure_modes[repr(type(data1.reason))] = data0.failure_modes.get(repr(type(data1.reason)),0) + 1
        else:
          data0.failure_modes["other reasons"] = data0.failure_modes.get("other reasons",0) + 1
      return data0
    if (data1.accept) :
      data0.n_accepted    += 1
      data0.completeness  += data1.completeness
      data0.completeness_predictions += data1.completeness_predictions
      data0.summed_N      += data1.summed_N
      data0.summed_weight += data1.summed_weight
      data0.summed_wt_I   += data1.summed_wt_I
      data0.ISIGI.extend(data1.ISIGI)
    else :
      data0.n_low_corr += 1
    data0.uc_values.add_cell(data1.indexed_cell,
      rejected=(not data1.accept))
    if not data0.params.short_circuit:
      data0.observations.append(data1.n_obs)
    if (data1.n_obs > 0) :
      frac_rejected = data1.n_rejected / data1.n_obs
      data0.rejected_fractions.append(frac_rejected)
      data0.d_min_values.append(data1.d_min)
    data0.corr_values.append(data1.corr)
    data0.wavelength.append(data1.wavelength)
    data0.finished_db_mgr.sequencer += data1.finished_db_mgr.sequencer
    return data0

  def validate(self):
    base_Script.validate(self)
    if (self.params.rescale_with_average_cell):
      raise Usage("""Rescaling_with_average_cell not supported with MPI
      (Would require a second round of scaling, inefficient).""")

  def run(self,comm,timing=True):
    from mpi4py import MPI
    rank = comm.Get_rank()
    size = comm.Get_size()
    merge_op = MPI.Op.Create(self.mpi_merge_op, commute=True)

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

    tar_file_names = transmitted_info["file_names"]

    if timing: print("SCALER_WORKERS START RANK=%d TIME=%f"%(rank,tt()))
    if rank == 0:
      for ix in range(len(tar_file_names)):
        if timing: print("SCALER_WORKER_RECV START=%d RANK=%d TIME=%f"%(ix,rank,tt()))
        rankreq = comm.recv(source=MPI.ANY_SOURCE)
        if timing: print("SCALER_WORKER_RECV START=%d RANK=%d TIME=%f"%(ix,rank,tt()))
        if timing: print("SCALER_WORKER_SEND START=%d RANK=%d,%d TIME=%f"%(ix,rank,rankreq,tt()))
        comm.send(ix,dest=rankreq)
        if timing: print("SCALER_WORKER_SEND END=%d RANK=%d,%d TIME=%f"%(ix,rank,rankreq,tt()))
      for rankreq in range(size-1):
        if timing: print("SCALER_WORKER_RECV_KILL START RANK=%d TIME=%f"%(rank,tt()))
        rankreq = comm.recv(source=MPI.ANY_SOURCE)
        if timing: print("SCALER_WORKER_RECV_KILL END RANK=%d TIME=%f"%(rank,tt()))
        if timing: print("SCALER_WORKER_SEND_KILL START RANK=%d,%d TIME=%f"%(rank,rankreq,tt()))
        comm.send('endrun',dest=rankreq)
        if timing: print("SCALER_WORKER_SEND_KILL END RANK=%d,%d TIME=%f"%(rank,rankreq,tt()))
      scaler_worker.finished_db_mgr = db_mgr

    else:
      while True:
        if timing: print("SCALER_WORKER_RANKSEND START RANK=%d TIME=%f"%(rank,tt()))
        comm.send(rank,dest=0)
        if timing: print("SCALER_WORKER_RANKSEND END RANK=%d TIME=%f"%(rank,tt()))
        if timing: print("SCALER_WORKER_IDXRECV START RANK=%d TIME=%f"%(rank,tt()))
        idx = comm.recv(source=0)
        if timing: print("SCALER_WORKER_IDXRECV END  RANK=%d TIME=%f"%(rank,tt()))

        if idx == 'endrun':
          scaler_worker.finished_db_mgr = db_mgr
          break
        if timing: print("SCALER_WORKER START=%s RANK=%d TIME=%f"%(str(tar_file_names[idx]), rank, tt()))
        scaler_worker._scale_all_serial([tar_file_names[idx],], db_mgr)
        if timing: print("SCALER_WORKER END=%s RANK=%d TIME=%f"%(str(tar_file_names[idx]), rank, tt()))

    if timing: print("SCALER_WORKERS END RANK=%d TIME=%f"%(rank, tt()))

    # might want to clean up a bit before returning
    del scaler_worker.log
    #del scaler_worker.params
    del scaler_worker.miller_set
    del scaler_worker.i_model
    del scaler_worker.reverse_lookup
    scaler_worker.myRank = rank

    if timing: print("SCALER_WORKERS_REDUCE START RANK=%d TIME=%f"%(rank, tt()))
    scaler_workers = comm.reduce(scaler_worker, op=merge_op, root=0)
    if timing: print("SCALER_WORKERS_REDUCE END RANK=%d TIME=%f"%(rank, tt()))

    MPI.Finalize()
    if rank == 0:
      if timing: print("SCALER_MASTER_ADD START RANK=%d TIME=%f"%(rank, tt()))
      scaler_master._add_all_frames(scaler_workers)
      if timing: print("SCALER_MASTER_ADD END RANK=%d TIME=%f"%(rank, tt()))

      for call_instance in scaler_workers.finished_db_mgr.sequencer:
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

  script = Script(refltable_scaling_manager_mpi)
  result = script.run(comm=comm,timing=False)
  if rank == 0:
    script.show_plot(result)

# LIBTBX_SET_DISPATCHER_NAME mpi.merge_cs
from __future__ import division
import sys,time,os
from xfel.command_line import cxi_merge
from libtbx.utils import Usage, multi_out
from cctbx.array_family import flex
from libtbx import easy_pickle

from xfel.merging.command_line.single_node_merge import get_observations

from xfel.merging.command_line.single_node_merge import scaling_manager as scaling_manager_base
class scaling_manager_mpi(scaling_manager_base):

  def insert_frame(self, **kwargs):
    order = []
    if self.params.postrefinement.enable==True and \
        self.params.postrefinement.algorithm in ["rs2","rs_hybrid"]:
      order_dict = {'wavelength': 1,
                  'beam_x': 2,
                  'beam_y': 3,
                  'distance': 4,
                  'G': 5,
                  'BFACTOR': 6,
                  'RS': 7,
                  'Astar_1': 8,
                  'Astar_2': 9,
                  'Astar_3': 10,
                  'Astar_4': 11,
                  'Astar_5': 12,
                  'Astar_6': 13,
                  'Astar_7': 14,
                  'Astar_8': 15,
                  'Astar_9': 16,
                  'thetax': 17,
                  'thetay': 18,
                  'unique_file_name': 19,
                  'c_c': 20}
    else:
      order_dict = {'wavelength': 1,
                  'beam_x': 2,
                  'beam_y': 3,
                  'distance': 4,
                  'c_c': 5,
                  'slope': 6,
                  'offset': 7,
                  'res_ori_1': 8,
                  'res_ori_2': 9,
                  'res_ori_3': 10,
                  'res_ori_4': 11,
                  'res_ori_5': 12,
                  'res_ori_6': 13,
                  'res_ori_7': 14,
                  'res_ori_8': 15,
                  'res_ori_9': 16,
                  'rotation100_rad': 17,
                  'rotation010_rad': 18,
                  'rotation001_rad': 19,
                  'half_mosaicity_deg': 20,
                  'wave_HE_ang': 21,
                  'wave_LE_ang': 22,
                  'domain_size_ang':23,
                  'unique_file_name': 24}
    for key in kwargs.keys():
      order.append(order_dict[key])
    parameters = [kwargs.values()]

    # Pick up the index of the row just added.  The file name is
    # assumed to to serve as a unique key.
    lastrowid_key = kwargs['unique_file_name']
    self._db_commands_queue.put(('frame', order, parameters, lastrowid_key))
    while True:
      item = self._db_results_queue.get()
      self._db_results_queue.task_done()
      if item[0] == kwargs['unique_file_name']:
        # Entry in the observation table is zero-based.
        return item[1] - 1
      else:
        # If the key does not match, put it back in the queue for
        # someone else to pick up.
        self._db_results_queue.put(item)

  def insert_observation(self, **kwargs):
    order = []
    order_dict = {'hkl_id_0_base': 0,
                  'i': 1,
                  'sigi': 2,
                  'detector_x': 3,
                  'detector_y': 4,
                  'frame_id_0_base': 5,
                  'overload_flag': 6,
                  'original_h': 7,
                  'original_k': 8,
                  'original_l': 9}
    for key in kwargs.keys():
      order.append(order_dict[key])
    parameters = zip(*kwargs.values())
    self._db_commands_queue.put(('observation', order, parameters, None))

  def mpi_initialize (self, file_names) :
    tar_list,self.integration_pickle_names = file_names
    self.t1 = time.time()

    assert self.params.backend == 'FS' # for prototyping rank=0 marshalling
    from xfel.cxi.merging_database_fs import manager
    db_mgr = manager(self.params)
    db_mgr.initialize_db(self.miller_set.indices())
    self.master_db_mgr = db_mgr

  def mpi_finalize (self) :
    t2 = time.time()
    print >> self.log, ""
    print >> self.log, "#" * 80
    print >> self.log, "FINISHED MERGING"
    print >> self.log, "  Elapsed time: %.1fs" % (t2 - self.t1)
    print >> self.log, "  %d of %d integration files were accepted" % (
      self.n_accepted, len(self.integration_pickle_names))
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
    assert checksum == len(self.integration_pickle_names)

    high_res_count = (self.d_min_values <= self.params.d_min).count(True)
    print >> self.log, "Of %d accepted images, %d accepted to %5.2f Angstrom resolution" % \
      (self.n_accepted, high_res_count, self.params.d_min)

    if self.params.raw_data.sdfac_refine:
      self.scale_errors()

    if self.params.raw_data.errors_from_sample_residuals:
      self.errors_from_residuals()

  def iterate_instances(self,item,rank,tt):
    count=0
    n_insert_frame = self.master_db_mgr.insert_frame

    for call_instance in item.finished_db_mgr.sequencer:
      print "INSTANCE=%d START RANK=%d TIME=%f"%(count,rank,tt())
      if call_instance["call"] == "insert_frame":
        print "insert_frame=%d START RANK=%d TIME=%f"%(count,rank,tt())
        frame_id_zero_base = n_insert_frame(**call_instance["data"])
        print "insert_frame=%d STOP RANK=%d TIME=%f"%(count,rank,tt())

      elif call_instance["call"] == "insert_observation":
        print "insert_observation=%d START RANK=%d TIME=%f"%(count,rank,tt())
        call_instance["data"]['frame_id_0_base'] = [frame_id_zero_base] * len(call_instance["data"]['frame_id_0_base'])
        self.master_db_mgr.insert_observation(**call_instance["data"])
        print "insert_observation=%d STOP RANK=%d TIME=%f"%(count,rank,tt())
      count+=1

class run_manager(object):
 def initialize(self,args):
  import iotbx.phil
  from xfel.command_line.cxi_merge import master_phil
  if ("--help" in args) :
    iotbx.phil.parse(master_phil).show(attributes_level=2)
    return
  phil = iotbx.phil.process_command_line(args=args, master_string=master_phil).show()
  work_params = phil.work.extract()
  from xfel.merging.phil_validation import application
  application(work_params)

  if ((work_params.d_min is None) or
      (work_params.data is None) or
      ( (work_params.model is None) and work_params.scaling.algorithm != "mark1") ) :
    command_name = os.environ["LIBTBX_DISPATCHER_NAME"]
    raise Usage(command_name + " "
                "d_min=4.0 "
                "data=~/scratch/r0220/006/strong/ "
                "model=3bz1_3bz2_core.pdb")
  if ((work_params.rescale_with_average_cell) and
      (not work_params.set_average_unit_cell)) :
    raise Usage("If rescale_with_average_cell=True, you must also specify "+
      "set_average_unit_cell=True.")
  if [work_params.raw_data.sdfac_auto, work_params.raw_data.sdfac_refine, work_params.raw_data.errors_from_sample_residuals].count(True) > 1:
    raise Usage("Specify only one of sdfac_auto, sdfac_refine or errors_from_sample_residuals.")

  # Read Nat's reference model from an MTZ file.  XXX The observation
  # type is given as F, not I--should they be squared?  Check with Nat!
  log = open("%s.log" % work_params.output.prefix, "w")
  out = multi_out()
  out.register("log", log, atexit_send_to=None)
  out.register("stdout", sys.stdout)
  print >> out, "I model"
  if work_params.model is not None:
    from xfel.merging.general_fcalc import run as run_fmodel
    i_model = run_fmodel(work_params)
    work_params.target_unit_cell = i_model.unit_cell()
    work_params.target_space_group = i_model.space_group_info()
    i_model.show_summary()
  else:
    i_model = None

  print >> out, "Target unit cell and space group:"
  print >> out, "  ", work_params.target_unit_cell
  print >> out, "  ", work_params.target_space_group
  from xfel.command_line.cxi_merge import consistent_set_and_model
  self.miller_set, self.i_model = consistent_set_and_model(work_params,i_model)
  self.work_params = work_params
  self.frame_files = get_observations(work_params)
  self.out = out

 def other(self,scaler):
  out = self.out
  work_params = self.work_params
  miller_set = self.miller_set
  if scaler.n_accepted == 0:
    return None
  scaler.show_unit_cell_histograms()
  if (work_params.rescale_with_average_cell) :
    average_cell_abc = scaler.uc_values.get_average_cell_dimensions()
    average_cell = uctbx.unit_cell(list(average_cell_abc) +
      list(work_params.target_unit_cell.parameters()[3:]))
    work_params.target_unit_cell = average_cell
    print >> out, ""
    print >> out, "#" * 80
    print >> out, "RESCALING WITH NEW TARGET CELL"
    print >> out, "  average cell: %g %g %g %g %g %g" % \
      work_params.target_unit_cell.parameters()
    print >> out, ""
    assert False,"must do this step again with MPI"
    scaler.reset()
    scaler.scale_all(frame_files)
    scaler.show_unit_cell_histograms()
  print >> out, "\n"

  # Sum the observations of I and I/sig(I) for each reflection.
  sum_I = flex.double(miller_set.size(), 0.)
  sum_I_SIGI = flex.double(miller_set.size(), 0.)
  for i in xrange(miller_set.size()) :
    index = miller_set.indices()[i]
    if index in scaler.ISIGI :
      for t in scaler.ISIGI[index]:
        sum_I[i] += t[0]
        sum_I_SIGI[i] += t[1]

  miller_set_avg = miller_set.customized_copy(
    unit_cell=work_params.target_unit_cell)
  table1 = cxi_merge.show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.completeness,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for all reflections",
    out=out,
    work_params=work_params)
  print >> out, ""
  if work_params.model is not None:
    n_refl, corr = scaler.get_overall_correlation(sum_I)
  else:
    n_refl, corr = ((scaler.completeness > 0).count(True), 0)
  print >> out, "\n"
  table2 = cxi_merge.show_overall_observations(
    obs=miller_set_avg,
    redundancy=scaler.summed_N,
    redundancy_to_edge=scaler.completeness_predictions,
    summed_wt_I=scaler.summed_wt_I,
    summed_weight=scaler.summed_weight,
    ISIGI=scaler.ISIGI,
    n_bins=work_params.output.n_bins,
    title="Statistics for reflections where I > 0",
    out=out,
    work_params=work_params)
  #from libtbx import easy_pickle
  #easy_pickle.dump(file_name="stats.pickle", obj=stats)
  #stats.report(plot=work_params.plot)
  #miller_counts = miller_set_p1.array(data=stats.counts.as_double()).select(
  #  stats.counts != 0)
  #miller_counts.as_mtz_dataset(column_root_label="NOBS").mtz_object().write(
  #  file_name="nobs.mtz")
  if work_params.data_subsubsets.subsubset is not None and work_params.data_subsubsets.subsubset_total is not None:
    easy_pickle.dump("scaler_%d.pickle"%work_params.data_subsubsets.subsubset, scaler)
  explanation = """
Explanation:
Completeness       = # unique Miller indices present in data / # Miller indices theoretical in asymmetric unit
Asu. Multiplicity  = # measurements / # Miller indices theoretical in asymmetric unit
Obs. Multiplicity  = # measurements / # unique Miller indices present in data
Pred. Multiplicity = # predictions on all accepted images / # Miller indices theoretical in asymmetric unit"""
  print >> out, explanation
  mtz_file, miller_array = scaler.finalize_and_save_data()
  #table_pickle_file = "%s_graphs.pkl" % work_params.output.prefix
  #easy_pickle.dump(table_pickle_file, [table1, table2])
  loggraph_file = os.path.abspath("%s_graphs.log" % work_params.output.prefix)
  f = open(loggraph_file, "w")
  f.write(table1.format_loggraph())
  f.write("\n")
  f.write(table2.format_loggraph())
  f.close()
  result = cxi_merge.scaling_result(
    miller_array=miller_array,
    plots=scaler.get_plot_statistics(),
    mtz_file=mtz_file,
    loggraph_file=loggraph_file,
    obs_table=table1,
    all_obs_table=table2,
    n_reflections=n_refl,
    overall_correlation=corr)
  easy_pickle.dump("%s.pkl" % work_params.output.prefix, result)
  return result

if (__name__ == "__main__"):
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()
  size = comm.Get_size()
  from time import time as tt

  # set things up
  if rank == 0:
    print "SETUP START RANK=%d TIME=%f"%(rank,tt())
    run = run_manager()
    run.initialize(args=sys.argv[1:])
    scaler_master = scaling_manager_mpi(
      miller_set=run.miller_set,
      i_model=run.i_model,
      params=run.work_params,
      log=run.out)
    scaler_master.mpi_initialize(run.frame_files)

    transmitted_info = dict(file_names=run.frame_files,
                            miller_set=run.miller_set,
                            model = run.i_model,
                            params = run.work_params )
    print "SETUP END RANK=%d TIME=%f"%(rank,tt())
  else:
    print "SETUP START RANK=%d TIME=%f"%(rank,tt())
    transmitted_info = None
    print "SETUP END RANK=%d TIME=%f"%(rank,tt())

  print "BROADCAST START RANK=%d TIME=%f"%(rank,tt())
  transmitted_info = comm.bcast(transmitted_info, root = 0)
  print "BROADCAST END RANK=%d TIME=%f"%(rank,tt())

  # now actually do the work
  print "SCALER_WORKER_SETUP START RANK=%d TIME=%f"%(rank,tt())
  scaler_worker = scaling_manager_mpi(transmitted_info["miller_set"],
                                      transmitted_info["model"],
                                      transmitted_info["params"],
                                      log = sys.stdout)
  print "SCALER_WORKER_SETUP END RANK=%d TIME=%f"%(rank,tt())

  assert scaler_worker.params.backend == 'FS' # only option that makes sense
  from xfel.cxi.merging_database_fs import manager2 as manager
  db_mgr = manager(scaler_worker.params)
  tar_file_names = transmitted_info["file_names"][0]

  print "SCALER_WORKERS START RANK=%d TIME=%f"%(rank,tt())
  if rank == 0:
    for ix in xrange(len(tar_file_names)):
      print "SCALER_WORKER_RECV START=%d RANK=%d TIME=%f"%(ix,rank,tt())
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      print "SCALER_WORKER_RECV START=%d RANK=%d TIME=%f"%(ix,rank,tt())
      print "SCALER_WORKER_SEND START=%d RANK=%d,%d TIME=%f"%(ix,rank,rankreq,tt())
      comm.send(ix,dest=rankreq)
      print "SCALER_WORKER_SEND END=%d RANK=%d,%d TIME=%f"%(ix,rank,rankreq,tt())
    for rankreq in range(size-1):
      print "SCALER_WORKER_RECV_KILL START RANK=%d TIME=%f"%(rank,tt())
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      print "SCALER_WORKER_RECV_KILL END RANK=%d TIME=%f"%(rank,tt())
      print "SCALER_WORKER_SEND_KILL START RANK=%d,%d TIME=%f"%(rank,rankreq,tt())
      comm.send('endrun',dest=rankreq)
      print "SCALER_WORKER_SEND_KILL END RANK=%d,%d TIME=%f"%(rank,rankreq,tt())
    scaler_worker.finished_db_mgr = db_mgr

  else:
    while True:
      print "SCALER_WORKER_RANKSEND START RANK=%d TIME=%f"%(rank,tt())
      comm.send(rank,dest=0)
      print "SCALER_WORKER_RANKSEND END RANK=%d TIME=%f"%(rank,tt())
      print "SCALER_WORKER_IDXRECV START RANK=%d TIME=%f"%(rank,tt())
      idx = comm.recv(source=0)
      print "SCALER_WORKER_IDXRECV END  RANK=%d TIME=%f"%(rank,tt())

      if idx == 'endrun':
        scaler_worker.finished_db_mgr = db_mgr
        break
      print "SCALER_WORKER START=%s RANK=%d TIME=%f"%(str([tar_file_names[idx],]), rank, tt())
      scaler_worker.tar_to_scale_frame_adapter(tar_list=[tar_file_names[idx],], db_mgr=db_mgr)
      print "SCALER_WORKER END=%s RANK=%d TIME=%f"%(str([tar_file_names[idx],]), rank, tt())

  print "SCALER_WORKERS END RANK=%d TIME=%f"%(rank, tt())

  # might want to clean up a bit before returning
  del scaler_worker.log
  del scaler_worker.params
  del scaler_worker.miller_set
  del scaler_worker.i_model
  del scaler_worker.reverse_lookup

  # gather reports and all add together
  print "GATHER START RANK=%d TIME=%f"%(rank, tt())
  reports = comm.gather(scaler_worker,root=0)
  print "GATHER END RANK=%d TIME=%f"%(rank, tt())

  if rank == 0:
    print "Processing reports from %d ranks"%(len(reports))
    ireport = 0
    for item in reports:
      print "SCALER_MASTER_ADD START RANK=%d TIME=%f"%(rank, tt())
      scaler_master._add_all_frames(item)
      print "SCALER_MASTER_ADD END RANK=%d TIME=%f"%(rank, tt())

      print "processing %d calls from report %d"%(len(item.finished_db_mgr.sequencer),ireport); ireport += 1
      scaler_master.iterate_instances(item,rank,tt)

    print "SCALER_MASTER_FINALISE START RANK=%d TIME=%f"%(rank, tt())
    scaler_master.master_db_mgr.join() # database written, finalize the manager
    scaler_master.mpi_finalize()
    print "SCALER_MASTER_FINALISE END RANK=%d TIME=%f"%(rank, tt())
    print "OK"
    run.other(scaler_master)

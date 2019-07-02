# LIBTBX_SET_DISPATCHER_NAME prime.mpi_run
"""
Find initial scaling factors for all integration results
"""
from __future__ import absolute_import, division, print_function
from mpi4py import MPI
import sys, os
from prime.postrefine.mod_input import process_input, read_pickles
from prime.postrefine.mod_util import intensities_scaler
from prime.postrefine.mod_merge_data import merge_data_handler
from prime.postrefine import postref_handler
from cctbx.array_family import flex
import time, math
from six.moves import range
from six.moves import zip

# setup mpi
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
assert size>1

def master(frame_token, iparams, activity):
  if activity == "scale":
    n_batch = 1
    frame_objects = frame_token
    indices = range(0, len(frame_objects), n_batch)
    for i in indices:
      i_end = i+n_batch if i+n_batch < len(frame_objects) else len(frame_objects)
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send((activity, (frame_objects[i:i_end], iparams)), dest=rankreq)
  elif activity == "pre_merge":
    frame_objects, avg_mode = frame_token
    n_batch = int(len(frame_objects)/(size*3))
    if n_batch < 10: n_batch = 10
    indices = range(0, len(frame_objects), n_batch)
    for i in indices:
      i_end = i+n_batch if i+n_batch < len(frame_objects) else len(frame_objects)
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send((activity, (frame_objects[i:i_end], iparams, avg_mode)), dest=rankreq)
  elif activity == "merge":
    frame_objects, avg_mode = frame_token
    its = intensities_scaler()
    cpo = its.combine_pre_merge(frame_objects, iparams)
    #assign at least 100k reflections at a time
    n_batch = int(1e5/(len(cpo[1])/cpo[0]))
    if n_batch < 1: n_batch = 1
    print("Merging with %d batch size"%(n_batch))
    indices = range(0, cpo[0], n_batch)
    for i in indices:
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      i_end = i+n_batch if i+n_batch < cpo[0] else cpo[0]
      sel = flex.bool([sel_l and sel_h for sel_l, sel_h in zip(cpo[1]>=i, cpo[1]<i_end)])
      batch_prep = [cpo_elem.select(sel) for cpo_elem in cpo[1:13]]
      batch_prep.insert(0, i_end-i)
      batch_prep[1] -= i
      batch_prep.append(cpo[13])
      batch_prep.append(cpo[14])
      batch_prep.append(cpo[15].select(sel))
      batch_prep.append("")
      comm.send((activity, (tuple(batch_prep), iparams, avg_mode)), dest=rankreq)
  elif activity == "postref":
    frame_files, miller_array_ref, results, avg_mode = frame_token
    #convert results to a dict obj so that pickle_filename is the key
    pres_dict = {}
    for pres in results:
      if pres: pres_dict[pres.pickle_filename]=pres
    for i in range(len(frame_files)):
      rankreq = comm.recv(source=MPI.ANY_SOURCE)
      comm.send((activity, (i, frame_files[i], iparams, miller_array_ref, pres_dict[frame_files[i]] if frame_files[i] in pres_dict else None, avg_mode)), dest=rankreq)
  print("Master for %s is completed. Time to stop all %d clients"%(activity, size-1))
  # stop clients
  for rankreq in range(size-1):
    rankreq = comm.recv(source=MPI.ANY_SOURCE)
    comm.send('endrun', dest=rankreq)

def client():
  result = []
  prh = postref_handler()
  while True:
    comm.send(rank, dest=0)
    msg = comm.recv(source=0)
    if str(msg) == 'endrun': break
    #receive contents for processing
    activity, act_params = msg
    if activity == "scale":
      frame_files, iparams = act_params
      for frame_index, frame_file in enumerate(frame_files):
        pres, _ = prh.scale_frame_by_mean_I(frame_index, frame_file, iparams, 0, 'average')
        result.append(pres)
    elif activity == "pre_merge":
      frame_results, iparams, avg_mode = act_params
      its = intensities_scaler()
      prep_output = its.prepare_output(frame_results, iparams, avg_mode)
      result.append(prep_output)
    elif activity == "merge":
      batch_prep, iparams, avg_mode = act_params
      its = intensities_scaler()
      mdh, _, txt_out_rejection  = its.calc_avg_I_cpp(batch_prep, iparams, avg_mode)
      result.append([mdh, txt_out_rejection])
    elif activity == "postref":
      frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode = act_params
      pres, _ = prh.postrefine_by_frame(frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode)
      result.append(pres)
  return result

def merge(pres_result, iparams, txt_out_prefix, output_prefix, avg_mode):
  #pre-merge task
  if rank == 0:
    pres_results = sum(pres_result, [])
    print("Scaling/post-refinement is done on %d cores for %d frames"%(size, len(pres_results)))
    master((pres_results, avg_mode), iparams, "pre_merge")
    premerge_result = []
  else:
    pres_results = None
    premerge_result = client()
  premerge_result = comm.gather(premerge_result, root=0)
  comm.Barrier()
  #merge task
  if rank == 0:
    print("Pre-merge is done on %d cores"%(len(premerge_result)))
    master((premerge_result, avg_mode), iparams, "merge")
    merge_result = []
  else:
    merge_result = client()
  #finalize merge
  merge_result = comm.gather(merge_result, root=0)
  comm.Barrier()
  if rank == 0:
    print("Merge completed on %d cores"%(len(merge_result)))
    merge_results = sum(merge_result, [])
    mdh = merge_data_handler()
    txt_out_rejection = ""
    for _mdh, _txt_out_rejection in merge_results:
      mdh.extend(_mdh)
      txt_out_rejection += _txt_out_rejection
    #selet only indices with non-Inf non-Nan stats
    selections = flex.bool([False if (math.isnan(r0) or math.isinf(r0) or math.isnan(r1) or math.isinf(r1)) else True for r0, r1  in zip(mdh.r_meas_div, mdh.r_meas_divisor)])
    mdh.reduce_by_selection(selections)
    its = intensities_scaler()
    mdh, txt_merge_mean_table = its.write_output(mdh, iparams, os.path.join(iparams.run_no, output_prefix), avg_mode)
    print(txt_merge_mean_table)
    #write log output
    with open(os.path.join(iparams.run_no, 'log.txt'), 'a') as f:
      f.write(txt_out_prefix+txt_merge_mean_table)
    with open(os.path.join(iparams.run_no, 'rejections.txt'), 'a') as f:
      f.write(txt_out_rejection)
  else:
    merge_results = None
    mdh = None
  return pres_results, mdh

def run(argv):
  comm.Barrier()
  start_time = MPI.Wtime()
  #broadcast parameters
  if rank == 0:
    iparams, txt_out_input = process_input(argv)
    iparams.flag_volume_correction = False
    iparams.flag_hush = True
    print(txt_out_input)
    frame_files = read_pickles(iparams.data)
  else:
    iparams = None
    frame_files = None
    txt_out_input = None
  comm.Barrier()
  #assign scaling task
  if rank == 0:
    master(frame_files, iparams, "scale")
    result = []
  else:
    result = client()
  result = comm.gather(result, root=0)
  comm.Barrier()
  _, mdh = merge(result, iparams, txt_out_input, "mean_scaled", "average")
  #postrefine task
  if rank == 0:
    n_postref_cycle = iparams.n_postref_cycle
    results = [None]*len(frame_files)
  else:
    n_postref_cycle = None
    results = None
  n_postref_cycle = comm.bcast(n_postref_cycle, root=0)
  avg_mode = 'weighted'
  for i_iter in range(n_postref_cycle):
    comm.Barrier()
    if i_iter == n_postref_cycle -1: avg_mode = 'final'
    if rank == 0:
      print("Start post-refinement cycle %d"%(i_iter+1))
      master((frame_files, mdh.miller_array_merge, results, avg_mode), iparams, "postref")
      result = []
    else:
      result = client()
    result = comm.gather(result, root=0)
    comm.Barrier()
    results, mdh = merge(result, iparams, "", "postref_cycle_%d"%(i_iter+1), avg_mode)
  #collect time profile
  comm.Barrier()
  end_time = MPI.Wtime()
  txt_time = 'Elapsed Time (s):%10.2f\n'%(end_time-start_time)
  MPI.Finalize()

if __name__ == "__main__":
  argv = sys.argv[1:] if len(sys.argv) > 1 else None
  run(argv)

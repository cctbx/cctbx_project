# LIBTBX_SET_DISPATCHER_NAME prime.cluster_isoform
""" cluster diffraction images by """
from __future__ import absolute_import, division, print_function
import six
__author__ = 'Monarin Uervirojnangkoorn, monarin@gmail.com'

from prime.isoform_cluster.mod_isoform_cluster import isoform_cluster_handler
from prime.index_ambiguity.mod_kmeans import kmeans_handler
from prime.postrefine.mod_mx import mx_handler
from six.moves import cPickle as pickle
from libtbx.easy_mp import pool_map
import numpy as np
import random, os, sys

def solve_with_mtz_mproc(args):
  frame_no, pickle_filename, main_obs, iparams, miller_array_ref_set = args
  isoch = isoform_cluster_handler()
  cluster_id, txt_out = isoch.calc_cc(frame_no, pickle_filename, main_obs, iparams, miller_array_ref_set)
  print(txt_out)
  return pickle_filename, cluster_id

def calculate_r_mproc(args):
  frame_no, pickle_filename, obs_main, obs_list = args
  isoch = isoform_cluster_handler()
  r_set, main_obs, txt_out = isoch.calc_r(frame_no, pickle_filename, obs_main, obs_list)
  print(txt_out)
  return pickle_filename, r_set, main_obs

def get_obs_mproc(args):
  frame_no, pickle_filename, iparams = args
  isoch = isoform_cluster_handler()
  obs_asu = isoch.get_observations(pickle_filename, iparams)
  return obs_asu, pickle_filename

def write_out_solutions(iparams, sol_pickle):
  cluster_files = [os.path.join(iparams.run_no,'cluster_'+str(i)+'.lst') for i in range(iparams.isoform_cluster.n_clusters)]
  cluster_text_set = ['']*iparams.isoform_cluster.n_clusters
  for frame, cluster_id in six.iteritems(sol_pickle):
    cluster_text_set[cluster_id] += frame + '\n'
  for i in range(iparams.isoform_cluster.n_clusters):
    with open(cluster_files[i], 'w') as f: f.write(cluster_text_set[i])

class isoform_cluster_run_handler(object):
  def __init__(self):
    self.module_name = "isoform_cluster"

  def get_observation_set(self, iparams, frame_files, n_sel_frames):
    n_frames = len(frame_files)
    frames = [(i, frame_files[i], iparams) for i in random.sample(range(n_frames), n_sel_frames)]
    #get observations list
    print("Reading observations")
    obs_results = pool_map(
          iterable=frames,
          func=get_obs_mproc,
          processes=iparams.n_processors)
    frame_files_sel = []
    obs_list = []
    for result in obs_results:
      if result:
        main_asu, pickle_filename = result
        if main_asu:
          frame_files_sel.append(pickle_filename)
          obs_list.append(main_asu)
    return frame_files_sel, obs_list

  def run(self, args):
    #read inputs
    from prime.postrefine.mod_input import process_input, read_pickles
    iparams, txt_out_input = process_input(args)
    print(txt_out_input)
    with open(os.path.join(iparams.run_no,self.module_name,'log.txt'), 'w') as f: f.write(txt_out_input)
    #read all integration pickles
    frame_files = read_pickles(iparams.data)
    n_frames = len(frame_files)
    if n_frames == 0:
      print("No integration pickle found. Exit program.")
      return None, iparams
    #start
    if iparams.isoform_cluster.isorefin:
      #get collection of iso. ref. reflection set.
      mxh = mx_handler()
      miller_array_ref_set = []
      for isorefin in iparams.isoform_cluster.isorefin:
        flag_ref_found, miller_array_ref = mxh.get_miller_array_from_reflection_file(isorefin)
        if flag_ref_found: miller_array_ref_set.append(miller_array_ref)
      #get observation list
      frame_files_sel, obs_list = self.get_observation_set(iparams, frame_files, n_frames)
      if miller_array_ref_set:
        frames = [(i, frame_files_sel[i], obs_list[i], iparams, miller_array_ref_set) for i in range(len(obs_list))]
        cc_results = pool_map(
          iterable=frames,
          func=solve_with_mtz_mproc,
          processes=iparams.n_processors)
        sol_pickle = {}
        for result in cc_results:
          pickle_filename, cluster_id = result
          sol_pickle[pickle_filename] = cluster_id
        write_out_solutions(iparams, sol_pickle)
        txt_out = "Cluster images with given "+str(len(miller_array_ref_set))+" mtz files completed. Use cluster_0.lst - cluster_k.lst (for k clusters) for merging.\n"
        print(txt_out)
        with open(os.path.join(iparams.run_no,self.module_name,'log.txt'), 'a') as f: f.write(txt_out)
      return

    #*************************************************
    #solve with Brehm & Diederichs - sample size n_sample_frames then bootstrap the rest
    txt_out = "Cluster images with B&D algorithms.\n"
    frame_files_sel, obs_list = self.get_observation_set(iparams, frame_files, iparams.isoform_cluster.n_sample_frames)
    frames = [(i, frame_files_sel[i], obs_list[i], obs_list) for i in range(len(frame_files_sel))]
    #calculate r
    print("Calculating R")
    calc_r_results = pool_map(
          iterable=frames,
          func=calculate_r_mproc,
          processes=iparams.n_processors)
    frame_files_sel = []
    r_matrix = []
    obs_list = []
    for result in calc_r_results:
      if result:
        pickle_filename, r_set, obs = result
        frame_files_sel.append(pickle_filename)
        obs_list.append(obs)
        if len(r_matrix) == 0:
          r_matrix = r_set
        else:
          r_matrix = np.append(r_matrix, r_set, axis=0)
    #choose groups with best R
    print("Selecting frames with best R")
    i_mean_r = np.argsort(np.mean(r_matrix, axis=1))[::-1]
    r_matrix_sorted = r_matrix[i_mean_r]
    frame_files_sorted = np.array(frame_files_sel)[i_mean_r]
    obs_list_sorted = np.array(obs_list)[i_mean_r]
    frame_files_sel = []
    obs_sel = []
    for frame_file, r_set, obs in zip(frame_files_sorted, r_matrix_sorted, obs_list_sorted):
      if frame_file not in frame_files_sel:
        frame_files_sel.append(frame_file)
        obs_sel.append(obs)
        print(frame_file, np.mean(r_set))
        if len(frame_files_sel) >= iparams.isoform_cluster.n_selected_frames:
          print('Found all %6.0f good frames'%(len(frame_files_sel)))
          break
    #Recalculate r for the new selected list
    frames = [(i, frame_files_sel[i], obs_sel[i], obs_sel) for i in range(len(frame_files_sel))]
    print("Re-calculating R")
    calc_r_results = pool_map(
          iterable=frames,
          func=calculate_r_mproc,
          processes=iparams.n_processors)
    frame_files_sel = []
    r_matrix = []
    obs_list = []
    for result in calc_r_results:
      if result:
        pickle_filename, r_set, obs = result
        frame_files_sel.append(pickle_filename)
        obs_list.append(obs)
        if len(r_matrix) == 0:
          r_matrix = r_set
        else:
          r_matrix = np.append(r_matrix, r_set, axis=0)
    print("Minimizing frame distance")
    isoch = isoform_cluster_handler()
    x_set = isoch.optimize(r_matrix, flag_plot=iparams.flag_plot)
    print("Clustering results")
    kmh = kmeans_handler()
    k = iparams.isoform_cluster.n_clusters
    centroids, labels = kmh.run(x_set, k, flag_plot=iparams.flag_plot)
    print("Get solution pickle and cluster files list")
    sol_pickle, cluster_files = isoch.assign_cluster(frame_files_sel, labels, k, \
        os.path.join(iparams.run_no,self.module_name))
    #if more frames found, merge the sample frames to get a reference set
    #that can be used for breaking the ambiguity.
    if n_frames > iparams.isoform_cluster.n_selected_frames:
      print("Assign cluster_id for the remaining images.")
      old_iparams_data = iparams.data[:]
      miller_array_ref_set = []
      from prime.command_line.postrefine import scale_frames, merge_frames
      for i in range(k):
        #generate a reference set from solved frames
        with open(cluster_files[i]) as f: frame_files_processed = f.read().split('\n')[:-1]
        scaled_pres_set = scale_frames(range(len(frame_files_processed)), frame_files_processed, iparams)
        mdh, txt_merge_out = merge_frames(scaled_pres_set, iparams, \
            mtz_out_prefix=os.path.join(self.module_name,'cluster_'+str(i)))
        miller_array_ref_set.append(mdh.miller_array_merge)
        txt_out += txt_merge_out
      #setup a list of remaining frames
      frame_files_remain = [frame for frame in frame_files if frame not in sol_pickle]
      frame_files_remain_sel, obs_remain_sel_list = self.get_observation_set(iparams, \
          frame_files_remain, len(frame_files_remain))
      frames = [(i, frame_files_remain_sel[i], obs_remain_sel_list[i], iparams, miller_array_ref_set) for i in range(len(obs_remain_sel_list))]
      cc_results = pool_map(
          iterable=frames,
          func=solve_with_mtz_mproc,
          processes=iparams.n_processors)
      for result in cc_results:
        if result:
          pickle_filename, cluster_id = result
          sol_pickle[pickle_filename] = cluster_id
      iparams.data = old_iparams_data[:]
    #write out solution pickle
    write_out_solutions(iparams, sol_pickle)
    #write out text output
    txt = "Cluster images completed. Use cluster_0.lst - cluster_k.lst (for k clusters) for merging.\n"
    txt_out += txt
    print(txt)
    with open(os.path.join(iparams.run_no,self.module_name,'log.txt'), 'a') as f: f.write(txt_out)

if (__name__ == "__main__"):
  isoc_runh = isoform_cluster_run_handler()
  isoc_runh.run(sys.argv[1:] if len(sys.argv) > 1 else None)

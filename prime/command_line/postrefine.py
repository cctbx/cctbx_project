from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.run
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Description : Main commandline for prime.
'''

from libtbx.easy_mp import pool_map
from cctbx.array_family import flex
from prime.command_line.solve_indexing_ambiguity import indexing_ambiguity_handler
from prime.postrefine import prepare_output, calc_avg_I_cpp, write_output
from prime.postrefine.mod_mx import mx_handler
import os, sys, math, glob
import numpy as np
from datetime import datetime, time

def determine_mean_I_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_file, iparams, avg_mode = args
  mean_I = prh.calc_mean_intensity(frame_file, iparams, avg_mode)
  return mean_I

def scale_frame_by_mean_I_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_no, frame_file, iparams, mean_of_mean_I, avg_mode = args
  pres = prh.scale_frame_by_mean_I(frame_no, frame_file, iparams, mean_of_mean_I, avg_mode)
  return pres

def postrefine_by_frame_mproc(args):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode = args
  pres = prh.postrefine_by_frame(frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode)
  return pres

def read_pickles(data):
  frame_files = []
  for p in data:
    if os.path.isdir(p) == False:
      if os.path.isfile(p):
        #check if list-of-pickle text file is given
        pickle_list_file = open(p,'r')
        pickle_list = pickle_list_file.read().split("\n")
      else:
        # p is a glob
        pickle_list = glob.glob(p)
      for pickle_filename in pickle_list:
        if os.path.isfile(pickle_filename):
          frame_files.append(pickle_filename)
    else:
      for pickle_filename in os.listdir(p):
        if pickle_filename.endswith('.pickle'):
          frame_files.append(p+'/'+pickle_filename)
  #check if pickle_dir is given in input file instead of from cmd arguments.
  if len(frame_files)==0:
    print 'No pickle files found.'
    exit()
  return frame_files

def scale_frames(frames, frame_files, iparams):
  """scale frames"""
  avg_mode = 'average'
  if iparams.flag_apply_b_by_frame:
    mean_of_mean_I = 0
  else:
    #Calculate <I> for each frame
    frame_args = [(frame_file, iparams, avg_mode) for frame_file in frame_files]
    determine_mean_I_result = pool_map(
            iterable=frame_args,
            func=determine_mean_I_mproc,
            processes=iparams.n_processors)
    frames_mean_I = flex.double()
    for result in determine_mean_I_result:
      if result is not None:
        mean_I, txt_out_result = result
        if mean_I is not None:
          frames_mean_I.append(mean_I)
    mean_of_mean_I = np.median(frames_mean_I)
  #use the calculate <mean_I> to scale each frame
  frame_args = [(frame_no, frame_file, iparams, mean_of_mean_I, avg_mode) for frame_no, frame_file in zip(frames, frame_files)]
  scale_frame_by_mean_I_result = pool_map(
          iterable=frame_args,
          func=scale_frame_by_mean_I_mproc,
          processes=iparams.n_processors)
  observations_merge_mean_set = []
  for result in scale_frame_by_mean_I_result:
    if result is not None:
      pres, txt_out_result = result
      if pres is not None:
        observations_merge_mean_set.append(pres)
  return observations_merge_mean_set

def merge_frames(pres_set, iparams, avg_mode='average', mtz_out_prefix='mean_scaled'):
  """merge frames using average as the default"""
  miller_array_ref, txt_out = (None,'')
  if pres_set:
    prep_output = prepare_output(pres_set, iparams, avg_mode)
    if prep_output is not None:
      cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort,  I_obs_all_sort, \
      sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort, rs_all_sort, wavelength_all_sort, \
      sin_sq_all_sort, SE_all_sort, uc_mean, wavelength_mean, pickle_filename_all_sort, txt_prep_out = prep_output
      calc_average_I_result = calc_avg_I_cpp(cn_group, group_id_list, miller_indices_all_sort,
          miller_indices_ori_all_sort, I_obs_all_sort,
          sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort,
          rs_all_sort, wavelength_all_sort, sin_sq_all_sort,
          SE_all_sort, avg_mode, iparams, pickle_filename_all_sort)
      miller_index, I_avg, sigI_avg, stat, I_avg_two_halves_tuple, txt_obs_out, txt_reject_out = calc_average_I_result
      I_avg_even, I_avg_odd, I_avg_even_h, I_avg_odd_h, I_avg_even_k, I_avg_odd_k, I_avg_even_l, I_avg_odd_l = I_avg_two_halves_tuple
      txt_out_rejection = txt_reject_out
      #selet only indices with non-Inf non-Nan stats
      s0, s1, s2, s3, s4 = stat
      selections = flex.bool([False if (math.isnan(_s0) or math.isinf(_s0) or math.isnan(_s1) or math.isinf(_s1)) else True for _s0, _s1  in zip(s0, s1)])
      stat_all = [[_s0,_s1,_s2,_s3,_s4] for _s0,_s1,_s2,_s3,_s4 in zip(s0.select(selections),s1.select(selections),
          s2.select(selections), s3.select(selections), s4.select(selections))]
      with open(iparams.run_no+'/rejections.txt', 'a') as f:
        f.write(txt_out_rejection)
      #merge all good indices
      miller_array_ref, txt_merge_mean_table = write_output(miller_index.select(selections),
          I_avg.select(selections), sigI_avg.select(selections),
          stat_all, (I_avg_even.select(selections), I_avg_odd.select(selections),
          I_avg_even_h.select(selections), I_avg_odd_h.select(selections),
          I_avg_even_k.select(selections), I_avg_odd_k.select(selections),
          I_avg_even_l.select(selections), I_avg_odd_l.select(selections)),
          iparams, uc_mean, wavelength_mean, iparams.run_no+'/'+mtz_out_prefix, avg_mode)
      print txt_merge_mean_table
      print txt_prep_out
      txt_out = txt_merge_mean_table + txt_prep_out
  return miller_array_ref, txt_out

def postrefine_frames(i_iter, frames, frame_files, iparams, pres_set, miller_array_ref, avg_mode):
  """postrefine given frames and previous postrefinement results"""
  miller_array_ref = miller_array_ref.generate_bijvoet_mates()
  txt_merge_postref = 'Post-refinement cycle '+str(i_iter+1)+' ('+avg_mode+')\n'
  txt_merge_postref += ' * R and CC show percent change.\n'
  print txt_merge_postref
  frame_args = [(frame_no, frame_file, iparams, miller_array_ref, pres_in, avg_mode) for frame_no, frame_file, pres_in in zip(frames, frame_files, pres_set)]
  postrefine_by_frame_result = pool_map(
      iterable=frame_args,
      func=postrefine_by_frame_mproc,
      processes=iparams.n_processors)
  postrefine_by_frame_good = []
  postrefine_by_frame_pres_list = []
  for results in postrefine_by_frame_result:
    if results is not None:
      pres, txt_out_result = results
      postrefine_by_frame_pres_list.append(pres)
      if pres is not None:
        postrefine_by_frame_good.append(pres)
    else:
      postrefine_by_frame_pres_list.append(None)
  return postrefine_by_frame_good, postrefine_by_frame_pres_list, txt_merge_postref

def run(argv):
  #capture starting time
  time_global_start=datetime.now()
  import logging
  logging.captureWarnings(True)
  formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
  console_handler = logging.StreamHandler()
  console_handler.setLevel(logging.ERROR)
  console_handler.setFormatter(formatter)
  logging.getLogger().addHandler(console_handler)
  logging.getLogger('py.warnings').addHandler(console_handler)
  logging.basicConfig(format='%(asctime)s\t%(levelname)s\t%(message)s', level=logging.DEBUG)
  #0.1 determine indexing ambiguity and setup iparams
  txt_indexing_ambiguity = "Determine if there is an indexing ambiguity on the dataset"
  print txt_indexing_ambiguity
  idah = indexing_ambiguity_handler()
  sol_fname, iparams = idah.run(argv)
  if sol_fname is None:
    print "No ambiguity."
    txt_indexing_ambiguity += "\nNo ambiguity."
  else:
    print "Ambiguity is solved. Solution file was saved to :"+str(sol_fname)
    txt_indexing_ambiguity += "Ambiguity is solved. Solution file was saved to :"+str(sol_fname)
    iparams.indexing_ambiguity.index_basis_in = sol_fname
  #0.2 setup parameters
  iparams.flag_volume_correction = False
  if iparams.partiality_model == "Lognormal":
    iparams.voigt_nu = 0.008 #use voigt_nu as lognpdf zero parameter
  #0.3 read frames
  frame_files = read_pickles(iparams.data)
  frames = range(len(frame_files))
  #1. prepare reference miller array
  txt_merge_mean = 'Generating a reference set (will not be used if hklrefin is set)'
  print txt_merge_mean
  #Always generate the mean-intensity scaled set.
  scaled_pres_set = scale_frames(frames, frame_files, iparams)
  miller_array_ref, _txt_merge_mean = merge_frames(scaled_pres_set, iparams)
  txt_merge_mean += '\n' + _txt_merge_mean
  if not iparams.n_postref_cycle:
    with open(iparams.run_no+'/log.txt', 'a') as f:
      f.write(txt_indexing_ambiguity + txt_merge_mean)
    raise Usage("No. of post-refinement cycle was set to 0. Exit without post-refinement.")
  if iparams.hklrefin is not None:
    mxh = mx_handler()
    _, miller_array_ref = mxh.get_miller_array_from_reflection_file(iparams.hklrefin)
  if miller_array_ref is None:
    raise Usage("Problem with the assigned reference set. Try setting hklrefin=None and rerun the program.")
  #2. Post-refinement
  txt_merge_postref = ''
  postref_pres_set = [None]*len(frames)
  avg_mode = 'weighted'
  for i_iter in xrange(iparams.n_postref_cycle):
    if i_iter == (iparams.n_postref_cycle-1): avg_mode = 'final'
    postref_good_pres_set, postref_pres_set, _txt_merge_postref = postrefine_frames(i_iter, frames, frame_files, iparams, postref_pres_set, miller_array_ref, avg_mode)
    if postref_good_pres_set:
      miller_array_ref, _txt_merge_postref = merge_frames(postref_good_pres_set, iparams,
          avg_mode=avg_mode, mtz_out_prefix='postref_cycle_'+str(i_iter+1))
      txt_merge_postref += _txt_merge_postref
    else:
      raise Usage("Problem with post-refinement. No images refined. Please check your input file.")
  #3. collect caculating time
  time_global_end=datetime.now()
  time_global_spent=time_global_end-time_global_start
  txt_out_time_spent = 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+ \
      ' seconds\nFinished: '+time_global_end.strftime("%A %d. %B %Y %H:%M:%S")+'\n'
  print txt_out_time_spent
  txt_out = txt_indexing_ambiguity + txt_merge_mean + txt_merge_postref + txt_out_time_spent
  with open(iparams.run_no+'/log.txt', 'a') as f:
    f.write(txt_out)
  return txt_out

if (__name__ == "__main__"):
  run(sys.argv[1:] if len(sys.argv) > 1 else None)

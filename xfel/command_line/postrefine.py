# LIBTBX_SET_DISPATCHER_NAME cxi.postrefine
from __future__ import division
from libtbx.utils import multi_out
import os, sys
from libtbx.easy_mp import pool_map
import numpy as np
from cctbx.array_family import flex
from datetime import date, datetime, time, timedelta

def determine_mean_I_mproc(frame_no, frame_files, iph):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  mean_I = prh.calc_mean_intensity(frame_files[frame_no], iph)
  return frame_no, mean_I

def scale_frame_by_mean_I_mproc(frame_no, frame_files, iph, mean_of_mean_I):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  pickle_filename, observations_sel, observations_weight = prh.scale_frame_by_mean_I(
        frame_files[frame_no], iph, mean_of_mean_I)
  return frame_no, observations_sel, observations_weight

def postrefine_by_frame_mproc(frame_no, frame_files, iph, refine_mode, miller_array_ref, G_set, B_factor_set, rotx_set, roty_set, ry_set, rz_set, refine_mode_codes=None):
  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()

  if refine_mode == 'scale_factor_final':
    refine_mode = 'scale_factor_final_'+str(refine_mode_codes[frame_no])
  pickle_filename, observations_refined, observations_weight, cc_r_results, lstsqr_results = prh.postrefine_by_frame(
        refine_mode, frame_files[frame_no], iph, miller_array_ref, scale_factors=(G_set[frame_no], B_factor_set[frame_no]),
        rotations=(rotx_set[frame_no], roty_set[frame_no]), reflecting_ranges=(ry_set[frame_no], rz_set[frame_no]))

  return frame_no, observations_refined, observations_weight, cc_r_results, lstsqr_results

def read_input(args):
  file_name_input = ''
  for i in range(len(args)):
    pair=args[i].split('=')
    if pair[0]=='input':
      file_name_input = pair[1]

  if file_name_input == '':
    print "Please provide input-parameters file (usage: input=yourinput.inp)"
    exit()

  from xfel.cxi.postrefine import postref_handler
  prh = postref_handler()
  iph = prh.read_input_parameters(file_name_input)

  #make run_no folder
  if not os.path.exists(iph.run_no):
    os.makedirs(iph.run_no)

  from xfel.cxi.postrefine import get_observations
  frame_files = get_observations(iph.pickle_dir, 0)

  return iph, frame_files


if (__name__ == "__main__"):
  #capture starting time
  time_global_start=datetime.now()

  #0 .read input parameters and frames (pickle files)
  iph, frame_files = read_input(args = sys.argv[1:])
  #frames = range(iph.frame_start, iph.frame_end)
  frames_best_100 = [198, 192, 208, 172, 41, 221, 684, 23, 246, 81, 382, 546, 446, 402, 308, 223, 53, 637, 343, 185, 140, 173, 211, 358, 160, 725, 346, 16, 169, 204, 293, 360, 170, 302, 516, 368, 189, 430, 432, 466, 183, 344, 148, 456, 476, 299, 244, 250, 231, 270, 119, 202, 251, 188, 200, 333, 574, 667, 59, 46, 245, 452, 284, 547, 92, 187, 328, 597, 696, 95, 209, 528, 711, 310, 156, 589, 641, 217, 233, 549, 721, 273, 353, 619, 48, 714, 228, 118, 491, 212, 214, 424, 289, 461, 659, 612, 268, 174, 523, 676]
  frames = frames_best_100[0:100]

  #1. prepare reference miller array
  if iph.file_name_ref_mtz == '':
    #if iso. ref. is not given, use the <I> to scale each frame.

    #calculate mean intensity for each frame and determine median of the
    #distribution
    def determine_mean_I_mproc_wrapper(arg):
      return determine_mean_I_mproc(arg, frame_files, iph)

    determine_mean_I_result = pool_map(
            args=frames,
            func=determine_mean_I_mproc_wrapper,
            processes=None)

    frames_mean_I = flex.double()
    for result in determine_mean_I_result:
      if result is not None:
        frames_mean_I.append(result[1])

    mean_of_mean_I = np.median(frames_mean_I)

    #use the calculate <mean_I> to scale each frame
    def scale_frame_by_mean_I_mproc_wrapper(arg):
      return scale_frame_by_mean_I_mproc(arg, frame_files, iph, mean_of_mean_I)

    scale_frame_by_mean_I_result = pool_map(
            args=frames,
            func=scale_frame_by_mean_I_mproc_wrapper,
            processes=None)

    observations_merge_mean_set = []
    observations_merge_mean_weight_set = []
    for result in scale_frame_by_mean_I_result:
      if result is not None:
        observations_merge_mean_set.append(result[1])
        observations_merge_mean_weight_set.append(result[2])

    from xfel.cxi.postrefine import merge_observations
    miller_array_ref, cc_merge_mean, slope_merge_mean, txt_merge_mean = merge_observations(observations_merge_mean_set,
        observations_merge_mean_weight_set, iph, iph.run_no+'/mean_scaled')

  else:
    miller_array_ref = iph.miller_array_iso

  #2. Post-refinement
  G_set = flex.double([1]*len(frame_files))
  B_factor_set = flex.double([0]*len(frame_files))
  rotx_set = flex.double([1]*len(frame_files))
  roty_set = flex.double([1]*len(frame_files))
  ry_set = flex.double([1]*len(frame_files))
  rz_set = flex.double([1]*len(frame_files))
  #2.1 Refine scale factors
  n_iters_scale = 1
  refine_mode = 'scale_factor'
  for i in range(n_iters_scale):
    def postrefine_by_frame_mproc_wrapper(arg):
      return postrefine_by_frame_mproc(arg, frame_files, iph, refine_mode, miller_array_ref,
          G_set, B_factor_set, rotx_set, roty_set, ry_set, rz_set)

    postrefine_by_frame_result = pool_map(
            args=frames,
            func=postrefine_by_frame_mproc_wrapper,
            processes=None)

    obs_pref_scale_set = []
    obs_pref_scale_weight_set = []
    r_int_pref_scale_set = flex.double()
    for result in postrefine_by_frame_result:
      if result is not None:
        frame_no = result[0]
        obs_pref_scale_set.append(result[1])
        obs_pref_scale_weight_set.append(result[2])
        cc_r_results = result[3]
        r_int_pref_scale_set.append(cc_r_results[5])
        scale_factors = result[4]
        G_set[frame_no] = scale_factors[0]
        B_factor_set[frame_no] = scale_factors[1]

    from xfel.cxi.postrefine import merge_observations
    miller_array_postref, cc_merge_postref, slope_merge_postref, txt_merge_postref = merge_observations(obs_pref_scale_set,
        obs_pref_scale_weight_set, iph, iph.run_no+'/postref_scale_cycle_'+str(i+1))

  #2.2 Refine crystal rotation (rotx, roty)
  n_iters_crystal_rotation = 1
  refine_mode = 'crystal_rotation'
  for i in range(n_iters_crystal_rotation):
    def postrefine_by_frame_mproc_wrapper(arg):
      return postrefine_by_frame_mproc(arg, frame_files, iph, refine_mode, miller_array_ref,
              G_set, B_factor_set, rotx_set, roty_set, ry_set, rz_set)

    postrefine_by_frame_result = pool_map(
            args=frames,
            func=postrefine_by_frame_mproc_wrapper,
            processes=None)

    obs_pref_rot_set = []
    obs_pref_rot_weight_set = []
    r_int_pref_rot_set = flex.double()
    for result in postrefine_by_frame_result:
      if result is not None:
        obs_pref_rot_set.append(result[1])
        obs_pref_rot_weight_set.append(result[2])
        cc_r_results = result[3]
        r_int_pref_rot_set.append(cc_r_results[5])
        rotations = result[4]
        rotx_set[frame_no] = rotations[0]
        roty_set[frame_no] = rotations[1]

    from xfel.cxi.postrefine import merge_observations
    miller_array_postref, cc_merge_postref, slope_merge_postref, txt_merge_postref = merge_observations(obs_pref_rot_set,
        obs_pref_rot_weight_set, iph, iph.run_no+'/postref_rotation_cycle_'+str(i+1))

  #2.3 Refine reflecting range (rs)
  n_iters_reflecting_range = 1
  refine_mode = 'reflecting_range'
  for i in range(n_iters_reflecting_range):
    def postrefine_by_frame_mproc_wrapper(arg):
      return postrefine_by_frame_mproc(arg, frame_files, iph, refine_mode, miller_array_ref,
              G_set, B_factor_set, rotx_set, roty_set, ry_set, rz_set)

    postrefine_by_frame_result = pool_map(
            args=frames,
            func=postrefine_by_frame_mproc_wrapper,
            processes=None)

    obs_pref_rs_set = []
    obs_pref_rs_weight_set = []
    r_int_pref_rs_set = flex.double()
    for result in postrefine_by_frame_result:
      if result is not None:
        obs_pref_rs_set.append(result[1])
        obs_pref_rs_weight_set.append(result[2])
        cc_r_results = result[3]
        r_int_pref_rs_set.append(cc_r_results[5])
        reflecting_ranges = result[4]
        ry_set[frame_no] = reflecting_ranges[0]
        rz_set[frame_no] = reflecting_ranges[1]


    from xfel.cxi.postrefine import merge_observations
    miller_array_postref, cc_merge_postref, slope_merge_postref, txt_merge_postref = merge_observations(obs_pref_rs_set,
        obs_pref_rs_weight_set, iph, iph.run_no+'/postref_reflecting_range_cycle_'+str(i+1))

  #2.4 Determine from Rint which refinement result to use in merging
  obs_sel_set = []
  obs_sel_weight_set = []
  pref_scale_final_code_set = flex.int([0]*len(frame_files))
  for i in range(len(frames)):
    r_int_all = flex.double([r_int_pref_scale_set[i], r_int_pref_rot_set[i], r_int_pref_rs_set[i]])
    perm = flex.sort_permutation(r_int_all)
    if perm[0] == 0:
      obs_sel_set.append(obs_pref_scale_set[i])
      obs_sel_weight_set.append(obs_pref_scale_weight_set[i])

    elif perm[0] == 1:
      obs_sel_set.append(obs_pref_rot_set[i])
      obs_sel_weight_set.append(obs_pref_rot_weight_set[i])
      pref_scale_final_code_set[frames[i]] = 1

    elif perm[0] == 2:
      obs_sel_set.append(obs_pref_rs_set[i])
      obs_sel_weight_set.append(obs_pref_rs_weight_set[i])
      pref_scale_final_code_set[frames[i]] = 2


  from xfel.cxi.postrefine import merge_observations
  miller_array_best, cc_merge_postref, slope_merge_postref, txt_merge_postref = merge_observations(obs_sel_set,
        obs_sel_weight_set, iph, iph.run_no+'/postref_best')

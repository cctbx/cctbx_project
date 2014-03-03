# LIBTBX_SET_DISPATCHER_NAME cxi.postrefine
from __future__ import division
import os, sys
from libtbx.easy_mp import pool_map
import numpy as np
from cctbx.array_family import flex
from datetime import datetime, time

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

  postref_results = prh.postrefine_by_frame(
        refine_mode, frame_files[frame_no], iph, miller_array_ref, scale_factors=(G_set[frame_no], B_factor_set[frame_no]),
        rotations=(rotx_set[frame_no], roty_set[frame_no]), reflecting_ranges=(ry_set[frame_no], rz_set[frame_no]))

  if postref_results is not None:
    pickle_filename, observations_refined, observations_weight, cc_r_results, lstsqr_results = postref_results
    return frame_no, observations_refined, observations_weight, cc_r_results, lstsqr_results
  else:
    return None


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
  frames = range(iph.frame_start, iph.frame_end)

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
    miller_array_mean_scaled, cc_merge_mean, slope_merge_mean, txt_merge_mean = merge_observations(observations_merge_mean_set,
        observations_merge_mean_weight_set, iph, iph.run_no+'/mean_scaled')
    miller_array_ref = miller_array_mean_scaled.generate_bijvoet_mates()
  else:
    miller_array_ref = iph.miller_array_iso
    txt_merge_mean = ''

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
  print 'Refine scale factors...'
  txt_merge_postref_scale = ''
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
    cn_result = 0
    for result in postrefine_by_frame_result:
      if result is not None:
        frame_no, obs_pref, obs_pref_weight, cc_r_results, scale_factors = result
        obs_pref_scale_set.append(obs_pref)
        obs_pref_scale_weight_set.append(obs_pref_weight)
        r_int_pref_scale_set.append(cc_r_results[5])
        G_set[frame_no] = scale_factors[0]
        B_factor_set[frame_no] = scale_factors[1]
        cn_result += 1
      else:
        obs_pref_scale_set.append(None)
        obs_pref_scale_weight_set.append(None)
        r_int_pref_scale_set.append(99)

    if cn_result > 0:
      from xfel.cxi.postrefine import merge_observations
      miller_array_postref_scale, cc_merge_postref, slope_merge_postref, txt_merge_postref_scale = merge_observations(
          obs_pref_scale_set, obs_pref_scale_weight_set, iph, iph.run_no+'/postref_scale_cycle_'+str(i+1))

    if iph.file_name_ref_mtz == '':
      miller_array_ref = miller_array_postref_scale.generate_bijvoet_mates()

  #2.2 Refine crystal rotation (rotx, roty)
  n_iters_crystal_rotation = 1
  refine_mode = 'crystal_rotation'
  print 'Refine crystal orientations...'
  txt_merge_postref_rot = ''
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
    cn_result = 0
    for result in postrefine_by_frame_result:
      if result is not None:
        frame_no, obs_pref, obs_pref_weight, cc_r_results, rotations = result
        obs_pref_rot_set.append(obs_pref)
        obs_pref_rot_weight_set.append(obs_pref_weight)
        r_int_pref_rot_set.append(cc_r_results[5])
        rotx_set[frame_no] = rotations[0]
        roty_set[frame_no] = rotations[1]
        cn_result += 1
      else:
        obs_pref_rot_set.append(None)
        obs_pref_rot_weight_set.append(None)
        r_int_pref_rot_set.append(99)

    if cn_result > 0:
      from xfel.cxi.postrefine import merge_observations
      miller_array_postref_rot, cc_merge_postref, slope_merge_postref, txt_merge_postref_rot = merge_observations(obs_pref_rot_set,
          obs_pref_rot_weight_set, iph, iph.run_no+'/postref_rotation_cycle_'+str(i+1))

    if iph.file_name_ref_mtz == '':
      miller_array_ref = miller_array_postref_rot.generate_bijvoet_mates()

  #2.3 Refine reflecting range (rs)
  n_iters_reflecting_range = 1
  refine_mode = 'reflecting_range'
  print 'Refine reflecting ranges...'
  txt_merge_postref_rs = ''
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
    cn_result = 0
    for result in postrefine_by_frame_result:
      if result is not None:
        frame_no, obs_pref, obs_pref_weight, cc_r_results, reflecting_ranges = result
        obs_pref_rs_set.append(obs_pref)
        obs_pref_rs_weight_set.append(obs_pref_weight)
        r_int_pref_rs_set.append(cc_r_results[5])
        ry_set[frame_no] = reflecting_ranges[0]
        rz_set[frame_no] = reflecting_ranges[1]
        cn_result += 1
      else:
        obs_pref_rs_set.append(None)
        obs_pref_rs_weight_set.append(None)
        r_int_pref_rs_set.append(99)

    if cn_result > 0:
      from xfel.cxi.postrefine import merge_observations
      miller_array_postref_rs, cc_merge_postref, slope_merge_postref, txt_merge_postref_rs = merge_observations(obs_pref_rs_set,
          obs_pref_rs_weight_set, iph, iph.run_no+'/postref_reflecting_range_cycle_'+str(i+1))

  #2.4 Determine from Rint which refinement result to use in merging
  txt_best = ''
  txt_merge_postref_best = ''
  if cn_result > 0:
    obs_sel_set = []
    obs_sel_weight_set = []
    pref_scale_final_code_set = flex.int([0]*len(frame_files))
    for i in range(len(frames)):
      if (r_int_pref_scale_set[i] == 99.0 and r_int_pref_rot_set[i] == 99.0 and r_int_pref_rs_set[i] == 99.0):
        txt_best += 'frame# %04d - not merged'%(frames[i])
      else:
        r_int_all = flex.double([r_int_pref_scale_set[i], r_int_pref_rot_set[i], r_int_pref_rs_set[i]])
        perm = flex.sort_permutation(r_int_all)
        r_int_all_sort = r_int_all.select(perm)
        if r_int_all_sort[0] <= iph.q_w_merge:
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
        else:
          txt_best += 'frame# %04d- not merged Rint=%6.4f'%(frames[i], r_int_all_sort[0])


    from xfel.cxi.postrefine import merge_observations
    miller_array_best, cc_merge_postref, slope_merge_postref, txt_merge_postref_best = merge_observations(obs_sel_set,
          obs_sel_weight_set, iph, iph.run_no+'/postref_best')
  else:
    print 'No frames refined'

  txt_out = iph.txt_out + txt_merge_mean + txt_merge_postref_scale + txt_merge_postref_rot + txt_merge_postref_rs + txt_merge_postref_best
  f = open(iph.run_no+'/log.txt', 'w')
  f.write(txt_out)
  f.close()

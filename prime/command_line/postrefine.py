from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME prime.postrefine
'''
Author      : Uervirojnangkoorn, M.
Created     : 7/13/2014
Description : Commands linked to prime.postrefine libraries.
'''
import os, sys
from libtbx.easy_mp import pool_map
import numpy as np
from cctbx.array_family import flex
from datetime import datetime, time
import math
from cctbx import miller
import matplotlib.pyplot as plt

def determine_mean_I_mproc(frame_no, frame_files, iparams, avg_mode):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  mean_I = prh.calc_mean_intensity(frame_files[frame_no], iparams, avg_mode)
  return mean_I

def scale_frame_by_mean_I_mproc(frame_no, frame_files, iparams, mean_of_mean_I, avg_mode):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  pres = prh.scale_frame_by_mean_I(frame_no,frame_files[frame_no], iparams, mean_of_mean_I, avg_mode)
  return pres

def postrefine_by_frame_mproc(frame_no, frame_files, iparams, miller_array_ref, pres_results, avg_mode):
  from prime.postrefine import postref_handler
  prh = postref_handler()
  if len(pres_results) == 0:
    pres_in = None
  else:
    pres_in = pres_results[frame_no]
  pres = prh.postrefine_by_frame(frame_no, frame_files[frame_no], iparams, miller_array_ref, pres_in, avg_mode)
  return pres

def read_pickles(data):
  frame_files = []
  for p in data:
    if os.path.isdir(p) == False:
      #check if list-of-pickle text file is given
      pickle_list_file = open(p,'r')
      pickle_list = pickle_list_file.read().split("\n")
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


if (__name__ == "__main__"):
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

  #0 .read input parameters and frames (pickle files)
  from prime.postrefine import read_input
  iparams, txt_out_input = read_input(sys.argv[:1])
  iparams.flag_volume_correction = False
  iparams.flag_apply_b_by_frame = True
  if iparams.flag_normalized:
    iparams.flag_apply_b_by_frame = False
    iparams.b_refine_d_min = -1
  print txt_out_input
  txt_out_verbose = 'Log verbose\n'+txt_out_input

  frame_files = read_pickles(iparams.data)
  frames = range(len(frame_files))

  #1. prepare reference miller array
  txt_merge_mean = 'Generating a reference set (will not be used if hklrefin is set)'
  print txt_merge_mean
  txt_merge_mean += '\n'
  txt_merge_mean_table = ''

  miller_array_ref = None
  avg_mode = 'average'
  #Always generate the mean-intensity scaled set.

  if iparams.flag_apply_b_by_frame:
    mean_of_mean_I = 0
  else:
    #Calculate <I> for each frame
    def determine_mean_I_mproc_wrapper(arg):
      return determine_mean_I_mproc(arg, frame_files, iparams, avg_mode)

    determine_mean_I_result = pool_map(
            iterable=frames,
            func=determine_mean_I_mproc_wrapper,
            processes=iparams.n_processors)

    frames_mean_I = flex.double()
    for result in determine_mean_I_result:
      if result is not None:
        mean_I, txt_out_result = result
        if mean_I is not None:
          frames_mean_I.append(mean_I)
    mean_of_mean_I = np.median(frames_mean_I)

  #use the calculate <mean_I> to scale each frame
  def scale_frame_by_mean_I_mproc_wrapper(arg):
    return scale_frame_by_mean_I_mproc(arg, frame_files, iparams, mean_of_mean_I, avg_mode)

  scale_frame_by_mean_I_result = pool_map(
          iterable=frames,
          func=scale_frame_by_mean_I_mproc_wrapper,
          processes=iparams.n_processors)

  observations_merge_mean_set = []
  for result in scale_frame_by_mean_I_result:
    if result is not None:
      pres, txt_out_result = result
      txt_merge_mean += txt_out_result
      if pres is not None:
        observations_merge_mean_set.append(pres)

  if len(observations_merge_mean_set) > 0:
    from prime.postrefine import prepare_output
    prep_output = prepare_output(observations_merge_mean_set, iparams, avg_mode)

    if prep_output is not None:
      cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort,  I_obs_all_sort, \
      sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort, rs_all_sort, wavelength_all_sort, \
      sin_sq_all_sort, SE_all_sort, uc_mean, wavelength_mean, pickle_filename_all_sort, txt_prep_out = prep_output

      if True:
        from prime.postrefine import calc_avg_I_cpp
        calc_average_I_result = calc_avg_I_cpp(cn_group, group_id_list, miller_indices_all_sort,
                                                 miller_indices_ori_all_sort, I_obs_all_sort,
                                                 sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort,
                                                 rs_all_sort, wavelength_all_sort, sin_sq_all_sort,
                                                 SE_all_sort, avg_mode, iparams, pickle_filename_all_sort)

        miller_index, I_avg, sigI_avg, stat, I_avg_two_halves_tuple, txt_obs_out, txt_reject_out = calc_average_I_result
        I_avg_even, I_avg_odd, I_avg_even_h, I_avg_odd_h, I_avg_even_k, I_avg_odd_k, I_avg_even_l, I_avg_odd_l = I_avg_two_halves_tuple
        miller_indices_merge = flex.miller_index()
        I_merge = flex.double()
        sigI_merge = flex.double()
        stat_all = []
        I_even = flex.double()
        I_odd = flex.double()
        I_even_h = flex.double()
        I_odd_h = flex.double()
        I_even_k = flex.double()
        I_odd_k = flex.double()
        I_even_l = flex.double()
        I_odd_l = flex.double()

        txt_out_verbose += 'Mean-scaled partiality-corrected set\n' + txt_obs_out
        txt_out_rejection = txt_reject_out
        for i in xrange(len(miller_index)):
          if math.isnan(stat[0][i]) or math.isinf(stat[0][i]) or math.isnan(stat[1][i]) or math.isinf(stat[1][i]):
            dummy = 0
          else:
            miller_indices_merge.append(miller_index[i])
            I_merge.append(I_avg[i])
            sigI_merge.append(sigI_avg[i])
            stat_all.append((stat[0][i],stat[1][i],stat[2][i],stat[3][i],stat[4][i]))
            I_even.append(I_avg_even[i])
            I_odd.append(I_avg_odd[i])
            I_even_h.append(I_avg_even_h[i])
            I_odd_h.append(I_avg_odd_h[i])
            I_even_k.append(I_avg_even_k[i])
            I_odd_k.append(I_avg_odd_k[i])
            I_even_l.append(I_avg_even_l[i])
            I_odd_l.append(I_avg_odd_l[i])


      f = open(iparams.run_no+'/rejections.txt', 'w')
      f.write(txt_out_rejection)
      f.close()

      from prime.postrefine import write_output
      miller_array_ref, txt_merge_mean_table = write_output(miller_indices_merge,
                                                                      I_merge, sigI_merge,
                                                                      stat_all, (I_even, I_odd, I_even_h, I_odd_h,
                                                                      I_even_k, I_odd_k, I_even_l, I_odd_l),
                                                                      iparams, uc_mean,
                                                                      wavelength_mean,
                                                                      iparams.run_no+'/mean_scaled',
                                                                      avg_mode)

      txt_merge_mean +=  txt_merge_mean_table + txt_prep_out
      print txt_merge_mean_table
      print txt_prep_out
      if iparams.flag_force_no_postrefine:
        txt_out = txt_out_input + txt_merge_mean
        f = open(iparams.run_no+'/log.txt', 'w')
        f.write(txt_out)
        f.close()

        if iparams.flag_output_verbose:
          f = open(iparams.run_no+'/log_verbose.txt', 'w')
          f.write(txt_out_verbose)
          f.close()
        exit()

  if iparams.hklrefin is not None:
    #In case ref. intensity set is given, overwrite the reference set.
    from iotbx import reflection_file_reader
    reflection_file_ref = reflection_file_reader.any_reflection_file(iparams.hklrefin)
    miller_arrays_ref=reflection_file_ref.as_miller_arrays()
    is_found_ref_as_intensity_array = False
    is_found_ref_as_amplitude_array = False
    for miller_array in miller_arrays_ref:
      if miller_array.is_xray_intensity_array():
        miller_array_ref = miller_array.deep_copy()
        is_found_ref_as_intensity_array = True
        break
      elif miller_array.is_xray_amplitude_array():
        is_found_ref_as_amplitude_array = True
        miller_array_converted_to_intensity = miller_array.as_intensity_array()
    if is_found_ref_as_intensity_array == False:
      if is_found_ref_as_amplitude_array:
        miller_array_ref = miller_array_converted_to_intensity.deep_copy()

  if miller_array_ref is None:
    print 'No reference intensity. Exit without post-refinement'
    exit()

  if iparams.flag_plot:
    #collect mean intensity of the reference set
    obs_ref = miller_array_ref.deep_copy()
    binner = obs_ref.setup_binner(n_bins=20)
    binner_indices = binner.bin_indices()
    obs_ref_mean = obs_ref.mean(use_binning=True)
    I_ref_mean = flex.double([0]*binner.n_bins_used())
    for i_bin in range(binner.n_bins_used()):
      if obs_ref_mean.data[i_bin+1] is not None:
        I_ref_mean[i_bin] = obs_ref_mean.data[i_bin+1]
    I_raw_set = []
    for pres in observations_merge_mean_set:
      obs = pres.observations.deep_copy()
      obs.use_binning_of(obs_ref)
      obs_mean = obs.mean(use_binning=True)
      I_mean = flex.double([0]*binner.n_bins_used())
      one_dsqr = flex.double()
      for i_bin in range(binner.n_bins_used()):
        one_dsqr.append(1/binner.bin_d_range(i_bin+1)[1]**2)
        if obs_mean.data[i_bin+1] is not None:
          I_mean[i_bin] = obs_mean.data[i_bin+1]
      I_raw_set.append(I_mean)

    for I_mean,cn_plot in zip(I_raw_set, range(len(I_raw_set))):
      plt.plot(one_dsqr, flex.log(I_mean), linestyle='-', linewidth=2.0, c='b')
    plt.title('Wilson plot before scaling')
    plt.xlabel('1/(d^2)')
    plt.ylabel('Log(<I>)')
    plt.grid()
    plt.show()

  #2. Generate Test Set
  ma_ref_r_free_flags = None
  if iparams.fraction_r_free > 0:
    ma_ref_r_free_flags = miller_array_ref.generate_r_free_flags(fraction=iparams.fraction_r_free, n_shells=iparams.n_bins)
    ma_indices_for_test = miller_array_ref.indices().select(ma_ref_r_free_flags.data())
    txt_r_free_out = ''
    for i in range(len(ma_indices_for_test)):
      txt_r_free_out += str(ma_indices_for_test[i][0])+' '+str(ma_indices_for_test[i][1])+' '+str(ma_indices_for_test[i][2])+'\n'
    f = open(iparams.run_no+'/rfree.txt', 'w')
    f.write(txt_r_free_out)
    f.close()

  #3. Post-refinement
  n_iters = iparams.n_postref_cycle
  txt_merge_postref = ''
  postrefine_by_frame_result = None
  postrefine_by_frame_pres_list = []
  for i_iter in range(n_iters):
    if i_iter == (n_iters-1):
      avg_mode = 'final'
    else:
      avg_mode = 'weighted'

    if iparams.postref.reflecting_range.flag_on and i_iter > 0:
      iparams.b_refine_d_min = -1

    _txt_merge_postref = 'Post-refinement cycle '+str(i_iter+1)+' ('+avg_mode+')\n'
    _txt_merge_postref += '* R and CC show percent change.\n'
    txt_merge_postref += _txt_merge_postref
    print _txt_merge_postref

    def postrefine_by_frame_mproc_wrapper(arg):
      return postrefine_by_frame_mproc(arg, frame_files, iparams,
                                       miller_array_ref, postrefine_by_frame_pres_list, avg_mode)

    postrefine_by_frame_result = pool_map(
            iterable=frames,
            func=postrefine_by_frame_mproc_wrapper,
            processes=iparams.n_processors)

    postrefine_by_frame_good = []
    postrefine_by_frame_pres_list = []
    for results in postrefine_by_frame_result:
      if results is not None:
        pres, txt_out_result = results
        txt_merge_postref += txt_out_result
        postrefine_by_frame_pres_list.append(pres)
        if pres is not None:
          postrefine_by_frame_good.append(pres)
      else:
        postrefine_by_frame_pres_list.append(None)

    if len(postrefine_by_frame_good) > 0:

      from prime.postrefine import prepare_output
      prep_output = prepare_output(postrefine_by_frame_good, iparams, avg_mode)

      if prep_output is not None:
        cn_group, group_id_list, miller_indices_all_sort, miller_indices_ori_all_sort, I_obs_all_sort, \
        sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort, rs_all_sort, wavelength_all_sort, \
        sin_sq_all_sort, SE_all_sort, uc_mean, wavelength_mean, pickle_filename_all_sort, txt_prep_out = prep_output

        if True:
          from prime.postrefine import calc_avg_I_cpp
          calc_average_I_result = calc_avg_I_cpp(cn_group, group_id_list, miller_indices_all_sort,
                                                   miller_indices_ori_all_sort, I_obs_all_sort,
                                                   sigI_obs_all_sort, G_all_sort, B_all_sort, p_all_sort,
                                                   rs_all_sort, wavelength_all_sort, sin_sq_all_sort,
                                                   SE_all_sort, avg_mode, iparams, pickle_filename_all_sort)

          miller_index, I_avg, sigI_avg, stat, I_avg_two_halves_tuple, txt_obs_out, txt_reject_out = calc_average_I_result
          I_avg_even, I_avg_odd, I_avg_even_h, I_avg_odd_h, I_avg_even_k, I_avg_odd_k, I_avg_even_l, I_avg_odd_l = I_avg_two_halves_tuple
          miller_indices_merge = flex.miller_index()
          I_merge = flex.double()
          sigI_merge = flex.double()
          stat_all = []
          I_even = flex.double()
          I_odd = flex.double()
          I_even_h = flex.double()
          I_odd_h = flex.double()
          I_even_k = flex.double()
          I_odd_k = flex.double()
          I_even_l = flex.double()
          I_odd_l = flex.double()

          txt_out_verbose += 'Post-refined merged set (cycle '+str(i_iter+1)+')\n' + txt_obs_out
          txt_out_rejection = txt_reject_out
          for i in xrange(len(miller_index)):
            if math.isnan(stat[0][i]) or math.isinf(stat[0][i]) or math.isnan(stat[1][i]) or math.isinf(stat[1][i]):
              dummy = 0
            else:
              miller_indices_merge.append(miller_index[i])
              I_merge.append(I_avg[i])
              sigI_merge.append(sigI_avg[i])
              stat_all.append((stat[0][i],stat[1][i],stat[2][i],stat[3][i],stat[4][i]))
              I_even.append(I_avg_even[i])
              I_odd.append(I_avg_odd[i])
              I_even_h.append(I_avg_even_h[i])
              I_odd_h.append(I_avg_odd_h[i])
              I_even_k.append(I_avg_even_k[i])
              I_odd_k.append(I_avg_odd_k[i])
              I_even_l.append(I_avg_even_l[i])
              I_odd_l.append(I_avg_odd_l[i])


        f = open(iparams.run_no+'/rejections.txt', 'a')
        f.write(txt_out_rejection)
        f.close()

        from prime.postrefine import write_output
        miller_array_pr, txt_merge_out = write_output(miller_indices_merge,
                                                               I_merge, sigI_merge, stat_all,
                                                               (I_even, I_odd, I_even_h, I_odd_h,
                                                               I_even_k, I_odd_k, I_even_l, I_odd_l), iparams, uc_mean,
                                                               wavelength_mean,
                                                               iparams.run_no+'/postref_cycle_'+str(i_iter+1), avg_mode)

        #calculate R and Rfree
        R_work, R_free = (0,0)
        if ma_ref_r_free_flags is not None:
          ma_ref_r_work_indices = ma_ref_r_free_flags.indices().select(ma_ref_r_free_flags.data()==False)
          ma_ref_r_free_indices = ma_ref_r_free_flags.indices().select(ma_ref_r_free_flags.data())

          matches = miller.match_multi_indices(
                  miller_indices_unique=miller_array_ref.indices(),
                  miller_indices=miller_array_pr.indices())
          I_ref_match = flex.double([miller_array_ref.data()[pair[0]] for pair in matches.pairs()])
          I_pr_match = flex.double([miller_array_pr.data()[pair[1]] for pair in matches.pairs()])
          miller_indices_match = flex.miller_index([miller_array_pr.indices()[pair[1]] for pair in matches.pairs()])

          matches = miller.match_multi_indices(
                  miller_indices_unique=miller_indices_match,
                  miller_indices=ma_ref_r_work_indices)
          I_ref_work = flex.double([I_ref_match[pair[0]] for pair in matches.pairs()])
          I_pr_work = flex.double([I_pr_match[pair[0]] for pair in matches.pairs()])

          matches = miller.match_multi_indices(
                  miller_indices_unique=miller_indices_match,
                  miller_indices=ma_ref_r_free_indices)
          I_ref_free = flex.double([I_ref_match[pair[0]] for pair in matches.pairs()])
          I_pr_free = flex.double([I_pr_match[pair[0]] for pair in matches.pairs()])

          R_work = math.sqrt(np.sum((I_ref_work - I_pr_work)**2)/np.sum(I_ref_work**2))
          R_free = math.sqrt(np.sum((I_ref_free - I_pr_free)**2)/np.sum(I_ref_free**2))
          txt_R_out =  ' Rwork:                   %12.2f%% REFL: %6.2f%%\n'%(R_work*100, (len(I_ref_work)/len(miller_indices_match))*100)
          txt_R_out += ' Rfree:                   %12.2f%% REFL: %6.2f%%\n'%(R_free*100, (len(I_ref_free)/len(miller_indices_match))*100)

        miller_array_ref = miller_array_pr.deep_copy()

        txt_merge_postref +=  txt_merge_out + txt_prep_out + txt_R_out
        print txt_merge_out
        print txt_prep_out
        print txt_R_out
    else:
      print "No frames merged as a reference set - exit without post-refinement"
      exit()

  #collect caculating time
  time_global_end=datetime.now()
  time_global_spent=time_global_end-time_global_start
  txt_out_time_spent = 'Total calculation time: '+'{0:.2f}'.format(time_global_spent.seconds)+' seconds\n'
  print txt_out_time_spent

  txt_out = txt_out_input + txt_merge_mean + txt_merge_postref + txt_out_time_spent
  f = open(iparams.run_no+'/log.txt', 'w')
  f.write(txt_out)
  f.close()

  if iparams.flag_output_verbose:
    f = open(iparams.run_no+'/log_verbose.txt', 'w')
    f.write(txt_out_verbose)
    f.close()

  #remove rejections.txt
  if os.path.isfile(iparams.run_no+'/rejections.txt'):
    os.remove(iparams.run_no+'/rejections.txt')

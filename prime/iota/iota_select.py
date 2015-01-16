from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 01/16/2015
Description : IOTA pickle selection module. Selects the best integration results from a
              set of pickles derived from a single image.
'''

import os
import numpy as np
import logging

import os,cPickle as pickle
from xfel.clustering.cluster import Cluster
from xfel.clustering.singleframe import SingleFrame

# Selects only integrated pickles that fit sg / uc parameters specified in the .phil
# file. Also checks that the low-res cutoff is beyond 10A.
def prefilter(gs_params, total_tmp_pickles):

  acceptable_pickles = []
  if gs_params.flag_prefilter == True:
    for tmp_pickle in total_tmp_pickles:

      tmp_pickle_name = os.path.basename(tmp_pickle)
      uc_tol = gs_params.target_uc_tolerance
      user_uc = []

      for prm in gs_params.target_unit_cell.parameters():
        user_uc.append(prm)

      # read pickle info and determine uc differences
      observations = SingleFrame(tmp_pickle, tmp_pickle_name).miller_array
      p_pg = observations.space_group_info()
      p_uc = observations.unit_cell().parameters()

      delta_a = abs(p_uc[0] - user_uc[0])
      delta_b = abs(p_uc[1] - user_uc[1])
      delta_c = abs(p_uc[2] - user_uc[2])
      delta_alpha = abs(p_uc[3] - user_uc[3])
      delta_beta = abs(p_uc[4] - user_uc[4])
      delta_gamma = abs(p_uc[5] - user_uc[5])

      # Determine if pickle satisfies sg / uc parameters within given
      # tolerance and low resolution cutoff
      if str(p_pg) == gs_params.target_pointgroup:
        if  (delta_a <= user_uc[0] * uc_tol and
              delta_b <= user_uc[1] * uc_tol and
              delta_c <= user_uc[2] * uc_tol and
              delta_alpha <= user_uc[3] * uc_tol and
              delta_beta <= user_uc[4] * uc_tol and
              delta_gamma <= user_uc[5] * uc_tol):
          acceptable_pickles.append(tmp_pickle)

  return acceptable_pickles


# Select integrated pickle with the most reflections with I / sigI over a specified limit
def best_by_strong(gs_params, acceptable_pickles):

  sel_pickle_list = []
  for pickle in acceptable_pickles:
    observations = SingleFrame(pickle, os.path.split(pickle)[1]).miller_array
    rl_observations = observations.resolution_filter(gs_params.selection_res_limit.d_max, gs_params.selection_res_limit.d_min)
    I_over_sigI = rl_observations.data() / rl_observations.sigmas()
    num_strong_obs = len([val for val in I_over_sigI if val >= gs_params.min_sigma])
    sel_pickle_list.append(num_strong_obs)

  best_file = acceptable_pickles[sel_pickle_list.index(max(sel_pickle_list))]

  return best_file


# Select integrated pickle with the most reflections within a set resolution limit
def best_by_reflections(gs_params, acceptable_pickles):

  sel_pickle_list = []
  for pickle in acceptable_pickles:
    observations = SingleFrame(pickle, os.path.split(pickle)[1]).miller_array
    rl_observations = observations.resolution_filter(gs_params.selection_res_limit.d_max, gs_params.selection_res_limit.d_min)
    sel_pickle_list.append(len(rl_observations.data()))

  best_file = acceptable_pickles[sel_pickle_list.index(max(sel_pickle_list))]

  return best_file

# Select integrated pickle with the closest unit cell to target
def best_by_uc(gs_params, acceptable_pickles):

  sel_pickle_list = []
  uc_tol = gs_params.target_uc_tolerance
  user_uc = []

  for prm in gs_params.target_unit_cell.parameters():
    user_uc.append(prm)

  for pickle in acceptable_pickles:
    p_uc = SingleFrame(pickle, os.path.split(pickle)[1]).miller_array.unit_cell().parameters()
    delta_a = abs(p_uc[0] - user_uc[0]) / user_uc[0]
    delta_b = abs(p_uc[1] - user_uc[1]) / user_uc[1]
    delta_c = abs(p_uc[2] - user_uc[2]) / user_uc[2]
    delta_alpha = abs(p_uc[3] - user_uc[3]) / user_uc[3]
    delta_beta = abs(p_uc[4] - user_uc[4]) / user_uc[4]
    delta_gamma = abs(p_uc[5] - user_uc[5]) / user_uc[5]
    uc_distance = np.mean([delta_a, delta_b, delta_c, delta_alpha, delta_beta, delta_gamma])

    sel_pickle_list.append(uc_distance)

  best_file = acceptable_pickles[sel_pickle_list.index(min(sel_pickle_list))]

  return best_file


# Select integrated pickle with the lowest x,y offset
def best_by_offset(acceptable_pickles):

  pickle_cluster = Cluster.from_files(acceptable_pickles)
  best_file = min(pickle_cluster.members, key=lambda im: im.spot_offset).path

  return best_file

# Main selection module. Looks through integrated pickles in a specified folder and
# copies the best ones to a file. Outputs a list to log file and marks the selected
# pickle file.
def best_file_selection(gs_params, output_entry, log_dir):

  ps_logger = logging.getLogger('ps_log')

  input_file = output_entry[0]
  abs_tmp_dir = output_entry[1]

  total_tmp_pickles = [os.path.join(abs_tmp_dir, tmp_pickle) for tmp_pickle in os.listdir(abs_tmp_dir) if ".pickle" in tmp_pickle]

  # apply prefilter if specified and make a list of acceptable pickles
  if len(total_tmp_pickles) == 0:
    ps_logger.info('No integrated images found ' \
                    'in {}:\n'.format(abs_tmp_dir))
    with open('{}/not_integrated.lst'.format(os.path.abspath(gs_params.output)), 'a') as no_int:
      no_int.write('{}\n'.format(input_file))
    acceptable_pickles = []

  else:
    acceptable_pickles = prefilter(gs_params, total_tmp_pickles)
    if len(acceptable_pickles) == 0:
      ps_logger.info('Discarded all {0} integrated pickles ' \
                      'in {1}:\n'.format(len(total_tmp_pickles), abs_tmp_dir))

      with open('{}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output)), 'a') as bad_int:
        bad_int.write('{}\n'.format(input_file))


    else:
      # Selection and copying of pickles, output of stats to log file
      ps_logger.info('Selecting from {0} out '\
                      'of {1} integrated pickles ' \
                      'in {2}:\n'.format(len(acceptable_pickles),
                      len(total_tmp_pickles), abs_tmp_dir))
      filename = str(os.path.split(acceptable_pickles[0])[1])
      categories = ' {:^{pwidth}}{:^16}{:^15}{:^45}' \
                    '{:^12}{:^10}'.format('Filename',
                    'resolution', 's.g.', 'unit cell',
                    'spots', 'strong', pwidth=len(filename)+5)
      line = ' {:-^{pwidth}}{:-^16}{:-^15}{:-^45}' \
              '{:-^12}{:-^10}'.format('', '', '', '',
               '', '', pwidth=len(filename)+5)
      ps_logger.info(categories)
      ps_logger.info(line)

      # Report pickle stats. Mark selected pickle with asterisk for posterity
      for pickle in acceptable_pickles:
        pickle_name = os.path.split(pickle)[1]

        observations = SingleFrame(pickle,
                        os.path.split(pickle)[1]).miller_array
        res = observations.d_max_min()
        pg = observations.space_group_info()
        uc = observations.unit_cell().parameters()
        rl_observations = observations.resolution_filter(gs_params.selection_res_limit.d_max, gs_params.selection_res_limit.d_min)
        I_over_sigI = rl_observations.data() / rl_observations.sigmas()
        ref = len(rl_observations.data())
        sref = len([val for val in I_over_sigI if val >= gs_params.min_sigma])

        info_line = '  {:<{pwidth}}{:>7.2f} - {:<5.2f}{:^15}{:>6.2f}, ' \
                     '{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                     '{:>6.2f}{:^12}{:^10}'.format(pickle_name, res[0],
                     res[1], pg, uc[0], uc[1], uc[2], uc[3], uc[4], uc[5],
                     ref, sref, pwidth=len(filename)+5)
        ps_logger.info(info_line)


      # Make selections & copy files
      selected_info = []

      # Total reflections
      sel_pickle = best_by_reflections(gs_params, acceptable_pickles)
      with open('{}/best_by_total.lst'.format(os.path.abspath(gs_params.output)), 'a') as best_file:
        best_file.write('{}\n'.format(sel_pickle))
      selected_info.append(['T', sel_pickle, os.path.split(sel_pickle)[1]])

      # Strong reflections
      sel_pickle = best_by_strong(gs_params, acceptable_pickles)
      with open('{}/best_by_strong.lst'.format(os.path.abspath(gs_params.output)), 'a') as best_file:
        best_file.write('{}\n'.format(sel_pickle))
      selected_info.append(['S', sel_pickle, os.path.split(sel_pickle)[1]])

      # Unit cell
      sel_pickle = best_by_uc(gs_params, acceptable_pickles)
      with open('{}/best_by_uc.lst'.format(os.path.abspath(gs_params.output)), 'a') as best_file:
        best_file.write('{}\n'.format(sel_pickle))
      selected_info.append(['U', sel_pickle, os.path.split(sel_pickle)[1]])

      # x,y offset
      sel_pickle = best_by_offset(acceptable_pickles)
      with open('{}/best_by_offset.lst'.format(os.path.abspath(gs_params.output)), 'a') as best_file:
        best_file.write('{}\n'.format(sel_pickle))
      selected_info.append(['O', sel_pickle, os.path.split(sel_pickle)[1]])


      # Output selected file information
      ps_logger.info('\nSelected:')
      for sel in selected_info:
        observations = SingleFrame(sel[1], sel[2]).miller_array
        res = observations.d_max_min()
        pg = observations.space_group_info()
        uc = observations.unit_cell().parameters()

        rl_observations = observations.resolution_filter(gs_params.selection_res_limit.d_max, gs_params.selection_res_limit.d_min)
        I_over_sigI = rl_observations.data() / rl_observations.sigmas()
        ref = len(rl_observations.data())
        sref = len([val for val in I_over_sigI if val >= gs_params.min_sigma])

        info_line = '{} {:<{pwidth}}{:>7.2f} - {:<5.2f}{:^15}{:>6.2f}, ' \
                    '{:>6.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}, ' \
                    '{:>6.2f}{:^12}{:^10}'.format(sel[0], sel[2], res[0],
                    res[1], pg, uc[0], uc[1], uc[2], uc[3], uc[4], uc[5],
                    ref, sref, pwidth=len(filename)+5)
        ps_logger.info(info_line)

    ps_logger.info('\n')

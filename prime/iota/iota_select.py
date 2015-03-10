from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 03/06/2015
Description : IOTA pickle selection module. Selects the best integration results from a
              set of pickles derived from a single image.
'''

import os
from prime.iota.iota_input import main_log

import os,cPickle as pickle
import dials.util.command_line as cmd


# Selects only integrated pickles that fit sg / uc parameters specified in the .phil
# file. Also checks that the low-res cutoff is beyond 10A.
def prefilter(gs_params, int_list):

  acceptable_results = []
  if gs_params.flag_prefilter == True:
    for int_entry in int_list:

      tmp_pickle_name = os.path.basename(int_entry[0])
      uc_tol = gs_params.target_uc_tolerance
      user_uc = []

      for prm in gs_params.target_unit_cell.parameters():
        user_uc.append(prm)

      # read integration info and determine uc differences
      p_pg = int_entry[3].replace(" ","")
      p_uc = int_entry[4], int_entry[5], int_entry[6],\
             int_entry[7], int_entry[8], int_entry[9]

      delta_a = abs(p_uc[0] - user_uc[0])
      delta_b = abs(p_uc[1] - user_uc[1])
      delta_c = abs(p_uc[2] - user_uc[2])
      delta_alpha = abs(p_uc[3] - user_uc[3])
      delta_beta = abs(p_uc[4] - user_uc[4])
      delta_gamma = abs(p_uc[5] - user_uc[5])

      # Determine if pickle satisfies sg / uc parameters within given
      # tolerance and low resolution cutoff
      if p_pg == gs_params.target_pointgroup.replace(" ",""):
        if  (delta_a <= user_uc[0] * uc_tol and
              delta_b <= user_uc[1] * uc_tol and
              delta_c <= user_uc[2] * uc_tol and
              delta_alpha <= user_uc[3] * uc_tol and
              delta_beta <= user_uc[4] * uc_tol and
              delta_gamma <= user_uc[5] * uc_tol):
          acceptable_results.append(int_entry)
  else:
    acceptable_results = total_tmp_pickles

  return acceptable_results


# Main selection module. Looks through integrated pickles in a specified folder and
# copies the best ones to a file. Outputs a list to log file and marks the selected
# pickle file.
def best_file_selection(gs_params, output_entry, log_dir, n_int):

  logfile = '{}/iota.log'.format(log_dir)

  abs_tmp_dir = output_entry[0]
  input_file = output_entry[1]
  result_file = os.path.join(abs_tmp_dir, output_entry[2])
  ps_log_output = []
  selection_result = []

  if not os.path.isfile(result_file):
    return

# Read integration record file and convert to list
  with open(result_file, 'r') as int_results:
    int_content = int_results.read()

  prelim_int_list = [int_res.split(',') for int_res in int_content.splitlines()]
  int_list = []
  for item in prelim_int_list:
    result_line = [item[0], int(item[1]), int(item[2]), item[3],
                   float(item[4]), float(item[5]), float(item[6]),
                   float(item[7]), float(item[8]), float(item[9]),
                   int(item[10]), float(item[11]), float(item[12]),
                   float(item[13])]
    int_list.append(result_line)


# apply prefilter if specified and make a list of acceptable pickles
  if not os.path.isfile(result_file):
    ps_log_output.append('No integrated images found ' \
                    'in {}:\n'.format(abs_tmp_dir))
    acceptable_results = []
    int_summary ='{} --     not integrated'.format(input_file)

    with open('{}/not_integrated.lst'.format(os.path.abspath(gs_params.output)), 'a') as no_int:
      no_int.write('{}\n'.format(input_file))

  else:
    acceptable_results = prefilter(gs_params, int_list)
    if len(acceptable_results) == 0:
      ps_log_output.append('Discarded all {0} integrated pickles ' \
                      'in {1}:\n'.format(len(int_list), abs_tmp_dir))
      with open('{}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output)), 'a') as bad_int:
        bad_int.write('{}\n'.format(input_file))

      int_summary ='{} --     failed prefilter'.format(input_file)


    else:
      # Selection and copying of pickles, output of stats to log file
      ps_log_output.append('Selecting from {0} out '\
                      'of {1} integration results for ' \
                      '{2}:\n'.format(len(acceptable_results),
                      len(int_list), input_file))
      categories = '{:^4}{:^4}{:^9}{:^8}{:^55}{:^12}{:^12}{:^12}'\
                   ''.format('H', 'A', 'RES', 'SG.',
                   'UNIT CELL', 'SPOTS', 'MOS', 'MQ')
      line = '{:-^4}{:-^4}{:-^9}{:-^8}{:-^55}{:-^12}{:-^12}{:^12}'\
             ''.format('', '', '', '', '','', '', '')
      ps_log_output.append(categories)
      ps_log_output.append(line)

      int_summary = '{} -- {:>3} successful integration '\
                         'results'.format(input_file, len(acceptable_results))


      mqs = []
      for acc in acceptable_results:
        cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
               ''.format(acc[4], acc[5], acc[6], acc[7], acc[8], acc[9])
        info_line = '{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^12.2f}{:^12.2f}'\
                    ''.format(acc[1], acc[2], acc[11], acc[3], cell,
                              acc[10], acc[12], acc[13])
        ps_log_output.append(info_line)
        #mosaicities.append(float(acc[12]))
        mqs.append(float(acc[13]))

      best = acceptable_results[mqs.index(min(mqs))]
      with open('{}/selected.lst'.format(os.path.abspath(gs_params.output)), \
                                                            'a') as sel_int:
        sel_int.write('{}\n'.format(input_file))

      selection_result = [input_file, best[1], best[1], best[2]]

      # Output selected file information
      ps_log_output.append('\nSelected:')

      cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
             ''.format(best[4], best[5], best[6], best[7], best[8], best[9])
      info_line = '{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^12.2f}{:^12.2f}'\
                  ''.format(best[1], best[2], best[11], best[3], cell,
                            best[10], best[12], best[13])
      ps_log_output.append(info_line)

    ps_log_output.append('\n')
    main_log(logfile, '\n'.join(ps_log_output))

  # Progress bar for selection
  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'a')) as prog_log:
    prog_log.write("{}\n".format(int_summary))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()
    prog_count = len(prog_content.splitlines())

  gs_prog = cmd.ProgressBar(title='PICKLE SELECTION', estimate_time=False, spinner=False)
  if prog_count >= n_int:
    gs_prog.finished()
  else:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)

  return selection_result

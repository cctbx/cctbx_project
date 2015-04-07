from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 04/01/2015
Description : IOTA pickle selection module. Selects the best integration results
              from grid search output.
'''

import os, sys, traceback
from prime.iota.iota_input import main_log
import numpy as np
import csv

import dials.util.command_line as cmd


# Selects only integrated pickles that fit sg / uc parameters specified in the .phil
# file. Also checks that the low-res cutoff is beyond 10A.
def prefilter(gs_params, int_list):

  acceptable_results = []
  if gs_params.flag_prefilter == True:
    for i in int_list:
      try:
        tmp_pickle_name = os.path.basename(i['img'])
        uc_tol = gs_params.target_uc_tolerance
        user_uc = [prm for prm in gs_params.target_unit_cell.parameters()]

        # read integration info and determine uc differences
        p_pg = i['sg'].replace(" ","")


        delta_a = abs(i['a'] - user_uc[0])
        delta_b = abs(i['b'] - user_uc[1])
        delta_c = abs(i['c'] - user_uc[2])
        delta_alpha = abs(i['alpha'] - user_uc[3])
        delta_beta = abs(i['beta'] - user_uc[4])
        delta_gamma = abs(i['gamma'] - user_uc[5])

        # Determine if pickle satisfies sg / uc parameters within given
        # tolerance and low resolution cutoff
        if p_pg == gs_params.target_pointgroup.replace(" ",""):
          if  (delta_a <= user_uc[0] * uc_tol and
                delta_b <= user_uc[1] * uc_tol and
                delta_c <= user_uc[2] * uc_tol and
                delta_alpha <= user_uc[3] * uc_tol and
                delta_beta <= user_uc[4] * uc_tol and
                delta_gamma <= user_uc[5] * uc_tol):
            acceptable_results.append(i)
      except TypeError:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))
  else:
    acceptable_results = int_list

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

# apply prefilter if specified and make a list of acceptable pickles
  if not os.path.isfile(result_file):
    ps_log_output.append('No integrated images found ' \
                    'in {}:\n'.format(abs_tmp_dir))
    acceptable_results = []
    int_summary ='{} --     not integrated'.format(input_file)

    with open('{}/not_integrated.lst'\
              ''.format(os.path.abspath(gs_params.output)), 'a') as no_int:
      no_int.write('{}\n'.format(input_file))

  else:

    with open(result_file, 'rb') as rf:
      reader = csv.DictReader(rf)
      int_list = [i for i in list(reader) if i['img'] != '' and i['a'] != None]

    types = [('img', str), ('sih', int), ('sph', int), ('spa', int),\
             ('sg', str), ('a', float), ('b', float), ('c', float),\
             ('alpha', float), ('beta', float), ('gamma', float),\
             ('strong', int), ('res', float), ('mos', float), ('mq', float)]

    for i in int_list:
      try:
        i.update((key, conv(i[key])) for (key, conv) in types)
      except ValueError:
        raise Exception("".join(traceback.format_exception(*sys.exc_info())))


    acceptable_results = prefilter(gs_params, int_list)
    if len(acceptable_results) == 0:
      ps_log_output.append('All {0} entries in {1} failed prefilter' \
                           '\n'.format(len(int_list), result_file))
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


      for acc in acceptable_results:
        cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
               ''.format(acc['a'], acc['b'], acc['c'],
                         acc['alpha'], acc['beta'], acc['gamma'])
        info_line = '{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^12.4f}{:^12.2f}'\
                    ''.format(acc['sph'], acc['spa'], acc['res'], acc['sg'],
                              cell, acc['strong'], acc['mos'], acc['mq'])
        ps_log_output.append(info_line)

      # Compute average values
      avg_res = np.mean([item['res'] for item in acceptable_results])
      avg_spots = np.mean([item['strong'] for item in acceptable_results])
      avg_mos = np.mean([item['mos'] for item in acceptable_results])
      avg_mq = np.mean([item['mq'] for item in acceptable_results])
      avg_cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
              ''.format(np.mean([item['a'] for item in acceptable_results]),
                        np.mean([item['b'] for item in acceptable_results]),
                        np.mean([item['c'] for item in acceptable_results]),
                        np.mean([item['alpha'] for item in acceptable_results]),
                        np.mean([item['beta'] for item in acceptable_results]),
                        np.mean([item['gamma'] for item in acceptable_results]))

      info_line = '\nAVG:    {:^9.2f}{:^8}{:^55}{:^12.2}{:^12.4f}{:^12.2f}'\
                    ''.format(avg_res, '', avg_cell, avg_spots, avg_mos, avg_mq)
      ps_log_output.append(info_line)

      # Compute standard deviations
      std_res = np.std([item['res'] for item in acceptable_results])
      std_spots = np.std([item['strong'] for item in acceptable_results])
      std_mos = np.std([item['mos'] for item in acceptable_results])
      std_mq = np.std([item['mq'] for item in acceptable_results])
      std_cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
              ''.format(np.std([item['a'] for item in acceptable_results]),
                        np.std([item['b'] for item in acceptable_results]),
                        np.std([item['c'] for item in acceptable_results]),
                        np.std([item['alpha'] for item in acceptable_results]),
                        np.std([item['beta'] for item in acceptable_results]),
                        np.std([item['gamma'] for item in acceptable_results]))

      info_line = 'STD:    {:^9.2f}{:^8}{:^55}{:^12.2}{:^12.4f}{:^12.2f}'\
                    ''.format(std_res, '', std_cell, std_spots, std_mos, std_mq)
      ps_log_output.append(info_line)


      # Select the 25% with lowest mosaicities, then select for most spots
      sorted_entries = sorted(acceptable_results, key=lambda i: i['mos'])
      subset = [j[1] for j in enumerate(sorted_entries) \
                if j[0] <= len(sorted_entries) * 0.25]
      sub_spots = [sp['strong'] for sp in subset]

      best = subset[np.argmax(sub_spots)]

      with open('{}/selected.lst'.format(os.path.abspath(gs_params.output)), \
                                                            'a') as sel_int:
        sel_int.write('{}\n'.format(input_file))

      selection_result = [input_file, best['sih'], best['sph'], best['spa'],
                          abs_tmp_dir]

      # Output selected file information
      ps_log_output.append('\nSelected:')
      cell = '{:>8.2f}, {:>8.2f}, {:>8.2f}, {:>6.2f}, {:>6.2f}, {:>6.2f}'\
             ''.format(best['a'], best['b'], best['c'],
                       best['alpha'], best['beta'], best['gamma'])
      info_line = '{:^4}{:^4}{:^9.2f}{:^8}{:^55}{:^12}{:^12.4f}{:^12.2f}'\
                  ''.format(best['sph'], best['spa'], best['res'], best['sg'],
                            cell, best['strong'], best['mos'], best['mq'])
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
  if prog_count < n_int:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step)
  else:
    gs_prog.finished()

  return selection_result

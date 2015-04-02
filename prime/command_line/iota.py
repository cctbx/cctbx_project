from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/01/2015
Description : IOTA command-line module. Version 1.21
'''

import os
import sys
from datetime import datetime
import numpy as np
from collections import Counter

from libtbx.easy_mp import parallel_map
import dials.util.command_line as cmd

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
import prime.iota.iota_select as ps

# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):
  return gs.gs_integration(current_img, len(mp_input_list), log_dir, gs_params)

# Multiprocessor wrapper for selection module
def selection_mproc_wrapper(output_entry):
  return ps.best_file_selection(gs_params, output_entry, log_dir, len(mp_output_list))

# Multiprocessor wrapper for final integration module
def final_mproc_wrapper(current_img):
  return gs.final_integration(current_img, len(sel_clean), log_dir, gs_params)


# ============================================================================ #

def run_grid_search(txt_out, gs_params, gs_range, input_dir_list,
                    input_list, mp_input_list):
  """ Runs grid search in multiprocessing mode.

      input: txt_out - text of *.param file, preserved for analysis
             gs_params - parameters from *.param file in PHIL format
             input_dir_list - list of input folder(s)
             input_list - list of input files
             mp_input_list - as above, but for multiprocessing
  """

 # Log starting info
  inp.main_log(logfile, '{:-^100} \n'.format(' GRID SEARCH AND PICKLE SELECTION '))

  inp.main_log(logfile, '\nSettings for this run:\n')
  inp.main_log(logfile, txt_out)

  with open(gs_params.target, 'r') as phil_file:
    phil_file_contents = phil_file.read()
  inp.main_log(logfile, "\nTarget file ({0}) contents:\n".format(gs_params.target))
  inp.main_log(logfile, phil_file_contents)

  inp.main_log(logfile, 'Found image files in the following folder(s):')
  for folder in input_dir_list:
    inp.main_log(logfile, str(os.path.abspath(folder)))
  inp.main_log(logfile, '\nSpot-finding parameter grid search: ' \
                '{0} input files, spot height: {1} - {2}, '\
                'spot area: {3} - {4} \n'.format(len(input_list),
                 gs_range[0], gs_range[1], gs_range[2], gs_range[3]))
  inp.main_log(logfile, '{:-^100} \n\n'.format(' STARTING GRID SEARCH '))

  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  # run grid search on multiple processes
  cmd.Command.start("Starting Grid Search")
  parallel_map(iterable=mp_input_list,
               func=index_mproc_wrapper,
               processes=gs_params.n_processors)
  cmd.Command.end("Finished Grid Search")

def run_pickle_selection(gs_params, mp_output_list):
  """ Runs pickle_selection in multiprocessing mode.

      input: gs_params - parameters from *.param file in PHIL format
             mp_output_list - list of output folders, for multiprocessing
  """
  inp.main_log(logfile, '\n\n{:-^100} \n'.format(' PICKLE SELECTION '))

  # clear list files from previous selection run
  if os.path.isfile("{}/not_integrated.lst".format(gs_params.output)):
    os.remove("{}/not_integrated.lst".format(gs_params.output))
  if os.path.isfile("{}/prefilter_fail.lst".format(gs_params.output)):
    os.remove("{}/prefilter_fail.lst".format(gs_params.output))
  if os.path.isfile("{}/selected.lst".format(gs_params.output)):
    os.remove("{}/selected.lst".format(gs_params.output))
  if os.path.isfile("{}/integrated.lst".format(gs_params.output)):
    os.remove("{}/integrated.lst".format(gs_params.output))
  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  if not gs_params.grid_search.flag_on:
    inp.main_log(logfile, '\nSettings for this run:\n')
    inp.main_log(logfile, txt_out)

  for output_dir in output_dir_list:
    inp.main_log(logfile, 'Found integrated pickles ' \
                    'under {0}'.format(os.path.abspath(output_dir)))

  if gs_params.flag_prefilter == True:
    prefilter = "ON"
  else:
    prefilter = "OFF"

  inp.main_log(logfile, 'Space group / unit cell prefilter '\
                 'turned {0} \n\n'.format(prefilter))
  inp.main_log(logfile, '{:-^100} \n'.format(' STARTING SELECTION '))

  # run pickle selection on multiple processes
  cmd.Command.start("Starting Pickle Selection")
  selection_results = parallel_map(iterable=mp_output_list,
                                   func=selection_mproc_wrapper,
                                   processes=gs_params.n_processors)
  cmd.Command.end("Finished Pickle Selection")

  return selection_results

def final_integration(sel_clean, gs_params):

  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  inp.main_log(logfile, "\n\n{:-^80}\n".format(' FINAL INTEGRATION '))

  cmd.Command.start("Integrating with selected spotfinding parameters")
  result_objects = parallel_map(iterable=sel_clean,
                                func=final_mproc_wrapper,
                                processes=gs_params.n_processors,
                                preserve_exception_message=True)
  cmd.Command.end("Integrating with selected spotfinding parameters -- DONE ")

  clean_results = [results for results in result_objects if results != []]

  return clean_results

def print_results(clean_results, gs_range):

  images = [results['img'] for results in clean_results]
  spot_heights = [results['sph'] for results in clean_results]
  spot_areas = [results['spa'] for results in clean_results]
  resolutions = [results['res'] for results in clean_results]
  num_spots = [results['strong'] for results in clean_results]
  mosaicities = [results['mos'] for results in clean_results]

  a = [results['a'] for results in clean_results]
  b = [results['b'] for results in clean_results]
  c = [results['c'] for results in clean_results]
  alpha = [results['alpha'] for results in clean_results]
  beta = [results['beta'] for results in clean_results]
  gamma = [results['gamma'] for results in clean_results]

  cons_h = Counter(spot_heights).most_common(1)[0][0]
  cons_a = Counter(spot_areas).most_common(1)[0][0]

  final_table = []
  final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))
  final_table.append("Total images:          {}".format(len(images)))
  final_table.append("Avg. spot height:      {:<8.3f}  std. dev:    {:<6.2f}"\
                     "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                     "".format(np.mean(spot_heights), np.std(spot_heights),
                               max(spot_heights), min(spot_heights), cons_h))
  final_table.append("Avg. spot areas:       {:<8.3f}  std. dev:    {:<6.2f}"\
                    "  max: {:<3}  min: {:<3}  consensus: {:<3}"\
                    "".format(np.mean(spot_areas), np.std(spot_areas),
                              max(spot_areas), min(spot_areas), cons_a))
  final_table.append("Avg. number of spots:  {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(num_spots), np.std(num_spots)))

  final_table.append("Avg. mosaicity:        {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(mosaicities), np.std(mosaicities)))
  final_table.append("Avg. unit cell:        "\
                     "{:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f}), "\
                     "{:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f}), "\
                     "{:<6.2f} ({:>4.2f}), {:<6.2f} ({:>4.2f})"
                    "".format(np.mean(a), np.std(a),
                              np.mean(b), np.std(b),
                              np.mean(c), np.std(c),
                              np.mean(alpha), np.std(alpha),
                              np.mean(beta), np.std(beta),
                              np.mean(gamma), np.std(gamma)))

  bad_mos_list = [item for item in clean_results if float(item['mos']) -\
                   np.mean(mosaicities) > np.std(mosaicities) * 2]

  for item in final_table:
      print item
      inp.main_log(logfile, item)

def print_summary(gs_params):
  """ Prints summary by reading contents of files listing
      a) images not integrated
      b) images that failed unit cell filter
      c) total images input
      d) final images successfully processed

      Appends summary to general log file. Also outputs some of it on stdout.

      input: gs_params - parameters from *.param file in PHIL format
  """

  summary = []
  int_fail_count = 0
  bad_int_count = 0
  final_count = 0

  print "\n\n{:-^80}\n".format('SUMMARY')
  inp.main_log(logfile, "\n\n{:-^80}\n".format('SUMMARY'))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as plog:
    prog_content = plog.read()

  if gs_params.advanced.random_sample.flag_on == True:
    summary.append('raw images processed:         {}'\
                   ''.format(gs_params.advanced.random_sample.number))
  else:
    summary.append('raw images processed:         {}'.format(len(input_list)))

  if os.path.isfile('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/not_integrated.lst'.format(os.path.abspath(gs_params.output)),
              'r') as int_fail_list:
      int_fail_list_contents = int_fail_list.read()
      int_fail_count = len(int_fail_list_contents.splitlines())

  if os.path.isfile('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/prefilter_fail.lst'.format(os.path.abspath(gs_params.output)),
            'r') as bad_int_list:
      bad_int_list_contents = bad_int_list.read()
      bad_int_count = len(bad_int_list_contents.splitlines())

  if os.path.isfile('{0}/selected.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/selected.lst'.format(os.path.abspath(gs_params.output)),
              'r') as sel_list:
      sel_list_contents = sel_list.read()
      sel_count = len(sel_list_contents.splitlines())

  if os.path.isfile('{0}/integrated.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/integrated.lst'.format(os.path.abspath(gs_params.output)),
              'r') as final_list:
      final_list_contents = final_list.read()
      final_count = len(final_list_contents.splitlines())

  summary.append('raw images not integrated:    {}'.format(int_fail_count))
  summary.append('images failed prefilter:      {}'.format(bad_int_count))
  summary.append('images in selection:          {}'.format(sel_count))
  summary.append('final integrated pickles:     {}'.format(final_count))
  summary.append('\n\nIOTA version {0}'.format(iota_version))

  for item in summary:
    print item
    inp.main_log(logfile, "{}".format(item))

  inp.main_log(logfile, "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now()))

# ============================================================================ #

if __name__ == "__main__":

  iota_version = '1.21'

  print "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
  print 'Starting IOTA ... \n\n'

  # read parameters from *.param file
  gs_params, txt_out = inp.process_input(sys.argv[1:])

  # If no parameter file specified, output blank parameter file
  if gs_params.input == None:
    print '{:-^100}\n\n'.format('IOTA Dry Run')
    print txt_out
    with open('iota.param', 'w') as default_settings_file:
      default_settings_file.write(txt_out)
    sys.exit()

  # generate input
  gs_range, input_list, input_dir_list, output_dir_list, log_dir, logfile,\
  mp_input_list, mp_output_list = inp.generate_input(gs_params)

  # debugging/experimental section - anything goes here
  if gs_params.advanced.debug:
    sys.exit()

  # run grid search
  if gs_params.grid_search.flag_on:
    run_grid_search(txt_out, gs_params, gs_range, input_dir_list,
                    input_list, mp_input_list)

  # run pickle selection
  selection_results = run_pickle_selection(gs_params, mp_output_list)
  sel_clean = [entry for entry in selection_results \
                     if entry != None and entry != []]

  # run final integration
  final_int = final_integration(sel_clean, gs_params)

  # print final integration results and summary
  print_results(final_int, gs_range)
  print_summary(gs_params)

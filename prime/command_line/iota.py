from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 03/06/2015
Description : IOTA command-line module. Version 0.9
'''

import os
import sys
import shutil
import random
from datetime import datetime
import numpy as np

from libtbx.easy_mp import parallel_map
import dials.util.command_line as cmd

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
import prime.iota.iota_select as ps
from prime.iota.iota_select import best_file_selection


# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):
  return gs.integrate_one_image(current_img, len(mp_input_list), log_dir,
                                False, gs_params)

# Multiprocessor wrapper for selection module
def selection_mproc_wrapper(output_entry):
  return ps.best_file_selection(gs_params, output_entry, log_dir, len(mp_output_list))

# Multiprocessor wrapper for final integration module
def final_mproc_wrapper(current_img):
  return gs.exp_integrate_one(current_img, log_dir, len(sel_clean), gs_params)

# Multiprocessor wrapper for single image integration module
def single_mproc_wrapper(current_img):
  return gs.int_single_image(current_img, log_dir, n_int, gs_params)

def experimental_mproc_wrapper(single_entry):
  return gs.integrate_selected_image(single_entry, log_dir, gs_params)

def generate_input(gs_params):
  """ This section generates input for grid search and/or pickle selection.

      parameters: gs_params - list of parameters from *.param file (in
      PHIL format)

      output: input_list - list of absolute paths to input files
              input_dir_list - list of absolute paths to input folder(s)
              output_dir_list - same for output folder(s)
              mp_input_list - for multiprocessing: filename + spotfinding paarams
              mp_output_list - same for output
  """

  # If no parameter file specified, output blank parameter file
  if gs_params.input == None:
    print '{:-^100}\n\n'.format('IOTA Dry Run')
    print txt_out
    with open('iota_default.phil', 'w') as default_settings_file:
      default_settings_file.write(txt_out)
    sys.exit()

  # Check for list of files; if exists, use as input list. If doesn't exist,
  # generate list from input directory
  if gs_params.input_list == None:
    input_list = inp.make_input_list(gs_params)
  else:
    with open(gs_params.input_list, 'r') as listfile:
      listfile_contents = listfile.read()
    input_list = listfile_contents.splitlines()

  # If grid-search turned on, check for existing output directory and remove
  if gs_params.grid_search.flag_on == True:

    # Remove old output if found
    if os.path.exists(os.path.abspath(gs_params.output)):
      cmd.Command.start("Found old folder {}... deleting...".format(gs_params.output))
      shutil.rmtree(os.path.abspath(gs_params.output))
      cmd.Command.end("Deleted old folder {} -- DONE".format(gs_params.output))

    # Make main output directory and log directory
    os.makedirs(os.path.abspath(gs_params.output))
    os.makedirs("{}/logs".format(os.path.abspath(gs_params.output)))

    # Make per-image output folders. ***May not be necessary!!
    cmd.Command.start("Generating output directory tree")
    inp.make_dirs(input_list, gs_params)
    cmd.Command.end("Generating output directory tree -- DONE")

  # If grid-search turned off, check that output directory exists, so that the
  # selection module has input
  else:
    if not os.path.exists(os.path.abspath(gs_params.output)):
      print "ERROR: No grid search results detected in"\
          "{}".format(os.path.abspath(gs_params.output))
      sys.exit()

  # Make log directory and input/output directory lists
  log_dir = "{}/logs".format(os.path.abspath(gs_params.output))
  input_dir_list, output_dir_list = inp.make_dir_lists(input_list, gs_params)

  # Initiate log file
  logfile = '{}/iota.log'.format(log_dir)
  inp.main_log_init(logfile)

  # Make input/output lists for multiprocessing
  cmd.Command.start("Generating multiprocessing input")
  mp_input_list, mp_output_list = inp.make_mp_input(input_list, gs_params)
  cmd.Command.end("Generating multiprocessing input -- DONE")

  return input_list, input_dir_list, output_dir_list, log_dir, logfile,\
         mp_input_list, mp_output_list

# =========================== EXPERIMENTAL SECTION =========================== #

def advanced_input(gs_params):
  """This is for various debugging / experimental stuff. Runs one image only.
  """
  print "EXPERIMENTAL / DEBUG MODE:"
  print "Trying out extracting relevant parameters from an image..."

  trial_list = []
  if gs_params.advanced.single_img:
    random_number = random.randrange(0, len(input_list))
    trial_list.append(input_list[random_number])
  elif gs_params.random_sample.flag_on:
    for i in range (0, gs_params.random_sample.number + 1):
      random_number = random.randrange(0, len(input_list))
      trial_list.append(input_list[random_number])
  else:
    trial_list = input_list

  mp_trial_list, mp_trial_output_list = inp.make_mp_input(trial_list, gs_params)

  #gs.debug_integrate_one(mp_trial_list[0], log_dir, gs_params)
  #gs.integrate_one_image(mp_trial_list[0], len(mp_input_list), log_dir, gs_params)

  # run grid search on multiple processes
  parallel_map(iterable=mp_trial_list,
               func=experimental_mproc_wrapper,
               processes=gs_params.n_processors,
               preserve_exception_message=False)

  print "END OF RUN"
# ============================================================================ #

def run_grid_search(txt_out, gs_params, input_dir_list, input_list, mp_input_list):
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
                gs_params.grid_search.h_min, gs_params.grid_search.h_max,
                gs_params.grid_search.a_min, gs_params.grid_search.a_max))
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
  if os.path.isfile("{}/selected.lst".format(gs_params.output)):
    os.remove("{}/prefilter_fail.lst".format(gs_params.output))
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
                                   processes=gs_params.n_processors,
                                   preserve_exception_message=True)
  cmd.Command.end("Finished Pickle Selection")

  return selection_results

def final_integration(sel_clean, gs_params):

  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  if gs_params.advanced.charts:
    inp.main_log(logfile, "\n\n{:-^80}\n".format('FINAL INTEGRATION WITH ALL PLOT PDFS'))
  elif gs_params.advanced.mosaicity_plot:
    inp.main_log(logfile, "\n\n{:-^80}\n".format('FINAL INTEGRATION WITH MOSAICITY PLOT PDFS'))
  else:
    inp.main_log(logfile, "\n\n{:-^80}\n".format('FINAL INTEGRATION, NO PLOTS'))

  cmd.Command.start("Integrating with selected spotfinding parameters")
  result_objects = parallel_map(iterable=sel_clean,
               func=final_mproc_wrapper,
               processes=gs_params.n_processors,
               preserve_exception_message=True)
  cmd.Command.end("Integrating with selected spotfinding parameters -- DONE ")

  clean_results = [results for results in result_objects if results != []]

  return clean_results

def print_results(clean_results):

  images = [results[0] for results in clean_results]
  spot_heights = [int(results[1]) for results in clean_results]
  spot_areas = [int(results[2]) for results in clean_results]
  resolutions = [float(results[3]) for results in clean_results]
  num_spots = [int(results[4]) for results in clean_results]
  dom_sizes = [float(results[5]) for results in clean_results]
  mosaicities = [float(results[6]) for results in clean_results]
  rmsds = [float(results[7]) for results in clean_results]

  a = [results[8][0] for results in clean_results]
  b = [results[8][1] for results in clean_results]
  c = [results[8][2] for results in clean_results]
  alpha = [results[8][3] for results in clean_results]
  beta = [results[8][4] for results in clean_results]
  gamma = [results[8][5] for results in clean_results]


  final_table = []
  final_table.append("\n\n{:-^80}\n".format('ANALYSIS OF RESULTS'))
  final_table.append("Total images:          {}".format(len(images)))
  final_table.append("Avg. spot height:      {:<8.3f}  std. dev:    {:<6.2f}"\
                     "  max: {:<3}  min: {:<3}".format(np.mean(spot_heights),
                     np.std(spot_heights), max(spot_heights), min(spot_heights)))
  final_table.append("Avg. spot areas:       {:<8.3f}  std. dev:    {:<6.2f}"\
                    "  max: {:<3}  min: {:<3}".format(np.mean(spot_areas),
                    np.std(spot_areas), max(spot_areas), min(spot_areas)))
  final_table.append("Avg. number of spots:  {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(num_spots), np.std(num_spots)))
  final_table.append("Avg. domain size:      {:<8.3f}  std. dev:    {:<6.2f}"\
                     "".format(np.mean(dom_sizes), np.std(dom_sizes)))
  final_table.append("Avg. mosaicity:        {:<8.3f}  std. dev:    {:<6.2f}"\
                    "".format(np.mean(mosaicities), np.std(mosaicities)))
  final_table.append("Avg. positional RMSD:  {:<8.3f}  std. dev:    {:<6.2f}"
                    "".format(np.mean(rmsds), np.std(rmsds)))
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

  bad_mos_list = [item for item in clean_results if float(item[6]) -\
                   np.mean(mosaicities) > np.std(mosaicities) * 2]

  if len(bad_mos_list) != 0:
    final_table.append("\nImages with poor mosaicity (> 2-sigma over mean):")
    for entry in bad_mos_list:
      final_table.append('{:<3}  {}: H = {:<3} A = {:<3}  RES = {:<4.2F}  SPOTS = {:<6}'\
                         ' DOM = {:<8.2f}  MOS = {:<6.4f}  RMSD = {:<6.2f}'\
                         ''.format(bad_mos_list.index(entry), entry[0],
                                   entry[1], entry[2], entry[3], entry[4],
                                   entry[5], entry[6], entry[7]))

  bad_rmsd_list = [item for item in clean_results if float(item[7]) - \
                     np.mean(rmsds) > np.std(rmsds) * 2]

  if len(bad_rmsd_list) != 0:
    i = 0
    final_table.append("\nImages with poor positional rmsd (> 2-sigma over mean):")
    for entry in bad_rmsd_list:
      final_table.append('{:<3}  {}: H = {:<3} A = {:<3}  RES = {:<4.2F}  SPOTS = {:<6} '\
                         'DOM = {:<8.2f}  MOS = {:<6.4f}  RMSD = {:<6.2f}'\
                         ''.format(bad_rmsd_list.index(entry), entry[0],
                                   entry[1], entry[2], entry[3], entry[4],
                                   entry[5], entry[6], entry[7]))

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

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()

  if gs_params.advanced.random_sample.flag_on == True:
    summary.append('raw images processed:         {}'.format(gs_params_random_sample.number))
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


def single_image_mode(gs_params):

  current_img = gs_params.advanced.single_img
  img_filename = os.path.basename(current_img).split('.')[0]

  # Remove old output if found
  if os.path.exists(os.path.abspath(gs_params.output)):
    cmd.Command.start("Found old folder {}... deleting...".format(gs_params.output))
    shutil.rmtree(os.path.abspath(gs_params.output))
    cmd.Command.end("Deleted old folder {} -- DONE".format(gs_params.output))

  # Make main output directory and log directory
  os.makedirs(os.path.abspath(gs_params.output))
  os.makedirs("{}/logs".format(os.path.abspath(gs_params.output)))

  output_dir = os.path.abspath(gs_params.output)
  tmp_output_dir = os.path.abspath("{0}/tmp_{1}"\
                                   "".format(output_dir, img_filename))
  if not os.path.exists(tmp_output_dir):
    os.makedirs(tmp_output_dir)
  print tmp_output_dir

  single_mp_list = []
  for sig_height in range(gs_params.grid_search.h_min,
                       gs_params.grid_search.h_max + 1):
   for spot_area in range (gs_params.grid_search.a_min,
                           gs_params.grid_search.a_max + 1):
     mp_entry = [current_img, sig_height, sig_height, spot_area]
     single_mp_list.append(mp_entry)

  cmd.Command.start("Processing single image ")
  int_results = parallel_map(iterable=single_mp_list,
                                func=single_mproc_wrapper,
                                processes=gs_params.n_processors)
  cmd.Command.end("Processing single image -- DONE ")

  single_img_results = [entry for entry in int_results if entry != []]

  return single_img_results


# ============================================================================ #

if __name__ == "__main__":

  iota_version = '1.03'

  print "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
  print 'Starting IOTA ... \n\n'

  # read parameters from *.param file
  gs_params, txt_out = inp.process_input(sys.argv[1:])

  if gs_params.advanced.single_img != None:
    log_dir = '{}/logs'.format(gs_params.output)
    logfile = '{}/iota.log'.format(log_dir)
    n_int = (gs_params.grid_search.a_max - gs_params.grid_search.a_min + 1) * \
            (gs_params.grid_search.h_max - gs_params.grid_search.h_min + 1)
    single_img_results = single_image_mode(gs_params)

    print_results(single_img_results)
    sys.exit()

  # generate input
  input_list, input_dir_list, output_dir_list, log_dir, logfile, \
  mp_input_list, mp_output_list = generate_input(gs_params)

  # debugging/experimental section - anything goes here
  if gs_params.advanced.debug:
    gs.integrate_one_image(mp_input_list[0], len(mp_input_list), log_dir,
                                False, gs_params)
    print "debugging... {}".format(mp_input_list[0])
    sys.exit()

  #run grid search
  if gs_params.grid_search.flag_on:
    run_grid_search(txt_out, gs_params, input_dir_list, input_list, mp_input_list)

  #run pickle selection
  selection_results = run_pickle_selection(gs_params, mp_output_list)
  sel_clean = [entry for entry in selection_results \
                     if entry != None and entry != []]

  final_int = final_integration(sel_clean, gs_params)       # final integration
  print_results(final_int)
  print_summary(gs_params)                                  # print summary

from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 02/19/2015
Description : IOTA command-line module. Version 0.9
'''

import os
import sys
import shutil
import random
import logging
from datetime import datetime

from libtbx.easy_mp import parallel_map
import dials.util.command_line as cmd

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
import prime.iota.iota_select as ps
from prime.iota.iota_select import best_file_selection


# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):
  return gs.integrate_one_image(current_img, len(mp_input_list), log_dir, gs_params)

# Multiprocessor wrapper for selection module
def selection_mproc_wrapper(output_entry):
  return ps.best_file_selection(gs_params, output_entry, log_dir, len(mp_output_list))

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

  # Make input/output lists for multiprocessing
  cmd.Command.start("Generating multiprocessing input")
  mp_input_list, mp_output_list = inp.make_mp_input(input_list, gs_params)
  cmd.Command.end("Generating multiprocessing input -- DONE")

  return input_list, input_dir_list, output_dir_list, log_dir,\
         mp_input_list, mp_output_list

def advanced_input(gs_params):
  """This is for various debugging / experimental stuff. Runs one image only.
  """
  print "DEBUG MODE:"
  random_number = random.randrange(0, len(mp_input_list))
  single_entry = mp_input_list[random_number]
  #gs.debug_integrate_one(single_entry, log_dir, gs_params)
  gs.integrate_one_image(single_entry, len(mp_input_list), log_dir, gs_params)
  print "END OF RUN"
  sys.exit()

def run_grid_search(txt_out, gs_params, input_dir_list, input_list, mp_input_list):
  """ Runs grid search in multiprocessing mode.

      input: txt_out - text of *.param file, preserved for analysis
             gs_params - parameters from *.param file in PHIL format
             input_dir_list - list of input folder(s)
             input_list - list of input files
             mp_input_list - as above, but for multiprocessing
  """

 # Log starting info
  logging.info('{:-^100} \n'.format(' GRID SEARCH AND PICKLE SELECTION '))

  logging.info('\nSettings for this run:\n')
  logging.info(txt_out)

  with open(gs_params.target, 'r') as phil_file:
    phil_file_contents = phil_file.read()
  logging.info("\nTarget file ({0}) contents:\n".format(gs_params.target))
  logging.info(phil_file_contents)

  logging.info('Found image files in the following folder(s):')
  for folder in input_dir_list:
    logging.info(str(os.path.abspath(folder)))
  logging.info('\nSpot-finding parameter grid search: ' \
                '{0} input files, spot height: {1} - {2}, '\
                'spot area: {3} - {4} \n'.format(len(input_list),
                gs_params.grid_search.h_min, gs_params.grid_search.h_max,
                gs_params.grid_search.a_min, gs_params.grid_search.a_max))
  logging.info('{:-^100} \n\n'.format(' STARTING GRID SEARCH '))


  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  # run grid search on multiple processes
  cmd.Command.start("Starting Grid Search")
  parallel_map(iterable=mp_input_list,
               func=index_mproc_wrapper,
               processes=gs_params.n_processors,
               preserve_exception_message=False)
  cmd.Command.end("Finished Grid Search")

def run_pickle_selection(gs_params, mp_output_list):
  """ Runs pickle_selection in multiprocessing mode.

      input: gs_params - parameters from *.param file in PHIL format
             mp_output_list - list of output folders, for multiprocessing
  """
  logging.info('\n\n{:-^100} \n'.format(' PICKLE SELECTION '))

  # clear list files from previous selection run
  if os.path.isfile("{}/best_by_strong.lst".format(gs_params.output)):
    os.remove("{}/best_by_strong.lst".format(gs_params.output))
  if os.path.isfile("{}/best_by_offset.lst".format(gs_params.output)):
    os.remove("{}/best_by_offset.lst".format(gs_params.output))
  if os.path.isfile("{}/best_by_uc.lst".format(gs_params.output)):
    os.remove("{}/best_by_uc.lst".format(gs_params.output))
  if os.path.isfile("{}/best_by_total.lst".format(gs_params.output)):
    os.remove("{}/best_by_total.lst".format(gs_params.output))
  if os.path.isfile("{}/not_integrated.lst".format(gs_params.output)):
    os.remove("{}/not_integrated.lst".format(gs_params.output))
  if os.path.isfile("{}/prefilter_fail.lst".format(gs_params.output)):
    os.remove("{}/prefilter_fail.lst".format(gs_params.output))
  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  if not gs_params.grid_search.flag_on:
    logging.info('\nSettings for this run:\n')
    logging.info(txt_out)

  for output_dir in output_dir_list:
    logging.info('Found integrated pickles ' \
                    'under {0}'.format(os.path.abspath(output_dir)))

  if gs_params.flag_prefilter == True:
    prefilter = "ON"
  else:
    prefilter = "OFF"

  logging.info('Space group / unit cell prefilter '\
                 'turned {0} \n\n'.format(prefilter))
  logging.info('{:-^100} \n'.format(' STARTING SELECTION '))

  #for output_entry in mp_output_list:
  #     best_file_selection(gs_params, output_entry, log_dir)

  # run pickle selection on multiple processes
  cmd.Command.start("Starting Pickle Selection")
  parallel_map(iterable=mp_output_list,
               func=selection_mproc_wrapper,
               processes=gs_params.n_processors)
  cmd.Command.end("Finished Pickle Selection")

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
  logging.info("\n\n{:-^80}\n".format('SUMMARY'))

  with (open ('{0}/logs/progress.log'.format(gs_params.output), 'r')) as prog_log:
    prog_content = prog_log.read()

  for item in prog_content.splitlines():
    print item
  print '\n\n'

  if gs_params.random_sample.flag_on == True:
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

  if os.path.isfile('{0}/best_by_strong.lst'.format(os.path.abspath(gs_params.output))):
    with open('{0}/best_by_strong.lst'.format(os.path.abspath(gs_params.output)),
              'r') as output_list:
      output_list_contents = output_list.read()
      final_count = len(output_list_contents.splitlines())

  summary.append('raw images not integrated:    {}'.format(int_fail_count))
  summary.append('images failed prefilter:      {}'.format(bad_int_count))
  summary.append('pickles in final selection:   {}'.format(final_count))
  summary.append('\n\nIOTA version {0}'.format(iota_version))

  for item in summary:
    print item
    logging.info("{}".format(item))

  logging.info("{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now()))


# ============================================================================ #

if __name__ == "__main__":

  iota_version = '0.92'

  print "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
  print 'Starting IOTA ... \n\n'

  # read parameters from *.param file
  gs_params, txt_out = inp.process_input(sys.argv[1:])

  # generate input
  input_list, input_dir_list, output_dir_list, log_dir, mp_input_list, \
  mp_output_list = generate_input(gs_params)

  # debugging/experimental section - anything goes here
  if gs_params.advanced.single_img == True:
    advanced_input()

  # Setup grid search logger
  logging.basicConfig(level=logging.INFO,
                      format='%(message)s',
                      filename='{0}/iota.log'.format(log_dir),
                      filemode='w')

  #run grid search
  if gs_params.grid_search.flag_on == True:
    run_grid_search(txt_out, gs_params, input_dir_list, input_list, mp_input_list)

  #run pickle selection
  run_pickle_selection(gs_params, mp_output_list)

  # print summary
  print_summary(gs_params)

from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/10/2015
Description : IOTA command-line module. Version 1.24
'''

import os
import sys
from datetime import datetime

from libtbx.easy_mp import parallel_map
import dials.util.command_line as cmd

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
import prime.iota.iota_select as ps
import prime.iota.iota_analysis as ia

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
  cmd.Command.start("Spotfinding Grid Search")
  parallel_map(iterable=mp_input_list,
               func=index_mproc_wrapper,
               processes=gs_params.n_processors)
  cmd.Command.end("Spotfinding Grid Search -- DONE")

def run_pickle_selection(gs_params, mp_output_list):
  """ Runs pickle_selection in multiprocessing mode.

      input: gs_params - parameters from *.param file in PHIL format
             mp_output_list - list of output folders, for multiprocessing

      output: selection_results - list of images w/ spotfinding params
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
  cmd.Command.start("Spotfinding Combination Selection")
  selection_results = parallel_map(iterable=mp_output_list,
                                   func=selection_mproc_wrapper,
                                   processes=gs_params.n_processors)
  cmd.Command.end("Spotfinding Combination Selection -- DONE")

  return selection_results

def final_integration(sel_clean, gs_params):
  """ Runs integration in multiprocessing mode. Uses spotfinding parameters
      determined by grid search and selection.

      input: gs_params - parameters from *.param file in PHIL format
             sel_clean - list of images w/ optimal spotfinding params

      output: clean_results - list of integrated pickles w/ integration data
  """

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

def dry_run():

    print '\n{:-^70}\n'.format('IOTA Usage')

    print 'Dry run mode'
    print 'Usage: prime.iota'
    print 'Generates default parameter file for IOTA (iota.param) and '
    print 'target file for cctbx.xfel (target.phil).\n'

    print 'Auto mode'
    print 'Usage: prime.iota path/to/raw/images'
    print 'Generates two files, parameter file for IOTA (iota.param) and'
    print 'target file for cctbx.xfel (target.phil). Integrates a random'
    print 'subset of images without target cell. Outputs basic analysis.'
    print 'The input files currently MUST be in pickle format. Modify by'
    print 'running cxi.image2pickle.\n'

    print 'Script mode'
    print 'Usage: prime.iota <script>.param'
    print 'Run using IOTA parameter file and target PHIL file generated from'
    print 'the dry run or auto mode. Make sure that IOTA parameter file has'
    print 'the path to the input image folder under "input". The input files'
    print 'currently MUST be in pickle format. Modify by running'
    print 'cxi.image2pickle.\n\n'

    help_out, txt_out = inp.print_params()

    print '\n{:-^70}\n'.format('IOTA Parameters')
    print help_out

    inp.write_defaults(os.path.abspath(os.path.curdir), txt_out)

# ============================================================================ #

if __name__ == "__main__":

  iota_version = '1.24'

  print "\n\n"
  print "     IIIIII          OOOO         TTTTTTTTTT           A              "
  print "       II           O    O            TT              A A             "
  print "       II           O    O            TT             A   A            "
  print ">------INTEGRATION--OPTIMIZATION------TRIAGE--------ANALYSIS--------->"
  print "       II           O    O            TT           A       A          "
  print "       II           O    O            TT          A         A         "
  print "     IIIIII          OOOO             TT         A           A   v{}"\
        "".format(iota_version)

  # read parameters from *.param file
  arg = sys.argv[1:]
  now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())

  # If no parameter file specified, output blank parameter file
  if arg == []:
    dry_run()
    sys.exit()
  else:
    print '\n{}\n'.format(now)
    carg = arg[0]
    if os.path.exists(os.path.abspath(carg)):

      # If user provided a parameter file
      if os.path.isfile(carg) and carg.endswith('.param'):
        gs_params, txt_out = inp.process_input(arg)

      # If user provided a data folder
      elif os.path.isdir(carg):
        print "IOTA will run in AUTO mode:\n"
        cmd.Command.start("Generating default parameters")
        gs_params, txt_out = inp.auto_mode(os.path.abspath(os.path.curdir),
                                           os.path.abspath(carg), now)
        cmd.Command.end("Generating default parameters -- DONE")

    # If user provided gibberish
    else:
      print "ERROR: Incorrect input! Need parameter filename or data folder."
      sys.exit()



  # generate input
  gs_range, input_list, input_dir_list, output_dir_list, log_dir, logfile,\
  mp_input_list, mp_output_list = inp.generate_input(gs_params)

  # debugging/experimental section - anything goes here
  #if gs_params.advanced.debug:
  #  sys.exit()

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
  ia.print_results(final_int, gs_range, logfile)
  ia.print_summary(gs_params, len(input_list), logfile, iota_version, now)

from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 05/05/2015
Description : IOTA command-line module. Version 1.31
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
  return gs.integration("grid", current_img, len(mp_input_list),
                        log_dir, gs_params)

# Multiprocessor wrapper for selection module after grid search
def sel_grid_mproc_wrapper(output_entry):
  return ps.best_file_selection("grid", gs_params, output_entry, log_dir,
                                len(mp_output_list))

# Multiprocessor wrapper for final integration module
def final_mproc_wrapper(current_img):
  return gs.integration("final", current_img, len(sel_clean),
                        log_dir, gs_params)


# ============================================================================ #

def iota_start(txt_out, gs_params, gs_range, input_dir_list,
                    input_list, mp_input_list):

  # Log starting info
  inp.main_log(logfile, '{:=^100} \n'.format(' IOTA MAIN LOG '))

  inp.main_log(logfile, '{:-^100} \n'.format(' SETTINGS FOR THIS RUN '))
  inp.main_log(logfile, txt_out)

  if gs_params.grid_search.flag_on:
    inp.main_log(logfile, '{:-^100} \n\n'.format(' TARGET FILE ({}) CONTENTS '
                                       ''.format(gs_params.target)))
    with open(gs_params.target, 'r') as phil_file:
      phil_file_contents = phil_file.read()
    inp.main_log(logfile, phil_file_contents)

    inp.main_log(logfile, '{:-^100} \n\n'.format(''))
    inp.main_log(logfile, 'Found image files in the following folder(s):')
    for folder in input_dir_list:
      inp.main_log(logfile, str(os.path.abspath(folder)))
    inp.main_log(logfile, '\nSpot-finding parameter grid search: ' \
                  '{0} input files, spot height: {1} - {2}, '\
                  'spot area: {3} - {4} \n'.format(len(input_list),
                   gs_range[0], gs_range[1], gs_range[2], gs_range[3]))


def run_integration(int_type, gs_params, mp_input_list):
  """ Runs grid search in multiprocessing mode.

      input: txt_out - text of *.param file, preserved for analysis
             gs_params - parameters from *.param file in PHIL format
             input_dir_list - list of input folder(s)
             input_list - list of input files
             mp_input_list - as above, but for multiprocessing
  """

  if os.path.isfile("{0}/logs/progress.log".format(gs_params.output)):
    os.remove("{0}/logs/progress.log".format(gs_params.output))

  if int_type == 'grid':
    inp.main_log(logfile, '{:-^100} \n\n'.format(' SPOTFINDING GRID SEARCH '))

    # run grid search on multiple processes
    cmd.Command.start("Spotfinding Grid Search")
    parallel_map(iterable=mp_input_list,
                 func=index_mproc_wrapper,
                 processes=gs_params.n_processors)
    cmd.Command.end("Spotfinding Grid Search -- DONE")

  elif int_type == 'final':
    inp.main_log(logfile, "\n\n{:-^80}\n".format(' FINAL INTEGRATION '))

    cmd.Command.start("Final integration")
    result_objects = parallel_map(iterable=mp_input_list,
                                  func=final_mproc_wrapper,
                                  processes=gs_params.n_processors,
                                  preserve_exception_message=True)
    cmd.Command.end("Final integration -- DONE ")

    clean_results = [results for results in result_objects if results != []]
    return clean_results

def run_selection(sel_type, gs_params, mp_output_list):
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
  if os.path.isfile("{}/gs_selected.lst".format(gs_params.output)):
    os.remove("{}/gs_selected.lst".format(gs_params.output))
  if os.path.isfile("{}/mos_selected.lst".format(gs_params.output)):
    os.remove("{}/mos_selected.lst".format(gs_params.output))
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
  if sel_type == 'grid':
    cmd.Command.start("Spotfinding Combination Selection")
    selection_results = parallel_map(iterable=mp_output_list,
                                     func=sel_grid_mproc_wrapper,
                                     processes=gs_params.n_processors)
    cmd.Command.end("Spotfinding Combination Selection -- DONE")

  return selection_results


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


def experimental(gs_params, log_dir):
  """EXPERIMENTAL SECTION: Contains stuff I just want to try out without
     running the whole darn thing
  """
  import prime.iota.iota_index as ix

  sample_img = gs_params.advanced.single_img
  print "Have image {}".format(sample_img)

  ix.spotfinding_param_search(sample_img, gs_params)


# ============================================================================ #

if __name__ == "__main__":

  iota_version = '1.31'

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
      print "ERROR: Invalid input! Need parameter filename or data folder."
      sys.exit()



  # generate input
  gs_range, input_list, input_dir_list, output_dir_list, log_dir, logfile,\
  mp_input_list, mp_output_list = inp.generate_input(gs_params)

  iota_start(txt_out, gs_params, gs_range, input_dir_list,
                    input_list, mp_input_list)

  # debugging/experimental section - anything goes here
  if gs_params.advanced.experimental:
    print "IOTA will run in EXPERIMENTAL mode:\n"
    experimental(gs_params, log_dir)
    sys.exit()

  # run grid search
  if gs_params.grid_search.flag_on:
    run_integration("grid", gs_params, mp_input_list)

  # run pickle selection
  selection_results = run_selection('grid', gs_params, mp_output_list)
  sel_clean = [entry for entry in selection_results \
                     if entry != None and entry != []]

  # run final integration
  final_int = run_integration("final", gs_params, sel_clean)

  # print final integration results and summary
  ia.print_results(final_int, gs_range, logfile)
  ia.print_summary(gs_params, len(input_list), logfile, iota_version, now)

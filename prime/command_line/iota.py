from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 05/19/2015
Description : IOTA command-line module. Version 1.50
'''

import os
import sys
import argparse
import shutil
from datetime import datetime

from libtbx.easy_mp import parallel_map

import prime.iota.iota_input as inp
import prime.iota.iota_gridsearch as gs
import prime.iota.iota_select as ps
import prime.iota.iota_analysis as ia
import prime.iota.iota_cmd as cmd
import prime.iota.iota_conversion as i2p

# Multiprocessor wrapper for grid search module
def index_mproc_wrapper(current_img):

  prog_count = mp_input_list.index(current_img)
  n_int = len(mp_input_list)

  gs_prog = cmd.ProgressBar(title='GRID SEARCH')
  if prog_count < n_int:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return gs.integration("grid", current_img, log_dir, gs_params)

# Multiprocessor wrapper for selection module after grid search
def sel_grid_mproc_wrapper(output_entry):

  prog_count = mp_output_list.index(output_entry)
  n_int = len(mp_output_list)

  gs_prog = cmd.ProgressBar(title='SELECTING')
  if prog_count < n_int:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return ps.best_file_selection("grid", gs_params, output_entry, log_dir)

# Multiprocessor wrapper for final integration module
def final_mproc_wrapper(current_img):

  prog_count = sel_clean.index(current_img)
  n_int = len(sel_clean)

  gs_prog = cmd.ProgressBar(title='INTEGRATING')
  if prog_count < n_int:
    prog_step = 100 / n_int
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return gs.integration("final", current_img, log_dir, gs_params)

def conversion_wrapper(img_entry):
  """ Multiprocessor wrapper, for raw image to pickle conversion
  """

  prog_count = raw_input_list.index(img_entry)
  n_img = len(raw_input_list)

  gs_prog = cmd.ProgressBar(title='CONVERTING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return i2p.convert_image(img_entry[0], img_entry[1], square, True)

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

  if int_type == 'grid':
    inp.main_log(logfile, '{:-^100} \n\n'.format(' SPOTFINDING GRID SEARCH '))

    # run grid search on multiple processes
    cmd.Command.start("Spotfinding Grid Search")
    parallel_map(iterable=mp_input_list,
                 func=index_mproc_wrapper,
                 processes=gs_params.n_processors,
                 preserve_order = True)
    cmd.Command.end("Spotfinding Grid Search -- DONE")

  elif int_type == 'final':
    inp.main_log(logfile, "\n\n{:-^80}\n".format(' FINAL INTEGRATION '))

    cmd.Command.start("Final integration")
    result_objects = parallel_map(iterable=mp_input_list,
                                  func=final_mproc_wrapper,
                                  processes=gs_params.n_processors,
                                  preserve_exception_message=False)
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


def output_cleanup(gs_params):

  int_list_file = os.path.abspath('{}/integrated.lst'.format(gs_params.output))
  dest_dir = os.path.abspath("{}/integrated".format(gs_params.output))
  os.makedirs(dest_dir)

  with open(int_list_file, 'r') as int_file:
    int_file_contents = int_file.read()
    int_list = int_file_contents.splitlines()

  for int_file in int_list:
    filename = os.path.basename(int_file)
    dest_file = "{}/{}".format(dest_dir, filename)
    shutil.copyfile(int_file, dest_file)

  for tmp_dir in os.listdir(os.path.abspath(gs_params.output)):
    if "tmp" in tmp_dir:
      rmd_dir = os.path.join(os.path.abspath(gs_params.output), tmp_dir)
      shutil.rmtree(rmd_dir)

  os.remove(int_list_file)

  for int_file in os.listdir(dest_dir):
    with open(int_list_file, 'a') as f_int:
        f_int.write('{}\n'.format(int_file))


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

def check_options(args):
  """ Checks command-line options and runs anything that is pertinent
  """

  if args.l:
    list_file = os.path.abspath("{}/input.lst".format(gs_params.output))
    print '\nIOTA will run in LIST INPUT ONLY mode'
    print 'Input list in {} \n\n'.format(list_file)
    with open(list_file, "a") as lf:
      for input_file in input_list:
        print input_file
        lf.write('{}\n'.format(input_file))
    print '\nExiting...\n\n'
    sys.exit()

  if args.c:
    sys.exit()

def experimental(mp_input_list, gs_params, log_dir):
  """EXPERIMENTAL SECTION: Contains stuff I just want to try out without
     running the whole darn thing
  """
  print "IT WORKS!"

# ============================================================================ #

if __name__ == "__main__":

  iota_version = '1.50'
  now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
  logo = "\n\n"\
   "     IIIIII          OOOO         TTTTTTTTTT           A              \n"\
   "       II           O    O            TT              A A             \n"\
   "       II           O    O            TT             A   A            \n"\
   ">------INTEGRATION--OPTIMIZATION------TRIAGE--------ANALYSIS--------->\n"\
   "       II           O    O            TT           A       A          \n"\
   "       II           O    O            TT          A         A         \n"\
   "     IIIIII          OOOO             TT         A           A   v{}"\
   "".format(iota_version)



  # Read arguments
  parser = argparse.ArgumentParser(prog = 'prime.iota',
            description=('Integration Optimization, Triage and Analysis'))
  parser.add_argument('path', type=str, nargs = '?', default = None,
            help = 'Path to data or file with IOTA parameters (blank = dry run)')
  parser.add_argument('--version', action = 'version',
            version = 'IOTA {}'.format(iota_version),
            help = 'Prints version info of IOTA')
  parser.add_argument('-l', action = 'store_true',
            help = 'Use to ONLY make and print to file a list of input images')
  parser.add_argument('-o', action = 'store_true',
            help = 'Split output into separate files based on space group')
  parser.add_argument('-c', action = 'store_true',
            help = 'Convert files only & output to temp input folder')

  args = parser.parse_args()
  print logo

  if args.path == None:
    dry_run()
    sys.exit()
  else:
    print '\n{}\n'.format(now)
    carg = os.path.abspath(args.path)
    if os.path.exists(carg):

      # If user provided a parameter file
      if os.path.isfile(carg) and os.path.basename(carg).endswith('.param'):
        gs_params, txt_out = inp.process_input([carg])

      # If user provided a data folder
      elif os.path.isdir(carg):
        print "\nIOTA will run in AUTO mode:\n"
        cmd.Command.start("Generating default parameters")
        gs_params, txt_out = inp.auto_mode(os.path.abspath(os.path.curdir),
                                           os.path.abspath(carg), now)
        cmd.Command.end("Generating default parameters -- DONE")

    # If user provided gibberish
    else:
      print "ERROR: Invalid input! Need parameter filename or data folder."
      sys.exit()

  initial_input_list = inp.make_input_list(gs_params)
  input_list = [i for i in initial_input_list if \
                i.endswith(('pickle', 'mccd', 'cbf', 'img'))]

  # Check if input needs to be converted to pickle format; also check if input
  # images need to be cropped / padded to be square, w/ beam center in the
  # center of image. If these steps are needed, carry them out
  img_check = i2p.check_image(input_list[0])
  if img_check == 'image' or img_check == 'raw pickle':
    square = gs_params.advanced.square_mode
    converted_img_list, input_folder = inp.make_raw_input(input_list, gs_params)
    raw_input_list = zip(input_list, converted_img_list)

    cmd.Command.start("Converting {} images".format(len(raw_input_list)))
    parallel_map(iterable=raw_input_list,
                 func=conversion_wrapper,
                 processes=gs_params.n_processors)
    cmd.Command.end("Converting {} images -- DONE ".format(len(raw_input_list)))

    input_list = converted_img_list
  elif img_check == 'converted pickle':
    input_folder = gs_params.input
  else:
    print "ERROR: Unknown image format. Please check your input"
    sys.exit()

  check_options(args)

  # generate general input
  gs_range, input_dir_list, output_dir_list, log_dir, logfile, mp_input_list,\
  mp_output_list = inp.generate_input(gs_params, input_list, input_folder)

  iota_start(txt_out, gs_params, gs_range, input_dir_list,
                    input_list, mp_input_list)

  # debugging/experimental section - anything goes here
  if gs_params.advanced.experimental:
    print "\nIOTA will run in EXPERIMENTAL mode:\n"
    experimental(mp_input_list, gs_params, log_dir)
    sys.exit()

  # run grid search
  if gs_params.grid_search.flag_on:
    n_int = len(mp_input_list)
    prog_count = 0
    run_integration("grid", gs_params, mp_input_list)

  # run pickle selection
  selection_results = run_selection('grid', gs_params, mp_output_list)
  sel_clean = [entry for entry in selection_results \
                     if entry != None and entry != []]
  if sel_clean == []:
    print "\nNO IMAGES INTEGRATED! Check input and try again.\n"
    inp.main_log(logfile,"\nNO IMAGES INTEGRATED!\n")
    sys.exit()

  # run final integration
  final_int = run_integration("final", gs_params, sel_clean)

  # print final integration results and summary
  ia.print_results(final_int, gs_range, logfile)
  ia.print_summary(gs_params, len(input_list), logfile, iota_version, now)

  if gs_params.advanced.clean_up_output:
    if os.path.isfile(os.path.abspath("{}/integrated.lst".format(gs_params.output))):
      output_cleanup(gs_params)

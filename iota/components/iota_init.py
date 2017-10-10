from __future__ import division

from builtins import range
from builtins import object
'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 04/13/2017
Description : Reads command line arguments. Initializes all IOTA starting
              parameters. Starts main log.
'''
from __future__ import print_function

import os
import sys
import argparse
import time

import iota.components.iota_input as inp
import iota.components.iota_cmd as cmd
import iota.components.iota_misc as misc
from iota.components.iota_utils import InputFinder

ginp = InputFinder()
pid = os.getpid()

try:
  user = os.getlogin()
except OSError:
  user = 'iota'

# --------------------------- Initialize IOTA -------------------------------- #


def parse_command_args(iver, help_message):
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog = 'iota.run',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(help_message),
            epilog=('\n{:-^70}\n'.format('')))
  parser.add_argument('path', type=str, nargs = '*', default = None,
            help = 'Path to data or file with IOTA parameters')
  parser.add_argument('--version', action = 'version',
            version = 'IOTA {}'.format(iver),
            help = 'Prints version info of IOTA')
  parser.add_argument('-l', '--list', action = 'store_true',
            help = 'Output a file (input.lst) with input image paths and stop')
  parser.add_argument('-w', '--watch', action = 'store_true',
            help = 'Run IOTA in watch mode - check for new images')
  parser.add_argument('-c', '--convert', action = 'store_true',
            help = 'Convert raw images to pickles and stop')
  parser.add_argument('-f', '--full', action = 'store_true',
            help = 'Run IOTA in "full-processing" mode (advanced)')
  parser.add_argument('-d', '--default', action = 'store_true',
            help = 'Generate default settings files and stop')
  parser.add_argument('-p', '--prefix', type=str, default="Auto",
            help = 'Specify custom prefix for converted pickles')
  parser.add_argument('-s', '--select', action = 'store_true',
            help = 'Selection only, no grid search')
  parser.add_argument('-r', type=int, nargs=1, default=0, dest='random',
            help = 'Run IOTA with a random subset of images, e.g. "-r 5"')
  parser.add_argument('-n', type=int, nargs=1, default=0, dest='nproc',
            help = 'Specify a number of cores for a multiprocessor run"')
  parser.add_argument('--mpi', type=str, nargs='?', const=None, default=None,
            help = 'Specify stage of process - for MPI only')
  parser.add_argument('--analyze', type=str, nargs='?', const=None, default=None,
            help = 'Use for analysis only; specify run number or folder')
  return parser

class InitAll(object):
  """ Class to initialize current IOTA run

      help_message = description (hard-coded)

  """

  def __init__(self, help_message):
    self.iver = misc.iota_version
    self.user_id = user
    self.now = misc.now
    self.logo = "\n\n"\
   "     IIIIII            OOOOOOO        TTTTTTTTTT          A                 \n"\
   "       II             O       O           TT             A A                \n"\
   "       II             O       O           TT            A   A               \n"\
   ">------INTEGRATION----OPTIMIZATION--------TRIAGE-------ANALYSIS------------>\n"\
   "       II             O       O           TT          A       A             \n"\
   "       II             O       O           TT         A         A            \n"\
   "     IIIIII            OOOOOOO            TT        A           A   v{}     \n"\
   "".format(self.iver)
    self.help_message = self.logo + help_message
    self.input_base = None
    self.conv_base = None
    self.obj_base = None
    self.int_base = None

  def make_input_list(self):
    """ Reads input directory or directory tree and makes lists of input images.
        Optional selection of a random subset
    """

    # Read input from provided folder(s) or file(s)
    cmd.Command.start("Reading input files")
    input_entries = [i for i in self.params.input if i != None]
    input_list = ginp.make_input_list(input_entries)
    cmd.Command.end("Reading input files -- DONE")

    if len(input_list) == 0:
      print("\nERROR: No data found!")
      sys.exit()

    # Pick a randomized subset of images
    if self.params.advanced.random_sample.flag_on and \
       self.params.advanced.random_sample.number < len(input_list):
      inp_list = self.select_random_subset(input_list)
    else:
      inp_list = input_list

    return inp_list


  def select_random_subset(self, input_list):
    """ Selects random subset of input entries """
    import random

    random_inp_list = []
    if self.params.advanced.random_sample.number == 0:
      if len(input_list) <= 5:
        random_sample_number = len(input_list)
      elif len(input_list) <= 50:
        random_sample_number = 5
      else:
        random_sample_number = int(len(input_list) * 0.1)
    else:
      random_sample_number = self.params.advanced.random_sample.number

    cmd.Command.start("Selecting {} random images out of {} found".format(random_sample_number, len(input_list)))
    for i in range(random_sample_number):
      random_number = random.randrange(0, len(input_list))
      if input_list[random_number] in random_inp_list:
        while input_list[random_number] in random_inp_list:
          random_number = random.randrange(0, len(input_list))
        random_inp_list.append(input_list[random_number])
      else:
        random_inp_list.append(input_list[random_number])
    cmd.Command.end("Selecting {} random images out of {} found -- DONE ".format(random_sample_number, len(input_list)))

    return random_inp_list


  def make_int_object_list(self):
    """ Generates list of image objects from previous grid search """
    from libtbx import easy_pickle as ep

    if self.params.cctbx.selection.select_only.grid_search_path == None:
      int_dir = misc.set_base_dir('integration', True)
    else:
      int_dir = self.params.cctbx.selection.select_only.grid_search_path

    img_objects = []

    cmd.Command.start("Importing saved grid search results")
    for root, dirs, files in os.walk(int_dir):
      for filename in files:
        found_file = os.path.join(root, filename)
        if found_file.endswith(('int')):
          obj = ep.load(found_file)
          img_objects.append(obj)
    cmd.Command.end("Importing saved grid search results -- DONE")

    # Pick a randomized subset of images
    if self.params.advanced.random_sample.flag_on and \
       self.params.advanced.random_sample.number < len(img_objects):
      gs_img_objects = self.select_random_subset(img_objects)
    else:
      gs_img_objects = img_objects

    return gs_img_objects


  def analyze_prior_results(self, analysis_source):
    """ Runs analysis of previous grid search / integration results, used in an
        analyze-only mode """

    from iota.components.iota_analysis import Analyzer

    if os.path.isdir(analysis_source):
      int_folder = os.path.abspath(analysis_source)
    else:
      try:
        int_folder = os.path.abspath(os.path.join(os.curdir,
                     'integration/{}/image_objects'.format(analysis_source)))
      except ValueError:
        print('Run #{} not found'.format(analysis_source))

    if os.path.isdir(int_folder):

      cmd.Command.start("Analyzing results in {}".format(int_folder))
      int_list = [os.path.join(int_folder, i) for i in os.listdir(int_folder)]
      img_objects = [ep.load(i) for i in int_list if i.endswith('.int')]
      cmd.Command.end("Analyzing results -- DONE")

      self.logfile = os.path.abspath(os.path.join(int_folder, 'iota.log'))
      self.viz_base = os.path.join('/'.join(int_folder.split('/')),
                                   'vizualization')

      analysis = Analyzer(self, img_objects, self.iver)
      analysis.print_results()
      analysis.unit_cell_analysis(write_files=False)
      analysis.print_summary(write_files=False)
    else:
      print('No results found in {}'.format(int_folder))


  def run(self):

    self.args, self.phil_args = parse_command_args(self.iver,
                              self.help_message).parse_known_args()

    # Check for type of input
    if len(self.args.path) == 0 or self.args.path is None:  # No input
      parse_command_args(self.iver, self.help_message).print_help()
      if self.args.default:                      # Write out default params and exit
        help_out, txt_out = inp.print_params()
        print('\n{:-^70}\n'.format('IOTA Parameters'))
        print(help_out)
        inp.write_defaults(os.path.abspath(os.path.curdir), txt_out)
      misc.iota_exit()
    elif len(self.args.path) > 1:  # If multiple paths / wildcards
      file_list = ginp.make_input_list(self.args.path)
      list_file = os.path.join(os.path.abspath(os.path.curdir), 'input.lst')
      with open(list_file, 'w') as lf:
        lf.write('\n'.join(file_list))
      msg = "\nIOTA will run in AUTO mode using wildcard datapath:\n" \
            "{} files found, compiled in {}\n".format(len(file_list), list_file)
      self.params, self.txt_out = inp.process_input(self.args,
                                                    self.phil_args,
                                                    list_file,
                                                    'auto',
                                                    self.now)
    else:                                   # If single path, check type
      carg = os.path.abspath(self.args.path[0])
      if os.path.exists(carg):

        # If user provided a parameter file
        if os.path.isfile(carg) and os.path.basename(carg).endswith('.param'):
          msg = ''
          self.params, self.txt_out = inp.process_input(self.args,
                                                        self.phil_args,
                                                        carg, 'file')

        # If user provided a list of input files
        elif os.path.isfile(carg) and os.path.basename(carg).endswith('.lst'):
          msg = "\nIOTA will run in AUTO mode using {}:\n".format(carg)
          self.params, self.txt_out = inp.process_input(self.args,
                                                        self.phil_args,
                                                        carg, 'auto', self.now)

        # If user provided a single filepath
        elif os.path.isfile(carg) and not os.path.basename(carg).endswith('.lst'):
          msg = "\nIOTA will run in SINGLE-FILE mode using {}:\n".format(carg)
          self.params, self.txt_out = inp.process_input(self.args,
                                                        self.phil_args,
                                                        carg, 'auto', self.now)

        # If user provided a data folder
        elif os.path.isdir(carg):
          msg = "\nIOTA will run in AUTO mode using {}:\n".format(carg)
          self.params, self.txt_out = inp.process_input(self.args,
                                                        self.phil_args,
                                                        carg, 'auto', self.now)
      # If user provided gibberish
      else:
        print(self.logo)
        print("ERROR: Invalid input! Need parameter filename or data folder.")
        misc.iota_exit()

    # Identify indexing / integration program
    if self.params.advanced.integrate_with == 'cctbx':
      prg = "                                                             with CCTBX.XFEL\n"
    elif self.params.advanced.integrate_with == 'dials':
      prg = "                                                                  with DIALS\n"

    self.logo += prg
    print(self.logo)
    print('\n{}\n'.format(self.now))
    if msg != '':
      print(msg)

    if self.args.analyze != None:
      print('ANALYSIS ONLY will be performed (analyzing run #{})'.format(
        self.args.analyze))
      self.analyze_prior_results('{:003d}'.format(int(self.args.analyze)))
      misc.iota_exit()

    if self.params.mp_method == 'mpi':
      rank, size = misc.get_mpi_rank_and_size()
      self.master_process = rank == 0
    else:
      self.master_process = True

    # Call function to read input folder structure (or input file) and
    # generate list of image file paths
    if self.params.cctbx.selection.select_only.flag_on:
      self.gs_img_objects = self.make_int_object_list()
      self.input_list = [i.conv_img for i in self.gs_img_objects]
    else:
      self.input_list = self.make_input_list()

    # Check for -l option, output list of input files and exit
    if self.args.list:
      list_file = os.path.abspath("{}/input.lst".format(os.curdir))

      # Check if other files of this name exist under the current folder
      list_folder = os.path.dirname(list_file)
      list_files = [i for i in os.listdir(list_folder) if i.endswith(".lst")]
      if len(list_files) > 0:
        list_file = os.path.join(list_folder,
                                 "input_{}.lst".format(len(list_files)))

      print('\nINPUT LIST ONLY option selected')
      print('Input list in {} \n\n'.format(list_file))
      with open(list_file, "w") as lf:
        for i, input_file in enumerate(self.input_list, 1):
          lf.write('{}\n'.format(input_file))
          print("{}: {}".format(i, input_file))
          lf.write('{}\n'.format(input_file))
      print('\nExiting...\n\n')
      misc.iota_exit()

    # If fewer images than requested processors are supplied, set the number of
    # processors to the number of images
    if self.params.n_processors > len(self.input_list):
      self.params.n_processors = len(self.input_list)

    # Generate base folder paths
    self.conv_base = misc.set_base_dir('converted_pickles', out_dir = self.params.output)
    self.int_base = misc.set_base_dir('integration', out_dir = self.params.output)
    self.obj_base = os.path.join(self.int_base, 'image_objects')
    self.fin_base = os.path.join(self.int_base, 'final')
    self.log_base = os.path.join(self.int_base, 'logs')
    self.viz_base = os.path.join(self.int_base, 'visualization')
    self.tmp_base = os.path.join('/tmp', '{}_{}'.format(os.getlogin(), time.time()))

    # Generate base folders
    os.makedirs(self.int_base)
    os.makedirs(self.obj_base)
    os.makedirs(self.fin_base)
    os.makedirs(self.log_base)
    os.makedirs(self.tmp_base)

    # Determine input base
    self.input_base = os.path.abspath(os.path.dirname(os.path.commonprefix(self.input_list)))

    # Initialize main log
    self.logfile = os.path.abspath(os.path.join(self.int_base, 'iota.log'))

    # Log starting info
    misc.main_log(self.logfile, '{:=^80} \n'.format(' IOTA MAIN LOG '))
    misc.main_log(self.logfile, '{:-^80} \n'.format(' SETTINGS FOR THIS RUN '))
    misc.main_log(self.logfile, self.txt_out)

    if self.params.advanced.integrate_with == 'cctbx':
      target_file = self.params.cctbx.target
    elif self.params.advanced.integrate_with == 'dials':
      target_file = self.params.dials.target
    misc.main_log(self.logfile, '{:-^80} \n\n'
                                ''.format(' TARGET FILE ({}) CONTENTS '
                                ''.format(target_file)))
    with open(target_file, 'r') as phil_file:
      phil_file_contents = phil_file.read()
    misc.main_log(self.logfile, phil_file_contents)


# ============================================================================ #
if __name__ == "__main__":

  iota_version = '1.0.001G'
  help_message = ""

  initialize = InitAll(iota_version, help_message)
  initialize.run()

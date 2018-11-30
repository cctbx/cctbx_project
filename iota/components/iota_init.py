from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 11/29/2018
Description : Reads command line arguments. Initializes all IOTA starting
              parameters. Starts main log.
'''

import os
import argparse
import time
assert time
from contextlib import contextmanager

import dials.util.command_line as cmd

import iota.components.iota_input as inp
from iota.components.iota_utils import InputFinder, get_mpi_rank_and_size
from iota.components.iota_base import InitBase


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
  parser.add_argument('-f', '--full', action = 'store_true',
            help = 'Run IOTA in "full-processing" mode (advanced)')
  parser.add_argument('-d', '--default', action = 'store_true',
            help = 'Generate default settings files and stop')
  parser.add_argument('--ha14', action = 'store_true',
            help = 'Run IOTA with old HA14 backend')
  parser.add_argument('-random', type=int, nargs=1, default=0,
            help = 'Size of randomized subset, e.g. "--random 10"')
  parser.add_argument('--range', type=str, nargs='?', default=None,
            help = 'Range of images, e.g."--range 1-5,25,200-250"')
  parser.add_argument('-n', '--nproc', type=int, nargs=1, default=0,
            help = 'Specify a number of cores for a multiprocessor run"')
  parser.add_argument('--mpi', type=str, nargs='?', const=None, default=None,
            help = 'Specify stage of process - for MPI only')
  parser.add_argument('--analyze', type=str, nargs='?', const=None, default=None,
            help = 'Use for analysis only; specify run number or folder')
  return parser

@contextmanager    # Will print start / stop messages around some processes
def prog_message(msg):
  cmd.Command.start(msg)
  yield
  cmd.Command.end('{} -- DONE'.format(msg))


class XInitAll(InitBase):
  """ Class to initialize IOTA parameters for command-line jobs """
  def __init__(self, help_message):
    InitBase.__init__(self)

    self.logo = "\n\n" \
                "     IIIIII            OOOOOOO        TTTTTTTTTT          A                 \n" \
                "       II             O       O           TT             A A                \n" \
                "       II             O       O           TT            A   A               \n" \
                ">------INTEGRATION----OPTIMIZATION--------TRIAGE-------ANALYSIS------------>\n" \
                "       II             O       O           TT          A       A             \n" \
                "       II             O       O           TT         A         A            \n" \
                "     IIIIII            OOOOOOO            TT        A           A   v{}     \n" \
                "".format(self.iver)
    self.help_message = self.logo + help_message

  def analyze_prior_results(self, analysis_source):
    """ Runs analysis of previous grid search / integration results, used in an
        analyze-only mode """
    #TODO: move to iota_run.py

    from iota.components.iota_analysis import Analyzer
    from libtbx import easy_pickle as ep

    if os.path.isdir(analysis_source):
      int_folder = os.path.abspath(analysis_source)
    else:
      try:
        int_folder = os.path.abspath(os.path.join(os.curdir,
                     'integration/{}/image_objects'.format(analysis_source)))
      except ValueError:
        int_folder = None
        print ('Run #{} not found'.format(analysis_source))

    if os.path.isdir(int_folder):
      with prog_message('Analyzing Results'):
        int_list = [os.path.join(int_folder, i) for i in os.listdir(int_folder)]
        img_objects = [ep.load(i) for i in int_list if i.endswith('.int')]

      self.logfile = os.path.abspath(os.path.join(int_folder, 'iota.log'))
      self.viz_base = os.path.join('/'.join(int_folder.split('/')),
                                   'vizualization')

      self.params.analysis.cluster_write_files=False

      analysis = Analyzer(self, img_objects, self.iver)
      analysis.print_results()
      analysis.unit_cell_analysis()
      analysis.print_summary(write_files=False)
    else:
      print ('No results found in {}'.format(int_folder))

  def initialize_interface(self):
    """ Initialize properties specific for command line

    :return: True if successful, False if fails
    """
    self.args, self.phil_args = parse_command_args(self.iver,
                                                   self.help_message).parse_known_args()
    ginp = InputFinder()

    # Check for type of input
    if not self.args.path:  # No input
      parse_command_args(self.iver, self.help_message).print_help()
      if self.args.default:  # Write out default params and exit
        help_out, txt_out = inp.print_params()
        print('\n{:-^70}\n'.format('IOTA Parameters'))
        print(help_out)
      return False, 'IOTA_XTERM_INIT: OUTPUT PARAMETERS ONLY'
    elif len(self.args.path) > 1:  # If multiple paths / wildcards
      file_list = ginp.make_input_list(self.args.path)
      list_file = os.path.join(os.path.abspath(os.path.curdir), 'input.lst')
      with open(list_file, 'w') as lf:
        lf.write('\n'.join(file_list))
      msg = "\nIOTA will run in AUTO mode using wildcard datapath:\n" \
            "{} files found, compiled in {}\n".format(len(file_list), list_file)
      self.iota_phil = inp.process_input(self.args, self.phil_args, list_file,
                                         'auto', self.now)
      self.params = self.iota_phil.extract()

    else:  # If single path, check type
      carg = os.path.abspath(self.args.path[0])
      if os.path.isfile(carg):
        ptype = ginp.get_file_type(carg)
        if ptype.lower() in ('raw image', 'image pickle'):
          msg = "\nIOTA will run in SINGLE-FILE mode using {}:\n".format(carg)
          mode = 'auto'
        elif ('iota' and 'settings' in ptype.lower()):
          msg = '\nIOTA will run in SCRIPT mode using {}:\n'.format(carg)
          mode = 'file'
        elif 'list' in ptype.lower():
          msg = "\nIOTA will run in AUTO mode using {}:\n".format(carg)
          mode = 'auto'
        else:
          pr = 'WARNING! File format not recognized. Proceed anyway? [Y/N] '
          unknown_file = raw_input(pr)
          if 'y' in unknown_file.lower():
            ftype = raw_input("File type? [image, list, or parameters] ")
            msg = "\nIOTA will run WITH DODGY input using {}:\n".format(carg)
            if 'par' in ftype:
              mode = 'file'
            else:
              mode = 'auto'
          else:
            print('Exiting...')
            return False, 'IOTA_XTERM_INIT_ERROR: Unrecognizable input!'
      elif os.path.isdir(carg):
        ptype = ginp.get_folder_type(carg)
        if ('image' and 'folder' in ptype.lower()):
          msg = "\nIOTA will run in AUTO mode using {}:\n".format(carg)
          mode = 'auto'
        else:
          msg = "IOTA_XTERM_INIT_ERROR: No images in {}!".format(carg)
          print(self.logo)
          print(msg)
          return False, msg

      # If user provided gibberish
      else:
        msg = "IOTA_XTERM_INIT_ERROR: Invalid input! Need parameter filename " \
              "or data folder."
        print(self.logo)
        print(msg)
        return False, msg

      # Initialize parameters for this command-line run
      self.iota_phil = inp.process_input(self.args, self.phil_args,
                                         carg, mode, self.now)
      self.params = self.iota_phil.extract()

    # Identify indexing / integration program and add to logo
    b_end = " with {}".format(str(self.params.advanced.processing_backend).upper())
    prg = "{:>{w}}".format(b_end, w=76)
    self.logo += prg
    print(self.logo)
    print('\n{}\n'.format(self.now))
    if msg != '':
      print(msg)

    if self.args.analyze is not None:
      print('ANALYSIS ONLY will be performed (analyzing run #{})'.format(
        self.args.analyze))
      self.analyze_prior_results('{:003d}'.format(int(self.args.analyze)))
      return False

    if self.params.mp.method == 'mpi':
      rank, size = get_mpi_rank_and_size()
      self.master_process = rank == 0
    else:
      self.master_process = True

    # Call function to read input folder structure (or input file) and
    # generate list of image file paths

    with prog_message("Reading input files"):
      self.input_list = self.make_input_list()

    # Select range of images/objects if turned on
    if self.params.advanced.image_range.flag_on:
      self.input_list = self.select_image_range(self.input_list)

    # Pick a randomized subset of images/objects if turned on
    if self.params.advanced.random_sample.flag_on and \
            self.params.advanced.random_sample.number < len(self.input_list):
      with prog_message("Selecting {} random images out of {} found"
                        "".format(self.params.advanced.random_sample.number,
                                  len(self.input_list))):
        self.input_list = self.select_random_subset(self.input_list)

      # Check for -l option, output list of input files and exit
    if self.args.list:
      list_file = os.path.abspath("{}/input.lst".format(os.curdir))

      # Check if other files of this name exist under the current folder
      list_folder = os.path.dirname(list_file)
      list_files = [i for i in os.listdir(list_folder) if i.endswith(".lst")]
      if len(list_files) > 0:
        list_file = os.path.join(list_folder,
                                 "input_{}.lst".format(len(list_files)))

      msg = 'IOTA_XTERM_INIT: INPUT LIST ONLY option selected'
      print ('\n{}'.format(msg))
      print ('Input list in {} \n\n'.format(list_file))
      with open(list_file, "w") as lf:
        for i, input_file in enumerate(self.input_list, 1):
          lf.write('{}\n'.format(input_file))
          print ("{}: {}".format(i, input_file))
          lf.write('{}\n'.format(input_file))
      print ('\nExiting...\n\n')
      return False, msg

    return True, 'IOTA_XTERM_INIT: Initialization complete!'

  def initialize_info_object(self):
    # Add input to info object
    self.info.img_list = [[i, len(self.input_list) + 1, j] for
                          i, j in enumerate(self.input_list, 1)]
    self.info.nref_list = [0] * len(self.info.img_list)
    self.info.nref_xaxis = [i[0] for i in self.info.img_list]
    self.info.res_list = [0] * len(self.info.img_list)

  def run(self):
    ''' Override of base INIT class to account for XTerm specifics

    :return: True if okay, False if not, and a message (ok or error)
    '''

    # Interface-specific options
    ui_init_good, msg = self.initialize_interface()
    if not ui_init_good:
      return False, msg

    # Create output file structure
    init_out, msg = self.initialize_output()
    if not init_out:
      return False, msg

    # Initialize IOTA and backend parameters
    init_param, msg = self.initialize_parameters()
    if not init_param:
      return False, msg

    # Initialize info object
    self.initialize_info_object()

    # Initalize main log (iota.log)
    init_log, msg = self.initialize_main_log()
    if not init_log:
      return False, msg

    # Return True for successful initialization
    return True, msg



# ============================================================================ #
if __name__ == "__main__":

  iota_version = '1.0.001G'
  help_message = ""

  initialize = InitBase()
  initialize.run()

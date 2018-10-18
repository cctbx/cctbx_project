from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 10/18/2018
Description : Reads command line arguments. Initializes all IOTA starting
              parameters. Starts main log.
'''

import os
import argparse
import time
assert time

import dials.util.command_line as cmd

import iota.components.iota_input as inp
from iota.components.iota_utils import InitAll, InputFinder, get_mpi_rank_and_size

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


class XInitAll(InitAll):
  """ Class to initialize IOTA parameters for command-line jobs """
  def __init__(self, help_message):
    InitAll.__init__(self)

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

      cmd.Command.start("Analyzing results in {}".format(int_folder))
      int_list = [os.path.join(int_folder, i) for i in os.listdir(int_folder)]
      img_objects = [ep.load(i) for i in int_list if i.endswith('.int')]
      cmd.Command.end("Analyzing results -- DONE")

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
    if len(self.args.path) == 0 or self.args.path is None:  # No input
      parse_command_args(self.iver, self.help_message).print_help()
      if self.args.default:  # Write out default params and exit
        help_out, txt_out = inp.print_params()
        print('\n{:-^70}\n'.format('IOTA Parameters'))
        print(help_out)
        inp.write_defaults(os.path.abspath(os.path.curdir), txt_out)
      return False
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
        if ('image' and 'file' in ptype.lower()):
          msg = "\nIOTA will run in SINGLE-FILE mode using {}:\n".format(carg)
          mode = 'auto'
        elif ('iota' and 'settings' in ptype.lower()):
          msg = ''
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
            return False
      elif os.path.isdir(carg):
        ptype = ginp.get_folder_type(carg)
        if ('image' and 'folder' in ptype.lower()):
          msg = "\nIOTA will run in AUTO mode using {}:\n".format(carg)
          mode = 'auto'
        else:
          print(self.logo)
          print("ERROR: No images in {}!".format(carg))
          return False

      # If user provided gibberish
      else:
        msg = "ERROR: Invalid input! Need parameter filename or data folder."
        print(self.logo)
        print(msg)
        return False

      # Initialize parameters for this command-line run
      self.iota_phil = inp.process_input(self.args, self.phil_args,
                                         carg, mode, self.now)
      self.params = self.iota_phil.extract()

    # Identify indexing / integration program and add to logo
    b_end = " with {}".format(str(self.params.advanced.integrate_with).upper())
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

    if self.params.mp_method == 'mpi':
      rank, size = get_mpi_rank_and_size()
      self.master_process = rank == 0
    else:
      self.master_process = True

    # Call function to read input folder structure (or input file) and
    # generate list of image file paths
    if self.params.cctbx.selection.select_only.flag_on:
      cmd.Command.start("Importing saved grid search results")
      self.gs_img_objects = self.make_int_object_list()
      self.input_list = [i.conv_img for i in self.gs_img_objects]
      cmd.Command.end("Importing saved grid search results -- DONE")
    else:
      cmd.Command.start("Reading input files")
      self.input_list = self.make_input_list()
      cmd.Command.end("Reading input files -- DONE")

    # Select range of images/objects if turned on
    if self.params.advanced.image_range.flag_on:
      self.input_list = self.select_image_range(self.input_list)

    # Pick a randomized subset of images/objects if turned on
    if self.params.advanced.random_sample.flag_on and \
            self.params.advanced.random_sample.number < len(self.input_list):
      cmd.Command.start("Selecting {} random images out of {} found"
                        "".format(self.params.advanced.random_sample.number,
                                  len(self.input_list)))
      self.input_list = self.select_random_subset(self.input_list)
      cmd.Command.end("Selecting {} random images out of {} found -- DONE"
                        "".format(self.params.advanced.random_sample.number,
                                  len(self.input_list)))

      # Check for -l option, output list of input files and exit
    if self.args.list:
      list_file = os.path.abspath("{}/input.lst".format(os.curdir))

      # Check if other files of this name exist under the current folder
      list_folder = os.path.dirname(list_file)
      list_files = [i for i in os.listdir(list_folder) if i.endswith(".lst")]
      if len(list_files) > 0:
        list_file = os.path.join(list_folder,
                                 "input_{}.lst".format(len(list_files)))

      print ('\nINPUT LIST ONLY option selected')
      print ('Input list in {} \n\n'.format(list_file))
      with open(list_file, "w") as lf:
        for i, input_file in enumerate(self.input_list, 1):
          lf.write('{}\n'.format(input_file))
          print ("{}: {}".format(i, input_file))
          lf.write('{}\n'.format(input_file))
      print ('\nExiting...\n\n')
      return False

    return True



# ============================================================================ #
if __name__ == "__main__":

  iota_version = '1.0.001G'
  help_message = ""

  initialize = InitAll()
  initialize.run()

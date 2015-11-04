from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 10/28/2015
Description : Reads command line arguments. Initializes all IOTA starting
              parameters. Starts main log.
'''

import os
import sys
import argparse
from datetime import datetime

import prime.iota.iota_input as inp
import prime.iota.iota_cmd as cmd
import prime.iota.iota_misc as misc

# --------------------------- Initialize IOTA -------------------------------- #


def parse_command_args(iver, help_message):
  """ Parses command line arguments (only options for now) """
  parser = argparse.ArgumentParser(prog = 'prime.iota',
            formatter_class=argparse.RawDescriptionHelpFormatter,
            description=(help_message),
            epilog=('\n{:-^70}\n'.format('')))
  parser.add_argument('path', type=str, nargs = '?', default = None,
            help = 'Path to data or file with IOTA parameters')
  parser.add_argument('--version', action = 'version',
            version = 'IOTA {}'.format(iver),
            help = 'Prints version info of IOTA')
  parser.add_argument('-l', '--list', action = 'store_true',
            help = 'Output a file (input.lst) with input image paths and stop')
  parser.add_argument('-c', '--convert', action = 'store_true',
            help = 'Convert raw images to pickles and stop')
  parser.add_argument('-d', '--default', action = 'store_true',
            help = 'Generate default iota.param and target.phil files and stop')
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

      iver = IOTA version (hard-coded)
      help_message = description (hard-coded)

  """

  def __init__(self, iver, help_message):
    from datetime import datetime
    self.iver = iver
    self.now = "{:%A, %b %d, %Y. %I:%M %p}".format(datetime.now())
    self.logo = "\n\n"\
   "     IIIIII            OOOOOOO        TTTTTTTTTT          A              \n"\
   "       II             O       O           TT             A A             \n"\
   "       II             O       O           TT            A   A            \n"\
   ">------INTEGRATION----OPTIMIZATION--------TRIAGE-------ANALYSIS--------->\n"\
   "       II             O       O           TT          A       A          \n"\
   "       II             O       O           TT         A         A         \n"\
   "     IIIIII            OOOOOOO            TT        A           A   v{}  \n"\
   "".format(iver)
    self.help_message = self.logo + help_message
    self.input_base = None
    self.conv_base = None
    self.obj_base = None
    self.int_base = None

  def make_input_list(self):
    """ Reads input directory or directory tree and makes lists of input images
        (in pickle format) using absolute path for each file. If a separate file
        with list of images is provided, parses that file and uses that as the
        input list. If random input option is selected, pulls a specified number
        of random images from the list and outputs that subset as the input list.
    """
    input_entries = [i for i in self.params.input if i != None]
    input_list = []

    # run through the list of multiple input entries (or just the one) and
    # concatenate the input list
    for input_entry in input_entries:
      if os.path.isfile(input_entry):
        if input_entry.endswith('.lst'):          # read from file list
          cmd.Command.start("Reading input list from file")
          with open(input_entry, 'r') as listfile:
            listfile_contents = listfile.read()
          input_list.extend(listfile_contents.splitlines())
          cmd.Command.end("Reading input list from file -- DONE")
        elif input_entry.endswith(('pickle', 'mccd', 'cbf', 'img')):
          input_list.append(input_entry)             # read in image directly

      elif os.path.isdir(input_entry):
        abs_inp_path = os.path.abspath(input_entry)

        cmd.Command.start("Reading files from data folder")
        for root, dirs, files in os.walk(abs_inp_path):
          for filename in files:
            found_file = os.path.join(root, filename)
            if found_file.endswith(('pickle', 'mccd', 'cbf', 'img')):
              input_list.append(found_file)
        cmd.Command.end("Reading files from data folder -- DONE")

    if len(input_list) == 0:
      print "\nERROR: No data found!"
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

    from prime.iota.iota_analysis import Analyzer
    from libtbx import easy_pickle as ep

    if os.path.isdir(analysis_source):
      int_folder = os.path.abspath(analysis_source)
    else:
      try:
        int_folder = os.path.abspath(os.path.join(os.curdir,
                     'integration/{}/image_objects'.format(analysis_source)))
      except ValueError:
        print 'Run #{} not found'.format(analysis_source)

    if os.path.isdir(int_folder):
      int_list = [os.path.join(int_folder, i) for i in os.listdir(int_folder)]
      img_objects = [ep.load(i) for i in int_list]

      analysis = Analyzer(img_objects, None, self.iver, self.now)
      analysis.print_results()
      analysis.unit_cell_analysis(self.params.analysis.cluster_threshold,
                                  int_folder, False)
      analysis.print_summary(None)
      analysis.show_heatmap()
    else:
      print 'No results found in {}'.format(int_folder)



  # Runs general initialization
  def run(self):

    self.args, self.phil_args = parse_command_args(self.iver,
                                self.help_message).parse_known_args()

    # Check for type of input
    if self.args.path == None:                   # No input
      parse_command_args(self.iver, self.help_message).print_help()
      if self.args.default:                      # Write out default params and exit
        help_out, txt_out = inp.print_params()
        print '\n{:-^70}\n'.format('IOTA Parameters')
        print help_out
        inp.write_defaults(os.path.abspath(os.path.curdir), txt_out)
      misc.iota_exit(self.iver)
    else:                                   # If input exists, check type
      carg = os.path.abspath(self.args.path)
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
        print self.logo
        print "ERROR: Invalid input! Need parameter filename or data folder."
        misc.iota_exit(self.iver)

    # Identify indexing / integration program
    if self.params.advanced.integrate_with == 'cctbx':
      prg = "                                                          with CCTBX.XFEL\n"
    elif self.params.advanced.integrate_with == 'dials':
      prg = "                                                               with DIALS\n"

    self.logo += prg
    print self.logo
    print '\n{}\n'.format(self.now)
    if msg != '':
      print msg

    # Check for -l option, output list of input files and exit
    if self.args.list:
      list_file = os.path.abspath("{}/input.lst".format(os.curdir))
      print '\nIOTA will run in LIST INPUT ONLY mode'
      print 'Input list in {} \n\n'.format(list_file)
      with open(list_file, "w") as lf:
        for i, input_file in enumerate(input_list, 1):
          print "{}: {}".format(i, input_file)
          lf.write('{}\n'.format(input_file))
      print '\nExiting...\n\n'
      misc.iota_exit(self.iver)

    if self.args.analyze != None:
      self.analyze_prior_results('{:003d}'.format(int(self.args.analyze)))
      misc.iota_exit(self.iver)

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

    # If fewer images than requested processors are supplied, set the number of
    # processors to the number of images
    if self.params.n_processors > len(self.input_list):
      self.params.n_processors = len(self.input_list)

    # Generate base folder paths
    self.conv_base = misc.set_base_dir('converted_pickles')
    self.int_base = misc.set_base_dir('integration')
    self.obj_base = os.path.join(self.int_base, 'image_objects')
    self.fin_base = os.path.join(self.int_base, 'final')
    self.tmp_base = os.path.join(self.int_base, 'tmp')
    if self.params.analysis.viz != 'None' or\
       self.params.analysis.heatmap != 'None' or\
       self.params.analysis.charts:
      self.viz_base = os.path.join(self.int_base, 'visualization')
    else:
      self.viz_base = None

    # Generate base folders
    os.makedirs(self.int_base)
    os.makedirs(self.obj_base)
    os.makedirs(self.fin_base)
    os.makedirs(self.tmp_base)

    # Determine input base
    self.input_base = os.path.abspath(os.path.dirname(os.path.commonprefix(self.input_list)))

    # Initialize main log
    self.logfile = os.path.abspath(os.path.join(self.int_base, 'iota.log'))

    # Log starting info
    misc.main_log(self.logfile, '{:=^100} \n'.format(' IOTA MAIN LOG '))
    misc.main_log(self.logfile, '{:-^100} \n'.format(' SETTINGS FOR THIS RUN '))
    misc.main_log(self.logfile, self.txt_out)

    misc.main_log(self.logfile, '{:-^100} \n\n'
                                ''.format(' TARGET FILE ({}) CONTENTS '
                                ''.format(self.params.target)))
    with open(self.params.target, 'r') as phil_file:
      phil_file_contents = phil_file.read()
    misc.main_log(self.logfile, phil_file_contents)


# ============================================================================ #
if __name__ == "__main__":

  iota_version = '2.00'
  help_message = ""

  initialize = InitAll(iota_version, help_message)
  initialize.run()

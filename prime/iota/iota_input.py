from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 03/24/2015
Description : IOTA I/O module. Reads PHIL input, creates output directories,
              creates input lists and organizes starting parameters
'''


import sys
import os
import shutil
import random
from cStringIO import StringIO


import iotbx.phil as ip
import dials.util.command_line as cmd

master_phil = ip.parse("""
description = Integration optimization and transfer app (IOTA) input file
  .type = str
  .help = Run description (optional).
  .multiple = False
  .optional = True
input = None
  .type = path
  .multiple = False
  .help = Path to folder with raw data in pickle format. Can be a tree w/ subfolders
  .optional = False
input_list = None
  .type = str
  .multiple = False
  .help = Filename with list of specific subset of files under input to process (overrides input_list and input_dir_list)
  .optional = False
output = iota_output
  .type = str
  .help = Base name for output folder
  .optional = False
target = target.phil
  .type = str
  .multiple = False
  .help = Target (.phil) file with integration parameters
  .optional = False
advanced
  .help = "Advanced options, mostly for debugging."
{
  single_img = None
    .type = str
    .help = Runs grid search on specified single image
  charts = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity, rmsd, etc.
  mosaicity_plot = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity, rmsd, etc.
  debug = False
    .type = bool
    .help = If True, activates the iPython command wherever it is.
  save_tmp_pickles = False
    .type = bool
    .help = If True, saves pickle for each integration attempt.
  pred_img
    .help = "Visualize spotfinding / integration results."
  {
    flag = False
      .type = bool
      .help = Choose to visualize the results of final integration
    cv_vectors = True
      .type = bool
      .help = Include x,y offset information
  }
  random_sample
    .help = "Random grid search."
  {
    flag_on = False
      .type = bool
      .help = Set to run grid search on a random set of images.
    number = 5
      .type = int
      .help = Number of random samples.
  }
}
grid_search
  .help = "Parameters for the grid search."
{
  flag_on = True
    .type = bool
    .help = Set to False to run selection of pickles only.
  a_avg = 5
    .type = int
    .help = Minimum spot area.
  a_std = 4
    .type = int
    .help = Maximum spot area.
  h_avg = 5
    .type = int
    .help = Minimum spot height.
  h_std = 4
    .type = int
    .help = Maximum spot height.
}
flag_prefilter = True
  .type = bool
  .help = Activate space group / unit cell pre-filter.
target_unit_cell = 79.4, 79.4, 38.1, 90.0, 90.0, 90.0
  .type = unit_cell
  .help = Target unit-cell parameters are used to discard outlier cells.
target_space_group = P 43 21 2
  .type = str
  .help = Target space group.
target_pointgroup = P 4
  .type = str
  .help = Target point group.
target_uc_tolerance = 0.05
  .type = float
  .help = Maximum allowed unit cell deviation from target
min_sigma = 5
  .type = int
  .help = minimum sigma cutoff for "strong spots"
n_processors = 32
  .type = int
  .help = No. of processing units
""")

class Capturing(list):
  def __enter__(self):
    self._stdout = sys.stdout
    sys.stdout = self._stringio = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio.getvalue().splitlines())
    sys.stdout = self._stdout

def process_input(input_file_list):
  """ Read and parse parameter file

      input: input_file_list - PHIL-format files w/ parameters

      output: params - PHIL-formatted parameters
              txt_output - plain text-formatted parameters
  """

  user_phil = [ip.parse(open(inp).read()) for inp in input_file_list]

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  #capture input read out by phil
  with Capturing() as output:
    working_phil.show()

  txt_out = ''
  for one_output in output:
    txt_out += one_output + '\n'

  return params, txt_out


def make_input_list (gs_params):
  """ Reads input directory or directory tree and makes lists of input images
      (in pickle format) using absolute path for each file. If a separate file
      with list of images is provided, parses that file and uses that as the
      input list. If random input option is selected, pulls a specified number
      of random images from the list and outputs that subset as the input list.

      input: gs_params - parameters in PHIL format
      output: inp_list - list of input files
  """

  input_list = []
  abs_inp_path = os.path.abspath(gs_params.input)

  if gs_params.input_list != None:
    with open(gs_params.input_list, 'r') as listfile:
      listfile_contents = listfile.read()
    input_list = listfile_contents.splitlines()
  else:
  # search for *.pickle files within the tree and record in a list w/
  # full absolute path and filanames
    for root, dirs, files in os.walk(abs_inp_path):
      for filename in files:
        if filename.endswith("pickle"):
          pickle_file = os.path.join(root, filename)
          input_list.append(pickle_file)

  if gs_params.advanced.random_sample.flag_on == True:
    random_inp_list = []
    for i in range(gs_params.advanced.random_sample.number):
      random_number = random.randrange(0, len(input_list))
      random_inp_list.append(input_list[random_number])

    inp_list = random_inp_list
  else:
    inp_list = input_list

  return inp_list


def make_dir_lists(input_list, gs_params):
  """ From the input list, makes a list of input and output folders, such that
      the output directory structure mirrors the input directory structure, in
      case of duplication of image filenames.

      input: input_list - list of input files (w/ absolute paths)
             gs_params - parameters in PHIL format

      output: input_dir_list - list of input folders
              output_dir_list - list of output folders
  """

  input_dir_list = []
  output_dir_list = []

  abs_inp_path = os.path.abspath(gs_params.input)
  abs_out_path = os.path.abspath(gs_params.output)

  # make lists of input and output directories and files
  for input_entry in input_list:
    path = os.path.dirname(input_entry)

    if os.path.relpath(path, abs_inp_path) == '.':  # in case of input in one dir
      input_dir = abs_inp_path
      output_dir = abs_out_path
    else:                                           # in case of input in tree
      input_dir = abs_inp_path + '/' + os.path.relpath(path, abs_inp_path)
      output_dir = abs_out_path + '/' + os.path.relpath(path, abs_inp_path)

    if input_dir not in input_dir_list:
      input_dir_list.append(os.path.normpath(input_dir))
    if output_dir not in output_dir_list:
      output_dir_list.append(os.path.normpath(output_dir))

  return input_dir_list, output_dir_list


def make_mp_input(input_list, gs_params, gs_range):
  """ Generates input for multiprocessor grid search and selection.

      input: input_list - list of input images (w/ absolute paths)
             gs_params - list of parameters in PHIL format
             gs_range - grid search limits

      output: mp_input - list of input entries for MP grid search:
                1. raw image file (absolute path)
                2-4. signal height, spot height & spot area parameters
                5. output folder for integration result (absolute path)
              mp_output - list of entries for MP selection
                1. output folder for integration result (absolute path)
                2. raw image file (absolute path)
                3. integration result file (filename only)

              (The reason for duplication of items in the two lists has to do
              with the user being able to run selection / re-integration witout
              repeating the time-consuming grid search.)
  """

  mp_item = []
  mp_input = []
  mp_output = []

  for current_img in input_list:
    # generate output folder tree
    path = os.path.dirname(current_img)
    img_filename = os.path.basename(current_img)

    if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
      output_dir = os.path.abspath(gs_params.output)
    else:
      output_dir = '{0}/{1}'.format(os.path.abspath(gs_params.output),
                                    os.path.relpath(path,
                                    os.path.abspath(gs_params.input)))

    current_output_dir = "{0}/tmp_{1}".format(output_dir,
                                              img_filename.split('.')[0])
    mp_output_entry = [current_output_dir, current_img,
                       "int_{}.lst".format(img_filename.split('.')[0])]
    mp_output.append(mp_output_entry)

    # Create input list w/ filename and spot-finding params
    h_min = gs_range[0]
    h_max = gs_range[1]
    a_min = gs_range[2]
    a_max = gs_range[3]

    for sig_height in range(h_min, h_max + 1):
      for spot_area in range (a_min, a_max + 1):
        mp_item = [current_img, sig_height, sig_height, spot_area,
                   current_output_dir]
        mp_input.append(mp_item)

  return mp_input, mp_output

def make_dirs (mp_output_list, gs_params):

  # If grid-search turned on, check for existing output directory and remove
  if os.path.exists(os.path.abspath(gs_params.output)):
    cmd.Command.start("Deleting old folder {}".format(gs_params.output))
    shutil.rmtree(os.path.abspath(gs_params.output))
    cmd.Command.end("Deleting old folder {} -- DONE".format(gs_params.output))

  # Make main output directory and log directory
  os.makedirs(os.path.abspath(gs_params.output))
  os.makedirs("{}/logs".format(os.path.abspath(gs_params.output)))

  # Make per-image output folders. ***May not be necessary!!
  cmd.Command.start("Generating output directory tree")

  output_folders = [op[0] for op in mp_output_list]

  for folder in output_folders:
    if not os.path.exists(folder):
      os.makedirs(folder)

  cmd.Command.end("Generating output directory tree -- DONE")


def main_log_init(logfile):
  """ Save log from previous run and initiate a new log (run once). This is only
      necessary if re-running selection after grid-search
  """

  log_count = 0
  log_dir = os.path.dirname(logfile)

  for item in os.listdir(log_dir):
    if item.find('iota') != -1:
      log_count += 1

  if log_count > 0:
    old_log_filename = logfile
    new_log_filename = '{0}/iota_{1}.log'.format(log_dir, log_count)
    os.rename(old_log_filename, new_log_filename)

  with open(logfile, 'w') as logfile:
    logfile.write("IOTA LOG\n\n")


def main_log(logfile, entry):
  """ Write main log (so that I don't have to repeat this every time). All this
      is necessary so that I don't have to use the Python logger module, which
      creates a lot of annoying crosstalk with other cctbx.xfel modules.
  """

  with open(logfile, 'a') as logfile:
    logfile.write('{}\n'.format(entry))

def generate_input(gs_params):
  """ This section generates input for grid search and/or pickle selection.

      parameters: gs_params - list of parameters from *.param file (in
      PHIL format)

      output: gs_range - grid search range from avg and std-dev
              input_list - list of absolute paths to input files
              input_dir_list - list of absolute paths to input folder(s)
              output_dir_list - same for output folder(s)
              log_dir - log directory
              logfile - the log filename
              mp_input_list - for multiprocessing: filename + spotfinding params
              mp_output_list - same for output
  """

  # Determine grid search range from average and std. deviation params
  gs_range = [gs_params.grid_search.h_avg - gs_params.grid_search.h_std,
              gs_params.grid_search.h_avg + gs_params.grid_search.h_std,
              gs_params.grid_search.a_avg - gs_params.grid_search.a_std,
              gs_params.grid_search.a_avg + gs_params.grid_search.a_std
             ]

  # Make input list
  input_list = make_input_list(gs_params)

  # Make log directory and input/output directory lists
  log_dir = "{}/logs".format(os.path.abspath(gs_params.output))
  cmd.Command.start("Reading data folder(s)")
  input_dir_list, output_dir_list = make_dir_lists(input_list, gs_params)
  cmd.Command.end("Reading data folder(s) -- DONE")

  # Make input/output lists for multiprocessing
  cmd.Command.start("Generating multiprocessing input")
  mp_input_list, mp_output_list = make_mp_input(input_list, gs_params, gs_range)
  cmd.Command.end("Generating multiprocessing input -- DONE")

  # If grid-search turned on, check for existing output directory and remove
  if gs_params.grid_search.flag_on == True:
    make_dirs(mp_output_list, gs_params)
  else:
    if not os.path.exists(os.path.abspath(gs_params.output)):
      print "ERROR: No grid search results detected in"\
          "{}".format(os.path.abspath(gs_params.output))
      sys.exit()

  # Initiate log file
  logfile = '{}/iota.log'.format(log_dir)
  main_log_init(logfile)

  return gs_range, input_list, input_dir_list, output_dir_list, log_dir,\
         logfile, mp_input_list, mp_output_list


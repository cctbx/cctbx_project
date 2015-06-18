from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 06/12/2015
Description : IOTA I/O module. Reads PHIL input, creates output directories,
              creates input lists and organizes starting parameters, also
              creates reasonable IOTA and PHIL defaults if selected
'''


import sys
import os
import shutil
import random
from cStringIO import StringIO

import iotbx.phil as ip
import iota_cmd as cmd

master_phil = ip.parse("""
description = Integration Optimization, Transfer and Analysis (IOTA)
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
grid_search
  .help = "Parameters for the grid search."
{
  flag_on = True
    .type = bool
    .help = Set to False to run selection and final integration only. (Requires grid search results to be present!)
  a_avg = 5
    .type = int
    .help = Minimum spot area.
  a_std = 2
    .type = int
    .help = Maximum spot area.
  h_avg = 4
    .type = int
    .help = Minimum spot height.
  h_std = 2
    .type = int
    .help = Maximum spot height.
}
flag_prefilter = False
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
advanced
  .help = "Advanced, debugging and experimental options."
{
  sig_height_search = False
    .type = bool
    .help = scan signal height in addition to spot height
  single_img = None
    .type = str
    .help = Runs grid search on specified single image
  charts = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity, rmsd, etc.
  mosaicity_plot = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity only.
  debug = False
    .type = bool
    .help = Used for various debugging purposes.
  experimental = False
    .type = bool
    .help = Set to true to run the experimental section of codes
  save_tmp_pickles = False
    .type = bool
    .help = If True, saves pickle for each integration attempt in grid search.
  clean_up_output = False
    .type = bool
    .help = Set to False to leave temporary folders and files in place'
  cluster_threshold = 5000
    .type = int
    .help = threshold value for unit cell clustering
  convert_images = True
    .type = bool
    .help = Set to False to force non-conversion of images
  square_mode = None *pad crop
    .type = choice
    .help = Method to generate square image
  beamstop = 0
    .type = float
    .help = Beamstop shadow threshold, zero to skip
  pred_img
    .help = "Visualize spotfinding and integration results."
  {
    flag = False
      .type = bool
      .help = Visualize the results of final integration
    cv_vectors = False
      .type = bool
      .help = Visualize x,y offset information
  }
  random_sample
    .help = "Random grid search."
  {
    flag_on = False
      .type = bool
      .help = Set to run grid search on a random set of images.
    number = 5
      .type = int
      .help = Number of random samples. Set to zero to select 10% of input.
  }
}
""")

class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """
  def __enter__(self):
    self._stdout = sys.stdout
    self._stderr = sys.stderr
    sys.stdout = self._stringio_stdout = StringIO()
    sys.stderr = self._stringio_stderr = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio_stdout.getvalue().splitlines())
    sys.stdout = self._stdout
    self.extend(self._stringio_stderr.getvalue().splitlines())
    sys.stderr = self._stderr

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

  if not gs_params.grid_search.flag_on:
    cmd.Command.start("Reading input list from file")
    with open('{}/input_images.lst'\
              ''.format(os.path.abspath(gs_params.output)), 'r') as listfile:
      listfile_contents = listfile.read()
    inp_list = listfile_contents.splitlines()
    cmd.Command.end("Reading input list from file -- DONE")
    return inp_list

  if gs_params.advanced.single_img != None:
    inp_list = [gs_params.advanced.single_img]
    return inp_list

  if gs_params.input_list != None:
    cmd.Command.start("Reading input list from file")
    with open(gs_params.input_list, 'r') as listfile:
      listfile_contents = listfile.read()
    input_list = listfile_contents.splitlines()
    cmd.Command.end("Reading input list from file -- DONE")
  else:
  # search for diffraction image files within the tree and record in a list w/
  # full absolute path and filenames
    cmd.Command.start("Generating input list")
    for root, dirs, files in os.walk(abs_inp_path):
      for filename in files:
        found_file = os.path.join(root, filename)
        if found_file.endswith(('pickle', 'mccd', 'cbf', 'img')):
          input_list.append(found_file)
    cmd.Command.end("Generating input list -- DONE")

    if len(input_list) == 0:
      print "\nERROR: No data found!"
      sys.exit()

  # Pick a randomized subset of images
  if gs_params.advanced.random_sample.flag_on and \
     gs_params.advanced.random_sample.number < len(input_list):
    random_inp_list = []

    if gs_params.advanced.random_sample.number == 0:
      if len(input_list) <= 5:
        random_sample_number = len(input_list)
      elif len(input_list) <= 50:
        random_sample_number = 5
      else:
        random_sample_number = int(len(input_list) * 0.1)
    else:
      random_sample_number = gs_params.advanced.random_sample.number

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

    inp_list = random_inp_list

  else:
    inp_list = input_list

  return inp_list

def make_raw_input(input_list, gs_params):

  raw_input_list = []
  converted_img_list = []

  conv_input_dir = os.path.abspath("{}/conv_pickles".format(os.curdir))

  if os.path.isdir(conv_input_dir):
    shutil.rmtree(conv_input_dir)
  os.makedirs(conv_input_dir)

  for raw_image in input_list:
    path = os.path.dirname(raw_image)
    img_filename = os.path.basename(raw_image)
    filename_no_ext = img_filename.split('.')[0]
    conv_filename = filename_no_ext + "_prep.pickle"

    if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
      dest_folder = conv_input_dir
    else:
      dest_folder = '{0}/{1}'.format(conv_input_dir,
                                     os.path.relpath(path,
                                     os.path.abspath(gs_params.input)))
    converted_img_list.append(os.path.join(dest_folder, conv_filename))

    if not os.path.isdir(dest_folder):
      os.makedirs(dest_folder)

  return converted_img_list, conv_input_dir


def make_dir_lists(input_list, input_folder, output_folder):
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

  abs_inp_path = os.path.abspath(input_folder)
  abs_out_path = os.path.abspath(output_folder)

  # make lists of input and output directories and files
  for input_entry in input_list:
    path = os.path.dirname(input_entry)

    if os.path.relpath(path, abs_inp_path) == '.':  # in case of input in one dir
      input_dir = abs_inp_path
      output_dir = abs_out_path
    else:                                           # in case of input in tree
      input_dir = abs_inp_path + '/' + os.path.relpath(path, abs_inp_path)
      output_dir = abs_out_path + '/' + os.path.relpath(path, abs_inp_path)

    input_dir = os.path.normpath(input_dir)
    output_dir = os.path.normpath(output_dir)

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

    input_folder = os.path.dirname(os.path.commonprefix(input_list))
    output_folder = os.path.abspath(gs_params.output)

    rel_path = os.path.normpath(os.path.relpath(path,
                                os.path.abspath(input_folder)))

    if rel_path == '.':
      output_dir = output_folder
    else:
      output_dir = os.path.normpath(os.path.join(output_folder, rel_path))

    current_output_dir = os.path.normpath(os.path.join(output_dir,
                                              img_filename.split('.')[0]))
    mp_output_entry = [current_output_dir, current_img,
                       img_filename.split('.')[0]]
    mp_output.append(mp_output_entry)

    # Generate spotfinding parameter ranges, make input list w/ filenames
    h_min = gs_range[0]
    h_max = gs_range[1]
    h_avg = gs_params.grid_search.h_avg
    h_std = gs_params.grid_search.h_std
    a_min = gs_range[2]
    a_max = gs_range[3]

    for spot_area in range(a_min, a_max + 1):
      for spot_height in range (h_min, h_max + 1):
        if gs_params.advanced.sig_height_search:
          if spot_height >= 1 + h_std:
            sigs = range(spot_height - h_std, spot_height + 1)
          elif spot_height < 1 + h_std:
            sigs = range(1, spot_height + 1)
          elif spot_height == 1:
            sigs = [1]
          for sig_height in sigs:
            mp_item = [current_img, current_output_dir, sig_height, spot_height,
                       spot_area]
            mp_input.append(mp_item)
        else:
          mp_item = [current_img, current_output_dir, spot_height, spot_height,
                     spot_area]
          mp_input.append(mp_item)

  return mp_input, mp_output

def make_dirs (mp_output_list, gs_params):
  """ Generates output directory tree, which mirrors the input directory tree
  """

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

  for item in os.listdir(os.curdir):
    if item.find('iota') != -1 and item.find('.log') != -1:
      log_count += 1

  if log_count > 0:
    old_log_filename = logfile
    new_log_filename = '{0}/iota_{1}.log'.format(os.curdir, log_count)
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


def generate_input(gs_params, input_list, input_folder):
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
              gs_params.grid_search.a_avg + gs_params.grid_search.a_std]


  # Make log directory and input/output directory lists
  log_dir = "{}/logs".format(os.path.abspath(gs_params.output))
  cmd.Command.start("Reading data folder(s)")
  input_dir_list, output_dir_list = make_dir_lists(input_list, input_folder,
                                                   gs_params.output)
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
  logfile = '{}/iota.log'.format(os.curdir)
  main_log_init(logfile)

  return gs_range, input_dir_list, output_dir_list, log_dir,\
         logfile, mp_input_list, mp_output_list

def write_defaults(current_path, txt_out):
  """ Generates list of default parameters for a reasonable target file
      (target.phil), which will be created in the folder from which IOTA is
      being run. Also writes out the IOTA parameter file.

      input: current_path - absolute path to current folder
             txt_out - IOTA parameters in text format
  """

  def_target_file = '{}/target.phil'.format(current_path)
  default_target = ['# -*- mode: conf -*-',
                    '# target_cell = 79.4 79.4 38.1 90 90 90  # insert your own target cell once known',
                    '# known_setting = 5                      # likewise known setting',
                    '# target_cell_centring_type = P          # likewise centring type',
                    'difflimit_sigma_cutoff = 0.01',
                    'force_method2_resolution_limit = 2.5',
                    'distl_highres_limit = 2.5',
                    'distl_lowres_limit=50.0',
                    'distl{',
                    '  #verbose=True',
                    '  res.outer=2.5',
                    '  res.inner=50.0',
                    '  peak_intensity_maximum_factor=1000',
                    '  spot_area_maximum_factor=20',
                    '  compactness_filter=False',
                    '  method2_cutoff_percentage=2.5',
                    '}',
                    'integration {',
                    '  background_factor=2',
                    '  model=user_supplied',
                    '  spotfinder_subset=spots_non-ice',
                    '  mask_pixel_value=-2',
                    '  detector_gain=0.32',
                    '  greedy_integration_limit=True',
                    '  combine_sym_constraints_and_3D_target=True',
                    '  spot_prediction=dials',
                    '  integration.guard_width_sq=4.',
                    '  mosaic {',
                    '    refinement_target=ML',
                    '    domain_size_lower_limit=4.',
                    '    enable_rotational_target_highsym=False',
                    '  }',
                    '}',
                    'mosaicity_limit=2.0',
                    'dist_permit_binning=False',
                    'distl_minimum_number_spots_for_indexing=16',
                    'distl_permit_binning=False',
                    'beam_search_scope=5'
                    ]
  with open(def_target_file, 'w') as targ:
    for line in default_target:
      targ.write('{}\n'.format(line))

  with open('{}/iota.param'.format(current_path), 'w') as default_settings_file:
    default_settings_file.write(txt_out)

def auto_mode(current_path, data_path, now):
  """ Automatically builds the IOTA parameter file and the target PHIL file and
      begins processing using reasonable default parameters.

      input:  current_path - absolute path to current directory
              data_path - provided path to folder with raw images
              now - current date and time

      output: gs_params - full list of parameters for IOTA
              txt_out - same but in text format suitable for printing
  """

  # Modify list of default IOTA parameters to include the absolute path to data
  # folder and a description with a time-stamp

  cmd.Command.start("Generating default parameters")

  gs_params = master_phil.extract()
  gs_params.description = 'IOTA parameters auto-generated on {}'.format(now)
  gs_params.input = data_path
  #gs_params.advanced.random_sample.flag_on = True
  #gs_params.advanced.random_sample.number = 0

  mod_phil = master_phil.format(python_object=gs_params)

  #capture input read out by phil
  with Capturing() as output:
    mod_phil.show()

  txt_out = ''
  for one_output in output:
    txt_out += one_output + '\n'

  # Write default parameter and target files
  write_defaults(current_path, txt_out)

  cmd.Command.end("Generating default parameters -- DONE")

  return gs_params, txt_out

def print_params():

  #capture input read out by phil
  with Capturing() as output:
    master_phil.show(attributes_level=1)

  help_out = ''
  for one_output in output:
    help_out += one_output + '\n'

  #capture input read out by phil
  with Capturing() as output:
    master_phil.show()

  txt_out = ''
  for one_output in output:
    txt_out += one_output + '\n'


  return help_out, txt_out

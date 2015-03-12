from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 03/06/2015
Description : IOTA I/O module. Reads PHIL input, creates output directories, etc.
'''


import sys
import os
import random

import iotbx.phil

master_phil = iotbx.phil.parse("""
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
n_processors = 32
  .type = int
  .help = No. of processing units
""")

def process_input(input_file_list):

  user_phil = []
  for input_file in input_file_list:
    user_phil.append(iotbx.phil.parse(open(input_file).read()))

  working_phil = master_phil.fetch(sources=user_phil)
  params = working_phil.extract()

  #capture input read out by phil
  from cStringIO import StringIO
  class Capturing(list):
    def __enter__(self):
      self._stdout = sys.stdout
      sys.stdout = self._stringio = StringIO()
      return self
    def __exit__(self, *args):
      self.extend(self._stringio.getvalue().splitlines())
      sys.stdout = self._stdout

  with Capturing() as output:
    working_phil.show()

  txt_out = ''
  for one_output in output:
    txt_out += one_output + '\n'

  return params, txt_out

# Read input directory tree (if any) and make lists of input folder, output folder and
# input files for use in everything
def make_input_list (gs_params):
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
    for i in range(gs_params.random_sample.number):
      random_number = random.randrange(0, len(input_list))
      random_inp_list.append(input_list[random_number])

    inp_list = random_inp_list
  else:
    inp_list = input_list

  return inp_list

def make_dir_lists(input_list, gs_params):

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

    #with open('{}/input_files.lst'.format(abs_out_path), 'a') as inp_list_file:
    #  inp_list_file.write('{0}, {1}\n'.format(input_entry, output_dir))

  return input_dir_list, output_dir_list

# Generates input list for MP grid seach
def make_mp_input(input_list, gs_params, gs_range):

  mp_item = []
  mp_input = []
  mp_output = []

  if gs_params.advanced.random_sample.flag_on == True:
    print "Selecting {0} samples from {1} images in {2}:"\
          "".format(gs_params.random_sample.number, len(input_list),
                    gs_params.input)
    random_inp_list = []
    for i in range(gs_params.random_sample.number):
      random_number = random.randrange(0, len(input_list))
      print input_list[random_number]
      random_inp_list.append(input_list[random_number])

    gs_params.grid_search.flag_on = True
    gs_params.grid_search.h_avg = 5
    gs_params.grid_search.h_std = 4
    gs_params.grid_search.a_avg = 5
    gs_params.grid_search.a_std = 4

    inp_list = random_inp_list
  else:
    inp_list = input_list


  for current_img in inp_list:
    # generate filenames, etc.
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

    h_min = gs_range[0]
    h_max = gs_range[1]
    a_min = gs_range[2]
    a_max = gs_range[3]

    # Create input list w/ filename and spot-finding params
    for sig_height in range(h_min, h_max + 1):
      for spot_area in range (a_min, a_max + 1):
        mp_item = [current_img, sig_height, sig_height, spot_area,
                   current_output_dir]
        mp_input.append(mp_item)

  return mp_input, mp_output

# Make output directories preserving the tree structure
def make_dirs (input_list, gs_params):

  for current_img in input_list:
    # generate filenames, etc.
    path = os.path.dirname(current_img)
    img_filename = os.path.basename(current_img)

    if os.path.relpath(path, os.path.abspath(gs_params.input)) == '.':
      output_dir = os.path.abspath(gs_params.output)
    else:
      output_dir = os.path.normpath('{0}/{1}'\
                                    ''.format(os.path.abspath(gs_params.output),
                                              os.path.relpath(path,
                                              os.path.abspath(gs_params.input))))

    current_output_dir = os.path.normpath("{0}/tmp_{1}".format(output_dir,
                                              img_filename.split('.')[0]))

#     # Make directories for output / log file for the image being integrated
    if not os.path.exists(current_output_dir):
      os.makedirs(current_output_dir)

# Save log from previous run and initiate a new log (run once)
# This is only necessary if re-running selection after grid-search
# Will likely be deprecated soon
def main_log_init(logfile):

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


# Write main log (so that I don't have to repeat this every time)
def main_log(logfile, entry):

  with open(logfile, 'a') as logfile:
    logfile.write('{}\n'.format(entry))


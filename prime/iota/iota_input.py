from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 01/26/2015
Description : IOTA I/O module. Reads PHIL input, creates output directories, etc.
'''


import sys
import os
import shutil
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
flag_inp_test = False
  .type = bool
  .help = Test of input list use for selection
pred_img = None
  .type = str
  .help = Output a PNG image of integration predictions overlayed on the raw image
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
grid_search
  .help = "Parameters for the grid search."
{
  flag_on = True
    .type = bool
    .help = Set to False to run selection of pickles only.
  a_min = 1
    .type = int
    .help = Minimum spot area.
  a_max = 10
    .type = int
    .help = Maximum spot area.
  h_min = 1
    .type = int
    .help = Minimum spot height.
  h_max = 10
    .type = int
    .help = Maximum spot height.
}
selection_res_limit
  .help = "Resolution limit for pickle selection."
{
  d_min = 6.0
    .type = float
    .help = Highest resolution limit.
  d_max = 15.0
    .type = float
    .help = Lowest resolution limit.
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
select_by = strong
        .type = str
  .help = Pickle selection method
min_sigma = 5.0
  .type = float
  .help = Minimum sigma for "strong" reflections
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
def make_input_list (input_dir):
  input_list = []
  abs_inp_path = os.path.abspath(input_dir)

  # search for *.pickle files within the tree and record in a list w/
  # full absolute path and filanames
  for root, dirs, files in os.walk(abs_inp_path):
    for filename in files:
      if filename.endswith("pickle"):
        pickle_file = os.path.join(root, filename)
        input_list.append(pickle_file)

  return input_list

def make_dir_lists(input_list, input_dir, output_dir):

  input_dir_list = []
  output_dir_list = []

  abs_inp_path = os.path.abspath(input_dir)
  abs_out_path = os.path.abspath(output_dir)

  # make lists of input and output directories and files
  for input_entry in input_list:
    path = os.path.dirname(input_entry)

    if os.path.relpath(path, abs_inp_path) == '.':  # in case of input in one dir
      input_dir = abs_inp_path
      if input_dir not in input_dir_list: input_dir_list.append(input_dir)
      output_dir = abs_out_path
    else:                                           # in case of input in tree
      input_dir = abs_inp_path + '/' + os.path.relpath(path, abs_inp_path)
      output_dir = abs_out_path + '/' + os.path.relpath(path, abs_inp_path)

    if input_dir not in input_dir_list: input_dir_list.append(input_dir)
    if output_dir not in output_dir_list: output_dir_list.append(output_dir)

    #with open('{}/input_files.lst'.format(abs_out_path), 'a') as inp_list_file:
    #  inp_list_file.write('{0}, {1}\n'.format(input_entry, output_dir))

  return input_dir_list, output_dir_list, abs_out_path + '/logs'

# Generates input list for MP grid seach
def make_mp_input(input_list, log_dir, gs_params):

  mp_item = []
  mp_input = []
  mp_output = []

  if gs_params.random_sample.flag_on == True:
    print "Selecting {0} samples from {1} images in {2}:".format(gs_params.random_sample.number, len(input_list), gs_params.input)
    random_inp_list = []
    for i in range(gs_params.random_sample.number):
      random_number = random.randrange(0, len(input_list))
      print input_list[random_number]
      random_inp_list.append(input_list[random_number])

    gs_params.grid_search.flag_on = True
    gs_params.grid_search.h_min = 1
    gs_params.grid_search.h_max = 10
    gs_params.grid_search.a_min = 1
    gs_params.grid_search.a_max = 10

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
    index_log_dir = "{0}/tmp_{1}".format(log_dir, img_filename.split('.')[0])
    mp_output_entry = [current_img, current_output_dir]
    mp_output.append(mp_output_entry)

    # Create input list w/ filename and spot-finding params
    for sig_height in range(gs_params.grid_search.h_min,
                          gs_params.grid_search.h_max + 1):
      for spot_area in range (gs_params.grid_search.a_min,
                              gs_params.grid_search.a_max + 1):
        mp_item = [current_img, sig_height, sig_height, spot_area]
        mp_input.append(mp_item)

  return mp_input, mp_output

# Make output directories preserving the tree structure
def make_dirs (input_list, output_dir_list, log_dir, gs_params):

  mp_item = []
  mp_input = []
  mp_output = []

  # Make output directory structure
  for output_dir in output_dir_list:
    if os.path.exists(output_dir):
      shutil.rmtree(output_dir)
      os.makedirs(output_dir)
    else:
      os.makedirs(output_dir)

  # make log folder (under main output folder regardless of tree structure)
  if os.path.exists(log_dir):
    shutil.rmtree(log_dir)
    os.makedirs(log_dir)
  else:
    os.makedirs(log_dir)

  for current_img in input_list:
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
    index_log_dir = "{0}/tmp_{1}".format(log_dir, img_filename.split('.')[0])

    # Make directories for output / log file for the image being integrated
    if not os.path.exists(current_output_dir):
      os.makedirs(current_output_dir)
    if not os.path.exists(index_log_dir):
      os.makedirs(index_log_dir)

from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 09/01/2015
Description : IOTA I/O module. Reads PHIL input, also creates reasonable IOTA
              and PHIL defaults if selected.
'''


import sys
import os
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
  .multiple = True
  .help = Path to folder with raw data in pickle format, list of files or single file
  .help = Can be a tree with folders
  .optional = False
target = target.phil
  .type = str
  .multiple = False
  .help = Target (.phil) file with integration parameters
  .optional = False
image_conversion
  .help = Parameters for raw image conversion to pickle format
{
  rename_pickle_prefix = Auto
    .type = str
    .help = Specify prefix (e.g. "HEWL_room_temp") to rename all input images
    .help = Set to None to keep original image filenames and directory tree
  convert_only = False
    .type = bool
    .help = Set to True (or use -c option) to convert and exit
  square_mode = None *pad crop
    .type = choice
    .help = Method to generate square image
  beamstop = 0
    .type = float
    .help = Beamstop shadow threshold, zero to skip
  beam_center
    .help = Alternate beam center coordinates (set to zero to leave the same)
  {
    x = 0
      .type = float
    y = 0
      .type = float
  }
}
image_triage
  .help = Check if images have diffraction using basic spotfinding (-t option)
{
  flag_on = True
    .type = bool
    .help = Set to true to activate
  min_Bragg_peaks = 10
    .type = int
    .help = Minimum number of Bragg peaks to establish diffraction
}
grid_search
  .help = "Parameters for the grid search."
{
  flag_on = True
    .type = bool
    .help = Set to false to only use median spotfinding parameters
  area_median = 5
    .type = int
    .help = Median spot area.
  area_range = 2
    .type = int
    .help = Plus/minus range for spot area.
  height_median = 4
    .type = int
    .help = Median spot height.
  height_range = 2
    .type = int
    .help = Plus/minus range for spot height.
  sig_height_search = False
    .type = bool
    .help = Set to true to scan signal height in addition to spot height
}
selection
  .help = Parameters for integration result selection
{
  select_only
    .help = set to True to re-do selection with previous
    {
    flag_on = False
      .type = bool
      .help = set to True to bypass grid search and just run selection
    grid_search_path = None
      .type = path
      .help = set if you want to use specific grid_search results
      .help = leave as None to use grid search results from previous run
    }
  min_sigma = 5
    .type = int
    .help = minimum I/sigma(I) cutoff for "strong spots"
  select_by = *mosaicity epv
    .type = choice
    .help = Use mosaicity or Ewald proximal volume for optimal parameter selection
  prefilter
    .help = Used to throw out integration results that do not fit user-defined unit cell information
  {
    flag_on = False
      .type = bool
      .help = Set to True to activate prefilter
    target_pointgroup = None
      .type = str
      .help = Target point group, e.g. "P4"
    target_unit_cell = None
      .type = unit_cell
      .help = In format of "a, b, c, alpha, beta, gamma", e.g. 79.4, 79.4, 38.1, 90.0, 90.0, 90.0
    target_uc_tolerance = 0.05
      .type = float
      .help = Maximum allowed unit cell deviation from target
    min_reflections = 0
      .type = int
      .help = Minimum integrated reflections per image
    min_resolution = None
      .type = float
      .help = Minimum resolution for accepted images
  }
}
analysis
  .help = "Analysis / visualization options."
{
  cluster_threshold = 5000
    .type = int
    .help = threshold value for unit cell clustering
  viz = *None integration cv_vectors
    .type = choice
    .help = Set to "integration" to visualize spotfinding and integration results.
    .help = Set to "cv_vectors" to visualize accuracy of CV vectors
  charts = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity, rmsd, etc.
  heatmap = None show *file both
    .type = choice
    .help = Show / output to file a heatmap of grid search results
}
advanced
  .help = "Advanced, debugging and experimental options."
{
  integrate_with = *cctbx dials
    .type = choice
    .help = Choose image processing software package
  debug = False
    .type = bool
    .help = Used for various debugging purposes.
  experimental = False
    .type = bool
    .help = Set to true to run the experimental section of codes
  random_sample
    .help = Use a randomized subset of images (or -r <number> option)
  {
    flag_on = False
      .type = bool
      .help = Set to run grid search on a random set of images.
    number = 0
      .type = int
      .help = Number of random samples. Set to zero to select 10% of input.
  }
}
n_processors = 32
  .type = int
  .help = No. of processing units
mp_method = *multiprocessing mpi
  .type = choice
  .help = Multiprocessing method
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

def process_input(args,
                  phil_args,
                  input_file,
                  mode='auto',
                  now=None):
  """ Read and parse parameter file

      input: input_file_list - PHIL-format files w/ parameters

      output: params - PHIL-formatted parameters
              txt_output - plain text-formatted parameters
  """

  from libtbx.phil.command_line import argument_interpreter
  from libtbx.utils import Sorry

  cmd.Command.start("Generating parameters")
  if mode == 'file':
    user_phil = [ip.parse(open(inp).read()) for inp in [input_file]]
    working_phil = master_phil.fetch(sources=user_phil)
    params = working_phil.extract()
  elif mode == 'auto':
    params = master_phil.extract()
    params.description = 'IOTA parameters auto-generated on {}'.format(now)
    params.input = [input_file]

  # Check for -r option and set random subset parameter
  if args.random > 0:
    params.advanced.random_sample.flag_on = True
    params.advanced.random_sample.number = args.random[0]

  # Check for -n option and set number of processors override
  # (for parallel map only, for now)
  if args.nproc > 0:
    params.n_processors = args.nproc[0]

  # Check for -c option and set flags to exit IOTA after raw image conversion
  if args.convert:
    params.image_conversion.convert_only = True

  # Check -p option to see if converted file prefix is supplied; will run
  # conversion automatically if prefix is supplied
  if str(args.prefix).lower() != "auto":
    params.image_conversion.convert_images = True
    params.image_conversion.rename_pickle_prefix = args.prefix

  #Check -s option to bypass grid search and run selection/integration only
  if args.select:
    params.selection.select_only.flag_on = True

  final_phil = master_phil.format(python_object=params)

  argument_interpreter = argument_interpreter(master_phil=master_phil)
  consume = []
  for arg in phil_args:
    try:
      command_line_params = argument_interpreter.process(arg=arg)
      final_phil = final_phil.fetch(sources=[command_line_params,])
      consume.append(arg)
    except Sorry,e:
      pass
  for item in consume:
    phil_args.remove(item)
  if len(phil_args) > 0:
    raise Sorry("Not all arguments processed, remaining: {}".format(phil_args))

  temp_phil = [final_phil]
  params = final_phil.extract()
  diff_phil = master_phil.fetch_diff(sources=temp_phil)

  with Capturing() as output:
    diff_phil.show()
  diff_out = ''
  for one_output in output:
    diff_out += one_output + '\n'

  if mode == 'auto':
    with Capturing() as diff_output:
      final_phil.show()
    txt_out = ''
    for one_output in diff_output:
      txt_out += one_output + '\n'
    write_defaults(os.path.abspath(os.curdir), txt_out)

  cmd.Command.end("Generating parameters -- DONE")

  return params, diff_out


def write_defaults(current_path, txt_out):
  """ Generates list of default parameters for a reasonable target file
      (target.phil), which will be created in the folder from which IOTA is
      being run. Also writes out the IOTA parameter file.

      input: current_path - absolute path to current folder
             txt_out - IOTA parameters in text format
  """

  def_target_file = '{}/target.phil'.format(current_path)
  default_target = ['# -*- mode: conf -*-',
                    '# target_cell = 79.4 79.4 38.1 90 90 90  # insert your own target unit cell if known',
                    '# known_setting = 9                      # Triclinic = 1, monoclinic = 2,',
                    '                                         # orthorhombic/rhombohedral = 5, tetragonal = 9,',
                    '                                         # hexagonal = 12, cubic = 22,',
                    '# target_cell_centring_type = *P C I R F',
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
                    '  guard_width_sq=4.',
                    '  mosaic {',
                    '    refinement_target=ML',
                    '    domain_size_lower_limit=4.',
                    '    enable_rotational_target_highsym=False',
                    '  }',
                    '}',
                    'mosaicity_limit=2.0',
                    'distl_minimum_number_spots_for_indexing=16',
                    'distl_permit_binning=False',
                    'beam_search_scope=5'
                    ]
  with open(def_target_file, 'w') as targ:
    for line in default_target:
      targ.write('{}\n'.format(line))

  with open('{}/iota.param'.format(current_path), 'w') as default_settings_file:
    default_settings_file.write(txt_out)

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

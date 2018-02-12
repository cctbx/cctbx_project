from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/12/2018
Description : IOTA I/O module. Reads PHIL input, also creates reasonable IOTA
              and PHIL defaults if selected.
'''


import sys
import os
from cStringIO import StringIO

import iotbx.phil as ip

master_phil = ip.parse("""
description = None
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
output = None
  .type = path
  .multiple = False
  .help = Base output directory, current directory in command-line, can be set in GUI
  .optional = True
image_conversion
  .help = Parameters for raw image conversion to pickle format
{
  rename_pickle = None keep_file_structure *auto_filename custom_filename
    .type = choice
    .help = "keep_file_structure" retains the input filenames w/ folder tree
    .help = "auto_filename" is <your_login_name>_<conversion_run>_<#####>
  rename_pickle_prefix = None
    .type = str
    .help = Will only be used in conjunction with "custom_filename"
  convert_only = False
    .type = bool
    .help = Set to True (or use -c option) to convert and exit
  square_mode = None no_modification *pad crop
    .type = choice
    .help = Method to generate square image
  mask = None
    .type = path
    .help = Mask for ignored pixels
  invert_boolean_mask = False
    .type = bool
    .help = Set to true if applying DIALS mask to Labelit processing
  beamstop = 0
    .type = float
    .help = Beamstop shadow threshold, zero to skip
  distance = 0
    .type = float
    .help = Alternate crystal-to-detector distance (set to zero to leave the same)
  beam_center
    .help = Alternate beam center coordinates (in PIXELS)
    .help = Set to zero to leave the same
  {
    x = 0
      .type = int
    y = 0
      .type = int
  }
}
image_triage
  .help = Check if images have diffraction using basic spotfinding (-t option)
{
  type = None no_triage *simple grid_search
    .type = choice
    .help = Set to None to attempt integrating all images
  min_Bragg_peaks = 10
    .type = int
    .help = Minimum number of Bragg peaks to establish diffraction
  grid_search
    .help = "Parameters for the grid search."
  {
    area_min = 6
      .type = int
      .help = Minimal spot area.
    area_max = 24
      .type = int
      .help = Maximal spot area.
    height_min = 2
      .type = int
      .help = Minimal spot height.
    height_max = 20
      .type = int
      .help = Maximal spot height.
    step_size = 4
      .type = int
      .help = Grid search step size
  }
}
cctbx
  .help = Options for CCTBX-based image processing
{
  target = None
    .type = str
    .multiple = False
    .help = Target (.phil) file with integration parameters
  resolution_limits
    .help = Sets several resolution limits in Labelit settings
  {
    low = 50.0
      .type = float
    high = 1.5
      .type = float
  }
  target_lattice_type = *None triclinic monoclinic orthorhombic tetragonal rhombohedral hexagonal cubic
    .type = choice
    .help = Target Bravais lattice type if known
  target_centering_type = *None P C I R F
    .type = choice
    .help = Target lattice centering type if known
  target_unit_cell = None
    .type = unit_cell
    .help = Target unit cell parameters (if known)
  grid_search
    .help = "Parameters for the grid search."
  {
    type = None no_grid_search *brute_force smart
      .type = choice
      .help = Set to None to only use median spotfinding parameters
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
    select_by = *epv mosaicity
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
      target_uc_tolerance = None
        .type = float
        .help = Maximum allowed unit cell deviation from target
      min_reflections = None
        .type = int
        .help = Minimum integrated reflections per image
      min_resolution = None
        .type = float
        .help = Minimum resolution for accepted images
    }
  }
}
dials
  .help = Options for DIALS-based image processing
  .help = This option is not yet ready for general use!
{
  target = None
    .type = str
    .multiple = False
    .help = Target (.phil) file with integration parameters for DIALS
  target_space_group = None
    .type = space_group
    .help = Target space (or point) group (if known)
  target_unit_cell = None
    .type = unit_cell
    .help = Target unit cell parameters (if known)
  use_fft3d = False
    .type = bool
    .help = Set to True to use FFT3D in indexing
  significance_filter
    .help = Set to True and add value to determine resolution based on I/sigI
  {
    flag_on = True
      .type = bool
      .help = Set to true to activate significance filter
    sigma = 1.0
      .type = float
      .help = Sigma level to determine resolution cutoff
  }
  determine_sg_and_reindex = True
    .type = bool
    .help = Will determine sg and reindex if no target space group supplied
  auto_threshold = True
    .type = bool
    .help = Set to True to estimate global threshold for each image
  filter
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
      target_uc_tolerance = None
        .type = float
        .help = Maximum allowed unit cell deviation from target
      min_reflections = None
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
  run_clustering = False
    .type = bool
    .help = Set to True to turn on hierarchical clustering of unit cells
  cluster_threshold = 5000
    .type = int
    .help = threshold value for unit cell clustering
  viz = None *no_visualization integration cv_vectors
    .type = choice
    .help = Set to "integration" to visualize spotfinding and integration results.
    .help = Set to "cv_vectors" to visualize accuracy of CV vectors
  charts = False
    .type = bool
    .help = If True, outputs PDF files w/ charts of mosaicity, rmsd, etc.
  summary_graphs = False
    .type = bool
    .help = If True: spot-finding heatmap, res. histogram and beamXY graph
}
advanced
  .help = "Advanced, debugging and experimental options."
{
  integrate_with = *cctbx dials
    .type = choice
    .help = Choose image processing software package
  estimate_gain = True
    .type = bool
    .help = Estimates detector gain (helps improve spotfinding in DIALS)
  flip_beamXY = True
    .type = bool
    .help = flip beamX and beamY parameters when modifying image for cctbx
  debug = False
    .type = bool
    .help = Used for various debugging purposes.
  experimental = False
    .type = bool
    .help = Set to true to run the experimental section of codes
  image_viewer = *dials.image_viewer cctbx.image_viewer distl.image_viewer cxi.view
    .type = choice
    .help = Select image viewer (GUI only)
  monitor_mode = False
    .type = bool
    .help = Set to true to keep watch for incoming images (GUI only)
  monitor_mode_timeout = False
    .type = bool
    .help = Set to true to auto-terminate continuous mode (GUI only)
  monitor_mode_timeout_length = 0
    .type = int
    .help = Timeout length in seconds (GUI only)
  prime_prefix = prime
    .type = str
    .help = Prefix for the PRIME script filename
  temporary_output_folder = None
    .type = path
    .help = If None, temp output goes to <output>/integration/###/tmp/
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
mp_method = *multiprocessing mpi lsf torq
  .type = choice
  .help = Multiprocessing method
mp_queue = None
  .type = str
  .help = Multiprocessing queue
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

  if mode == 'file':
    user_phil = [ip.parse(open(inp).read()) for inp in [input_file]]
    working_phil = master_phil.fetch(sources=user_phil)
    params = working_phil.extract()
  elif mode == 'auto':
    params = master_phil.extract()
    params.description = 'IOTA parameters auto-generated on {}'.format(now)
    params.input = [input_file]

  final_phil = master_phil.format(python_object=params)

  # Parse in-line params into phil
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

  # Perform command line check and modify params accordingly
  params = final_phil.extract()

  if mode == 'auto':
    output_dir = os.path.abspath(os.curdir)
    if params.advanced.integrate_with == 'dials':
      params.dials.target = os.path.join(output_dir, 'dials.phil')
    elif params.advanced.integrate_with == 'cctbx':
      params.cctbx.target = os.path.join(output_dir, 'cctbx.phil')

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
    params.cctbx.selection.select_only.flag_on = True

  # Check if grid search is turned off; if so, set everything to zero
  if str(params.cctbx.grid_search.type).lower() == 'none':
    params.cctbx.grid_search.area_range = 0
    params.cctbx.grid_search.height_range = 0
    params.cctbx.grid_search.sig_height_search = False

  final_phil = master_phil.format(python_object=params)

  temp_phil = [final_phil]
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
    write_defaults(os.path.abspath(os.curdir), txt_out, params.advanced.integrate_with)

  return params, diff_out


def write_defaults(current_path=None, txt_out=None, method='cctbx',
                   write_target_file=True, write_param_file=True):
  """ Generates list of default parameters for a reasonable target file:
        - if cctbx.xfel, target.phil will be created in the folder from which IOTA is
          being run.
        - if DIALS, dials.phil will be created in the folder from which IOTA is being run.
        - also writes out the IOTA parameter file.

      input: current_path - absolute path to current folder
             txt_out - IOTA parameters in text format
  """

  if method == 'cctbx':
    def_target_file = '{}/cctbx.phil'.format(current_path)
    default_target = ['# -*- mode: conf -*-',
                      'difflimit_sigma_cutoff = 0.01',
                      'distl{',
                      '  peak_intensity_maximum_factor=1000',
                      '  spot_area_maximum_factor=20',
                      '  compactness_filter=False',
                      '  method2_cutoff_percentage=2.5',
                      '}',
                      'integration {',
                      '  background_factor=2',
                      '  enable_one_to_one_safeguard=True',
                      '  model=user_supplied',
                      '  spotfinder_subset=spots_non-ice',
                      '  mask_pixel_value=-2',
                      '  greedy_integration_limit=True',
                      '  combine_sym_constraints_and_3D_target=True',
                      '  spot_prediction=dials',
                      '  guard_width_sq=4.',
                      '  mosaic {',
                      '    refinement_target=LSQ *ML',
                      '    domain_size_lower_limit=4.',
                      '    enable_rotational_target_highsym=True',
                      '  }',
                      '}',
                      'mosaicity_limit=2.0',
                      'distl_minimum_number_spots_for_indexing=16',
                      'distl_permit_binning=False',
                      'beam_search_scope=0.5'
                      ]
  elif method == 'dials':
    def_target_file = '{}/dials.phil'.format(current_path)
    default_target = ['verbosity=10',
                      'spotfinder {',
                      '  threshold {',
                      '    dispersion {',
                      '      gain = 1',
                      '      sigma_strong = 3',
                      '      global_threshold = 0',
                      '      }',
                      '   }',
                      '}',
                      'geometry {',
                      '  detector {',
                      '    distance = None',
                      '    slow_fast_beam_centre = None',
                      '  }',
                      '}',
                      'indexing {',
                      '  refinement_protocol {',
                      '    n_macro_cycles = 1',
                      '    d_min_start = 2.0',
                      '  }',
                      '  basis_vector_combinations.max_combinations = 10',
                      '  stills { ',
                      '    indexer = stills',
                      '    method_list = fft1d real_space_grid_search',
                      '  }',
                      '}',
                      'refinement {',
                      '  parameterisation {',
                      '    beam.fix = *all in_spindle_plane out_spindle_plane wavelength',
                      '    detector  {',
                      '      fix = all position orientation',
                      '      hierarchy_level=0',
                      '    }',
                      '    auto_reduction {',
                      '      action=fix',
                      '      min_nref_per_parameter=1',
                      '    }',
                      '  }',
                      '  reflections {',
                      '    outlier.algorithm=null',
                      '    weighting_strategy  {',
                      '      override=stills',
                      '      delpsi_constant=1000000',
                      '    }',
                      '  }',
                      '}',
                      'integration {',
                      '  integrator=stills',
                      '  profile.fitting=False',
                      '  background {',
                      '    simple {',
                      '      outlier {',
                      '        algorithm = null',
                      '      }',
                      '    }',
                      '  }',
                      '}',
                      'profile {',
                      '  gaussian_rs {',
                      '    min_spots.overall = 0',
                      '  }',
                      '}'
                      ]

  if write_target_file:
    with open(def_target_file, 'w') as targ:
      for line in default_target:
        targ.write('{}\n'.format(line))

  if write_param_file:
    with open('{}/iota.param'.format(current_path), 'w') as default_param_file:
      default_param_file.write(txt_out)

  return default_target, txt_out

def print_params():
  """ Print out master parameters for IOTA """

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

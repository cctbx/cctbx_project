from __future__ import division, print_function, absolute_import

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 12/10/2018
Description : IOTA I/O module. Reads PHIL input, also creates reasonable IOTA
              and PHIL defaults if selected.
'''


import os
import iotbx.phil as ip

from iota.components.iota_utils import convert_phil_to_text
from iota.components.iota_processing import cctbx_str

master_phil_str = """
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
output = $PWD
  .type = path
  .multiple = False
  .help = Base output directory
image_import
  .help = Parameters for image importing (& triage)
{
  image_triage = True
    .type = bool
    .help = Perform a round of spotfinding to see if image has diffraction
  minimum_Bragg_peaks = 10
    .type = int
    .help = Minimum number of Bragg peaks to establish diffraction
  strong_sigma = 5.0
    .type = float
    .help = Sigma level to define a strong reflection
  mask = None
    .type = path
    .help = Mask for ignored pixels
  estimate_gain = False
    .type = bool
    .help = Estimates detector gain (helps improve spotfinding)
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
analysis
  .help = "Analysis / visualization options."
{
  run_clustering = False
    .type = bool
    .help = Set to True to turn on hierarchical clustering of unit cells
  cluster_threshold = 5000
    .type = int
    .help = threshold value for unit cell clustering
  cluster_n_images = None
    .type = int
    .help = How many images to cluster? (Useful for huge datasets.)
  cluster_limit = 5
    .type = int
    .help = "Major clusters" are defined as clusters with over n members
  cluster_write_files = True
    .type = bool
    .help = Set to True to write lists of images belonging to major clusters
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
  processing_backend = *cctbx.xfel ha14
    .type = choice
    .help = Choose image processing backend software
  debug = False
    .type = bool
    .help = Used for various debugging purposes.
  experimental = False
    .type = bool
    .help = Set to true to run the experimental section of codes
  prime_prefix = prime
    .type = str
    .help = Prefix for the PRIME script filename
  temporary_output_folder = None
    .type = path
    .help = If None, temp output goes to <output>/integration/###/tmp/
  image_range
    .help = Use a range of images, e.g. 5 - 1000
  {
    flag_on = False
      .type = bool
    range = None
      .type = str
      .help = Can input multiple ranges, e.g. "5, 60-100, 200-1500"
  }
  random_sample
    .help = Use a randomized subset of images (or -r <number> option)
  {
    flag_on = False
      .type = bool
    number = 0
      .type = int
      .help = Number of images in random sample
  }
}
mp
  .help = Multiprocessing options
{
  n_processors = 0
    .type = int
    .help = No. of processing units
  method = *multiprocessing mpi lsf torq custom
    .type = choice
    .help = Multiprocessing method
  queue = None
    .type = str
    .help = Multiprocessing queue
  submit_command = None
    .type = str
    .help = Command to submit IOTA job to a queue
  kill_command = None
    .type = str
    .help = Command to kill the current IOTA job submitted to a queue
}
"""

master_phil = ip.parse(master_phil_str + cctbx_str, process_includes=True)

def get_input_phil(input_source=None, phil_args=None, mode='auto',
                   ha14=False, gui=False):
  ''' Generate PHIL from file, master, and/or command arguments
  :param args: command line arguments
  :param phil_args: PHIL settings as command line arguments
  :param input_source: file with input settings in PHIL format
  :param mode:
  :return:
  '''
  from libtbx.phil.command_line import argument_interpreter
  from libtbx.utils import Sorry

  # Depending on mode, either read input from file, or generate defaults
  if mode == 'file':
    with open(input_source, 'r') as inpf:
      user_phil = ip.parse(inpf.read())
    phil_fixer = PHILFixer()
    working_phil = phil_fixer.run(old_phil=user_phil, write_file=True)
  else:
    if ha14:
      from iota.components.iota_cctbx_ha14 import ha14_str
      working_phil = ip.parse(master_phil_str + ha14_str, process_includes=True)
    else:
      working_phil = master_phil

  if gui:
    from libtbx.phil import find_scope
    if not find_scope(working_phil, 'gui'):
      from iota.components.iota_ui_base import gui_phil
      working_phil.adopt_scope(gui_phil)

  # Parse in-line params into phil
  if phil_args:
    argument_interpreter = argument_interpreter(master_phil=working_phil)
    for arg in phil_args:
      try:
        command_line_params = argument_interpreter.process(arg=arg)
        working_phil = working_phil.fetch(sources=[command_line_params,])
        phil_args.remove(arg)
      except Sorry:
        pass
    if len(phil_args) > 0:
      raise Sorry("Not all arguments processed, remaining: {}".format(phil_args))

  # Self-fetch to resolve variables
  return working_phil.fetch(source=working_phil)

def process_ui_input(args, phil_args, input_source, mode='auto'):
  ''' Read and parse parameter file and/or command-line args for IOTA GUI

  :param args: command-line arguments upon launch
  :param phil_args: command-line arguments pertaining to IOTA parameters
  :param input_source: text file with IOTA parameters
  :return:
  '''
  working_phil = get_input_phil(phil_args=phil_args, input_source=input_source,
                                mode=mode, gui=True)
  params = working_phil.extract()

  # Check for -r option and set random subset parameter
  if args.random > 0:
    params.advanced.random_sample.flag_on = True
    params.advanced.random_sample.number = args.random[0]

  if args.watch > 0:
    params.gui.monitor_mode = True
    params.gui.monitor_mode_timeout = True
    params.gui.monitor_mode_timeout_length = args.watch[0]

  if args.tmp is not None:
    params.advanced.temporary_output_folder = args.tmp[0]

  # Check for -n option and set number of processors override
  # (for parallel map only, for now)
  from multiprocessing import cpu_count
  max_proc = cpu_count() - 2
  if args.nproc > 0:
    if args.nproc >= max_proc:
      params.mp.n_processors = max_proc
    else:
      params.mp.n_processors = args.nproc[0]
  elif params.mp.method == 'multiprocessing':
    if (params.mp.n_processors >= max_proc or
            params.mp.n_processors == 0):
      params.mp.n_processors = int(max_proc / 2)

  return working_phil.format(python_object=params)

def process_input(args, phil_args, input_source, mode='auto', now=None):
  """ Read and parse parameter file and/or command-line args; if none found,
  create a default parameter object

  :param args: command-line arguments upon launch
  :param phil_args: command-line arguments pertaining to IOTA parameters
  :param input_source: text file with IOTA parameters (if 'file' mode) or
                       source of images (if 'auto' mode)
  :param mode: Mode of XtermIOTA run. See the InitAll base class
  :param now: date / time stamp
  :return: PHIL-formatted parameters
  """

  working_phil = get_input_phil(phil_args=phil_args, ha14=args.ha14,
                                mode=mode, input_source=input_source)

  # Perform command line check and modify params accordingly
  params = working_phil.extract()
  if args.ha14:
    params.advanced.processing_backend = 'ha14'

  if mode == 'auto':
    output_dir = os.path.abspath(os.curdir)
    target_filename = '{}.phil'.format(
      params.advanced.processing_backend.replace('.', '_'))
    params.cctbx_xfel.target = os.path.join(output_dir, target_filename)
    params.description = 'IOTA parameters auto-generated on {}'.format(now)
    params.input = [input_source]
    params.output = output_dir

  # Check for -r option and set random subset parameter
  if args.random > 0:
    params.advanced.random_sample.flag_on = True
    params.advanced.random_sample.number = args.random[0]

  # Check for -n option and set number of processors override
  # (for parallel map only, for now)
  from multiprocessing import cpu_count
  max_proc = cpu_count() - 2
  if args.nproc > 0:
    if args.nproc >= max_proc:
      params.mp.n_processors = max_proc
    else:
      params.mp.n_processors = args.nproc[0]
  elif params.mp.method == 'multiprocessing':
    if (params.mp.n_processors >= max_proc or
            params.mp.n_processors == 0):
      params.mp.n_processors = int(max_proc / 2)

  final_phil = working_phil.format(python_object=params)

  if mode == 'auto':
    write_defaults(os.path.abspath(os.curdir), final_phil.as_str(),
                   params.advanced.processing_backend)

  return final_phil


def write_defaults(current_path=None, txt_out=None, method='cctbx_xfel',
                   write_target_file=True, write_param_file=True,
                   filepath=None):
  """ Generates list of default parameters for a reasonable target file:
        - old cctbx.xfel (HA14) now deprecated, but can still be used
        - also writes out the IOTA parameter file.

  :param current_path: path to current output folder
  :param txt_out: text of IOTA parameters
  :param method: which backend is used (default = cctbx.xfel AB18)
  :param write_target_file: write backend parameters to file
  :param write_param_file: write IOTA parameters to file
  :return: default backend settings as list, IOTA parameters as string
  """

  if filepath:
    def_target_file = filepath
  else:
    def_target_file = '{}/{}.phil' \
                      ''.format(current_path, method.replace('.', '_'))

  if method.lower() in ('cctbx', 'ha14'):
    # This is a deprecated backend from 2014; no longer recommended, may suck
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
  elif method.lower() in ('dials', 'cctbx.xfel'):
    default_target_str = '''verbosity=10
                            spotfinder {
                              threshold {
                                dispersion {
                                  gain = 1
                                  sigma_strong = 3
                                  global_threshold = 0
                                  }
                               }
                            }
                            geometry {
                              detector {
                                distance = None
                                slow_fast_beam_centre = None
                              }
                            }
                            indexing {
                              refinement_protocol {
                                n_macro_cycles = 1
                                d_min_start = 2.0
                              }
                              basis_vector_combinations.max_combinations = 10
                              stills {
                                indexer = stills
                                method_list = fft1d real_space_grid_search
                              }
                            }
                            refinement {
                              parameterisation {
                                beam.fix = *all in_spindle_plane out_spindle_plane wavelength
                                detector  {
                                  fix = all position orientation
                                  hierarchy_level=0
                                }
                                auto_reduction {
                                  action=fix
                                  min_nref_per_parameter=1
                                }
                              }
                              reflections {
                                outlier.algorithm=null
                                weighting_strategy  {
                                  override=stills
                                  delpsi_constant=1000000
                                }
                              }
                            }
                            integration {
                              integrator=stills
                              profile.fitting=False
                              background {
                                simple {
                                  outlier {
                                    algorithm = null
                                  }
                                }
                              }
                            }
                            profile {
                              gaussian_rs {
                                min_spots.overall = 0
                              }
                            }
                            '''

  target_phil = ip.parse(default_target_str)
  if write_target_file:
    with open(def_target_file, 'w') as targ:
      targ.write('{}'.format(target_phil.as_str()))

  if write_param_file:
    with open('{}/iota.param'.format(current_path), 'w') as default_param_file:
      default_param_file.write(txt_out)

  return target_phil.as_str(), txt_out

def print_params():
  """ Print out master parameters for IOTA """

  help_out = convert_phil_to_text(master_phil, att_level=1)
  txt_out = convert_phil_to_text(master_phil)

  return help_out, txt_out

# -------------------------- Backwards Compatibility ------------------------- #

class PHILFixer():
  ''' Class for backwards compatibility; will read old IOTA parameters and
  convert to current IOTA parameters (PHIL format ONLY). Can optionally
  output a new param file. Will run automatically most of the time, but can
  also be called from command-line '''

  def __init__(self):

    self.diffs = [
      ('image_triage.min_Bragg_peaks', 'image_import.miminum_Bragg_peaks'),
      ('cctbx_xfel.minimum_Bragg_peaks', 'image_import.miminum_Bragg_peaks'),
      ('cctbx_xfel.strong_sigma', 'image_import.strong_sigma'),
      ('cctbx_xfel.selection.min_sigma', 'image_import.strong_sigma'),
      ('image_conversion', 'image_import'),
      ('integrate_with', 'processing_backend'),
      ('n_processors', 'mp.n_processors'),
      ('mp_method', 'mp.method'),
      ('mp_queue', 'mp.queue'),
      ('advanced.image_viewer', 'gui.image_viewer'),
      ('advanced.monitor_mode', 'gui.monitor_mode')
    ]

    self.word_diffs = [
      ('cctbx_ha14.image_conversion.square_mode', "None", 'no_modification'),
      ('cctbx_ha14.image_triage.type', "None", 'no_triage'),
      ('cctbx_ha14.grid_search.type', "None", 'no_grid_search'),
      ('analysis.viz', "None", "no_visualization"),
      ('processing_backend', 'dials', 'cctbx.xfel'),
      ('processing_backend', 'cctbx', 'ha14')
    ]

  def read_in_phil(self, old_phil, current_phil):
    """ Extract parts of incoming PHIL not recognizable to current PHIL """
    new_phil, wrong_defs = current_phil.fetch(old_phil,
                                              track_unused_definitions=True)
    return new_phil, wrong_defs

  def make_changes(self, wrong_defs, backend='cctbx_xfel'):
    """ Find paths that have changed and change them to new paths

    :return: 'user phil' object with corrected paths
    """
    from libtbx.phil import strings_from_words
    fix_txt = ''
    fixed_defs = []

    # Change old paths to new paths
    for l in wrong_defs:
      path = l.path
      for fix in self.diffs:
        if fix[0] in path:
          new_path = path.replace(fix[0], fix[1])
          value = strings_from_words(l.object.words)
          value = self.check_values(new_path, value)
          if type(value) == list:
            value = ' '.join(value)
          entry = '{} = {}\n'.format(new_path, value)
          fix_txt += entry
          fixed_defs.append(l)

    # Change backend to the appropriate version
    remaining_defs = list(set(wrong_defs) - set(fixed_defs))
    if backend == 'ha14':
      backend_defs = ['cctbx', 'cctbx_ha14']
    else:
      backend_defs = ['dials']
    for r in remaining_defs:
      path = r.path
      for bd in backend_defs:
        if bd in path:
          new_path = path.replace(bd, 'cctbx_xfel')
          value = strings_from_words(r.object.words)
          value = self.check_values(new_path, value)
          if type(value) == list:
            value = ' '.join(value)
          entry = '{} = {}\n'.format(new_path, value)
          fix_txt += entry

    return ip.parse(fix_txt)

  def check_values(self, path, value):
    for wd in self.word_diffs:
      if wd[0] in path:
        if type(value) == list:
          for v in value:
            if wd[1] == v.replace('*', ''):
              value[value.index(v)]= v.replace(wd[1], wd[2])
              break
        else:
          value = wd[2]
    return value

  # def fix_ha14_params(self, params):
  #   # Renaming of imported images
  #   if not hasattr(params.cctbx_ha14.image_conversion, 'rename_pickle'):
  #     prefix = params.cctbx_ha14.image_conversion.rename_pickle_prefix
  #     if 'none' in str(prefix).lower():
  #       params.cctbx_ha14.image_conversion.__inject__('rename_pickle',
  #                                          'keep_file_structure')
  #       params.cctbx_ha14.image_conversion.rename_pickle_prefix = None
  #     elif 'auto' in str(prefix).lower():
  #       params.cctbx_ha14.image_conversion.__inject__('rename_pickle', 'auto_filename')
  #       params.cctbx_ha14.image_conversion.rename_pickle_prefix = None
  #     else:
  #       params.cctbx_ha14.image_conversion.__inject__('rename_pickle', 'custom_filename')
  #       if hasattr(prefix, '__iter__'):
  #         if len(prefix) > 1:
  #           params.cctbx_ha14.image_conversion.rename_pickle_prefix = \
  #             [i for i in prefix if "*" in i][0].replace('*', '')
  #         else:
  #           params.cctbx_ha14.image_conversion.rename_pickle_prefix = prefix[0]
  #
  #   if str(params.cctbx_ha14.image_conversion.square_mode).lower() == 'none':
  #     params.cctbx_ha14.image_conversion.square_mode = 'no_modification'
  #
  #   if str(params.cctbx_ha14.image_triage.type).lower() == 'none':
  #     params.cctbx_ha14.image_triage.type = 'no_triage'
  #
  #   if str(params.cctbx_ha14.grid_search.type).lower() == 'none':
  #     params.cctbx_ha14.grid_search.type = 'no_grid_search'
  #
  # def fix_old_phil(self, phil):
  #   """ Backwards compatibility: convert settings from old format to new """
  #
  #   temp_phil = master_phil.fetch(source=phil)
  #   params = temp_phil.extract()
  #
  #   if str(params.analysis.viz).lower() == 'none':
  #     params.analysis.viz = 'no_visualization'
  #
  #   fixed_phil = temp_phil.format(python_object=params)
  #
  #   return fixed_phil

  def determine_backend(self, phil=None, prm=None):
    # Try to get parameter for backend type (for backwards compatibility)
    if phil and not prm:
      prm = phil.extract()

    backend = None
    try:
      backend = prm.advanced.processing_backend
    except AttributeError:
      backend = prm.advanced.integrate_with
    finally:
      if backend and type(backend) == list:
        for b in backend:
          if b.startswith('*'):
            backend = b[1:]
            break
      return backend

  def create_current_phil(self, phil):
    # Create master PHIL based on backend and settings
    if phil.__class__.__name__ == 'scope_extract':
      backend = self.determine_backend(prm=phil)
    elif phil.__class__.__name__ == 'scope':
      backend = self.determine_backend(phil=phil)
    else:
      return None

    if not backend:
      return None

    if backend.lower() in ('ha14', 'cctbx'):
      from iota.components.iota_cctbx_ha14 import ha14_str
      new_phil = ip.parse(master_phil_str + ha14_str)
    else:
      new_phil = master_phil

    # Append GUI PHIL if it should be there:
    from libtbx.phil import find_scope
    gui_phil = find_scope(phil, 'gui')
    if gui_phil:
      from iota.components.iota_ui_base import gui_phil
      new_phil.adopt_scope(gui_phil)

    return new_phil

  def run(self, old_phil, current_phil=None, write_file=False):

    if not current_phil:
      current_phil = self.create_current_phil(phil=old_phil)

    # Get mismatched definitions
    new_phil, wrong_defs = self.read_in_phil(old_phil=old_phil,
                                             current_phil=current_phil)

    # If no mismatches found, return PHIL, else, fix it
    if not wrong_defs:
      fixed_phil = new_phil
    else:
      backend = self.determine_backend(phil=current_phil)
      phil_fixes = self.make_changes(wrong_defs, backend=backend)
      fixed_phil = new_phil.fetch(source=phil_fixes)

      if write_file:
        fixed_filepath = '{}/iota_fixed.param'.format(
          os.path.abspath(os.curdir))
        with open(fixed_filepath, 'w') as default_param_file:
          default_param_file.write(fixed_phil.as_str())

    return fixed_phil

from __future__ import absolute_import, division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 10/12/2014
Last Changed: 09/20/2019
Description : Interprets command line arguments. Initializes all IOTA starting
              parameters. Starts main log. Options for a variety of running
              modes, including resuming an aborted run.
'''

import os
import time
assert time
from contextlib import contextmanager

import dials.util.command_line as cmd

import iota.components.iota_input as inp
from iota.components import iota_utils as util
from iota.components.iota_base import ProcInfo


@contextmanager    # Will print start / stop messages around some processes
def prog_message(msg, mode='xterm'):
  if mode == 'xterm':
    cmd.Command.start(msg)
    yield
    cmd.Command.end('{} -- DONE'.format(msg))
  elif mode == 'verbose':
    print ('IOTA: {}'.format(msg))
    yield
    print ('IOTA: {} -- DONE'.format(msg))
  else:
    yield


def initialize_interface(args, phil_args=None, gui=False):
  """ Read and process input, create PHIL """

  msg = []
  input_dict = util.ginp.process_mixed_input(args.path)
  if input_dict and not gui and not input_dict['imagefiles']:
      return None, None, "IOTA_INIT_ERROR: No readable image files in path(s)!"

  # Move args that were included in paths and not processed into phil_args,
  # to try and interpret them as PHIL args
  if input_dict['badpaths']:
    phil_args.extend(input_dict['badpaths'])

  # Read in parameters, make IOTA PHIL
  iota_phil, bad_args = inp.process_input(args=args,
                                          phil_args=phil_args,
                                          paramfile=input_dict['paramfile'],
                                          gui=gui)

  # Check if any PHIL args not read into the PHIL were in fact bad paths
  if bad_args:
    input_dict['badpaths'] = [a for a in bad_args if a in input_dict['badpaths']]
    if input_dict['badpaths']:
      msg.append('Files or directories not found:')
      for badpath in input_dict['badpaths']:
        msg += '\n{}'.format(badpath)
    bad_args = [a for a in bad_args if a not in input_dict['badpaths']]
    if bad_args:
      msg += '\nThese arguments could not be interpreted: '
      for arg in bad_args:
        msg += '\n{}'.format(arg)

  return input_dict, iota_phil, msg


def initialize_new_run(phil, input_dict=None, target_phil=None):
  ''' Create base integration folder; safe phil, input, and info to file '''
  try:
    params = phil.extract()
    int_base, run_no = util.set_base_dir(dirname='integration',
                                         out_dir=params.output,
                                         get_run_no=True)
    if not os.path.isdir(int_base):
      os.makedirs(int_base)

    # Create input list file and populate param input line
    if input_dict:
      if len(input_dict['imagepaths']) >= 10:
        input_list_file = os.path.join(int_base, 'input.lst')
        with open(input_list_file, 'w') as lf:
          for f in input_dict['imagefiles']:
            lf.write('{}\n'.format(f))
          params.input = [input_list_file]
      else:
        # FIXME: This will break down if user specifies > 10 individual
        #  imagefiles from a variety of paths
        input_list_file = None
        if len(input_dict['imagefiles']) <= 10:
          params.input = input_dict['imagefiles']
        else:
          params.input = input_dict['imagepaths']
    else:
      input_list_file = None

    # Generate default backend PHIL, write to file, and update params
    target_fp = os.path.join(int_base,'target.phil')
    if target_phil:
      target_phil = inp.write_phil(phil_str=target_phil,
                                   dest_file=target_fp,
                                   write_target_file=True)
    else:
      if params.cctbx_xfel.target:
        target_phil = inp.write_phil(phil_file=params.cctbx_xfel.target,
                                     dest_file=target_fp,
                                     write_target_file=True)
      else:
        method = params.advanced.processing_backend
        target_phil, _ = inp.write_defaults(method=method,
                                            write_param_file=False,
                                            filepath=target_fp)
    params.cctbx_xfel.target = target_fp

    # Save PHIL for this run in base integration folder
    paramfile = os.path.join(int_base, 'iota_r{}.param'.format(run_no))
    phil = phil.format(python_object=params)
    with open(paramfile, 'w') as philf:
      philf.write(phil.as_str())

    # Initialize main log
    logfile = os.path.abspath(os.path.join(int_base, 'iota.log'))

    # Initialize proc.info object and save to file
    info = ProcInfo.from_args(
      iota_phil=phil.as_str(),
      target_phil=target_phil.as_str(),
      int_base=int_base,
      input_list_file=input_list_file,
      info_file=os.path.join(int_base, 'proc.info'),
      cluster_info_file=os.path.join(int_base, 'cluster.info'),
      paramfile=paramfile,
      logfile=logfile,
      run_number=run_no,
      description=params.description,
      status='initialized',
      have_results=False,
      errors=[],
      init_proc=False)
    info.export_json()
    return True, info, 'IOTA_XTERM_INIT: Initialization complete!'
  except Exception as e:
    import traceback
    traceback.print_exc()

    msg = 'IOTA_INIT_ERROR: Could not initialize run! {}'.format(e)
    return False, None, msg


def initialize_processing(paramfile, run_no):
  ''' Initialize processing for a set of images
  :param paramfile: text file with IOTA parameters
  :param run_no: number of the processing run
  :return: info: INFO object
           params: IOTA params
  '''
  try:
    phil, _ = inp.get_input_phil(paramfile=paramfile)
  except Exception as e:
    msg = 'IOTA_PROC_ERROR: Cannot import IOTA parameters! {}'.format(e)
    return False, msg
  else:
    params = phil.extract()

  # Reconstruct integration base path and get info object
  int_base = os.path.join(params.output, 'integration/{:03d}'.format(run_no))
  try:
    info_file = os.path.join(int_base, 'proc.info')
    info = ProcInfo.from_json(filepath=info_file)
  except Exception as e:
    msg = 'IOTA_PROC_ERROR: Cannot import INFO object! {}'.format(e)
    return False, msg

  # Generate input list and input base
  if not hasattr(info, 'input_list'):
    info.generate_input_list(params=params)
  common_pfx = os.path.abspath(
    os.path.dirname(os.path.commonprefix(info.input_list)))

  input_base = common_pfx
  if os.path.isdir(os.path.abspath(params.input[0])):
    new_common_pfx = os.path.commonprefix([os.path.abspath(params.input[0]), common_pfx])
    if new_common_pfx not in ('', '.'):
      input_base = new_common_pfx

  # Generate subfolder paths
  paths = dict(obj_base=os.path.join(int_base, 'image_objects'),
               fin_base=os.path.join(int_base, 'final'),
               log_base=os.path.join(int_base, 'logs'),
               viz_base=os.path.join(int_base, 'visualization'),
               tmp_base=os.path.join(int_base, 'tmp'),
               input_base=input_base)
  for bkey, bvalue in paths.items():
    if bkey == "input_base":
      continue
    if not os.path.isdir(bvalue):
      os.makedirs(bvalue)
  info.update(paths)

  # Generate filepaths for various info files
  info_files = dict(
    obj_list_file=os.path.join(info.tmp_base, 'finished_objects.lst'),
    idx_file=os.path.join(info.int_base, 'observations.pickle')
  )
  info.update(info_files)

  # Initialize stat containers
  info = generate_stat_containers(info=info, params=params)

    # Initialize main log
  util.main_log(info.logfile, '{:*^80} \n'.format(' IOTA MAIN LOG '))
  util.main_log(info.logfile, '{:-^80} \n'.format(' SETTINGS FOR THIS RUN '))
  util.main_log(info.logfile, info.iota_phil)
  util.main_log(info.logfile, '{:-^80} \n'.format('BACKEND SETTINGS'))
  util.main_log(info.logfile, info.target_phil)

  info.export_json()

  return info, params


def resume_processing(info):
  ''' Initialize run parameters for an existing run (e.g. for resuming a
      terminated run or re-submitting with new images)
  :param info: INFO object
  :return: info: Updated INFO object
           params: IOTA params
  '''

  if not info.init_proc:
    return initialize_processing(info.paramfile, info.run_number)
  else:
    try:
      phil, _ = inp.get_input_phil(paramfile=info.paramfile)
    except Exception as e:
      return None, None
    else:
      info.status = 'processing'
      return info, phil.extract()

def initialize_single_image(img, paramfile, output_file=None, output_dir=None,
                            min_bragg=10):

  phil, _ = inp.get_input_phil(paramfile=paramfile)
  params = phil.extract()

  params.input = [img]
  params.mp.n_processors = 1
  params.data_selection.image_triage.minimum_Bragg_peaks = min_bragg
  phil = phil.format(python_object=params)

  info = ProcInfo.from_args(iota_phil=phil.as_str(),
                            paramfile=paramfile)

  # Initialize output
  if output_file is not None:
    if output_dir is not None:
      output = os.path.join(os.path.abspath(output_dir), output_file)
    else:
      output = os.path.abspath(output_file)
  else:
    output = None
  info.obj_list_file = output

  info.generate_input_list(params=params)
  info = generate_stat_containers(info=info, params=params)

  return info, params

def generate_stat_containers(info, params):
  # Generate containers for processing information
  info.update(
    bookmark=0,
    merged_indices={},
    b_factors=[],
    final_objects=[],
    finished_objects=[],
    status_summary={
      'nonzero': [],
      'names': [],
      'patches': []
    },
    cluster_iterable=[],
    clusters=[],
    prime_info=[],
    user_sg='P1',
    best_pg=None,
    best_uc=None,
    msg='',
    categories=dict(
      total=(info.input_list, 'images read in', 'full_input.lst', None),
      have_diffraction=([], 'have diffraction', 'have_diffraction.lst', None),
      failed_triage=([], 'failed triage', 'failed_triage.lst', '#d73027'),
      failed_spotfinding=([], 'failed spotfinding', 'failed_spotfinding.lst', '#f46d43'),
      failed_indexing=([], 'failed indexing', 'failed_indexing.lst', '#fdae61'),
      failed_grid_search=([], 'failed grid search', 'failed_integration.lst', '#fee090'),
      failed_integration=([], 'failed integration', 'failed_integration.lst', '#fee090'),
      failed_filter=([], 'failed filter', 'failed_filter.lst', '#ffffbf'),
      integrated=([], 'integrated', 'integrated.lst', '#4575b4'),
      not_processed=(info.input_list, 'not processed', 'not_processed.lst', '#e0f3f8')),
    stats={},
    pixel_size=None,
    status='processing',
    init_proc=True,
    have_results=False
  )

  # Grid search stats dictionary (HA14 - deprecated)
  if params.advanced.processing_backend == 'ha14':
    gs_stat_keys = [('s', 'signal height', 'Signal Height'),
                    ('h', 'spot height', 'Spot Height'),
                    ('a', 'spot area', 'Spot Area')]
    info.gs_stats = {}
    for key in gs_stat_keys:
      k = key[0]
      l = key[2]
      info.stats[k] = dict(lst=[], mean=0, std=0, max=0, min=0, cons=0, label=l)

  # Statistics dictionary
  stat_keys = [('res', 'Resolution'),
               ('lres', 'Low Resolution'),
               ('strong', 'Number of spots'),
               ('mos', 'Mosaicity'),
               ('wavelength', 'X-ray Wavelength'),
               ('distance', 'Detector Distance'),
               ('beamX', 'BeamX (mm)'),
               ('beamY', 'BeamY (mm)')]
  for key in stat_keys:
    k = key[0]
    l = key[1]
    info.stats[k] = dict(lst=[], median=0, mean=0, std=0,
                         max=0, min=0, cons=0, label=l)

  return info

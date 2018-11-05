from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/18/2018
Last Changed: 11/05/2018
Description : IOTA base classes
'''

import os
import math
try:  # for Py3 compatibility
    import itertools.ifilter as filter
except ImportError:
    pass

from dxtbx import datablock as db, load as data_load
from libtbx import easy_pickle as ep
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp

from libtbx.easy_mp import parallel_map

from iota.components.iota_processing import Integrator
from threading import Thread
import iota.components.iota_utils as util

# -------------------------------- Image Base -------------------------------- #

class SingleImageBase(object):
  ''' Base class for containing all info (including data) for single image;
  can also save to file if necessary or perform conversion if necessary '''

  def __init__(self, imgpath, idx=None):
    self.img_path = imgpath
    self.obj_path = None
    self.obj_file = None
    self.int_path = None
    self.int_file = None
    self.viz_path = None
    self.int_log = None
    self.datablock = None

    if idx:
      self.img_index = idx
    else:
      self.img_index = 0


    # Processing properties
    self.status = None
    self.fail = None
    self.log_info = []

    # Final statistics (may add stuff to dictionary later, as needed)
    self.final = {'img': None,                        # Image filepath
                  'final': None,                      # Integrated pickle
                  'info': '',                         # Information (general)
                  'a': 0, 'b': 0, 'c': 0,             # Cell edges
                  'alpha': 0, 'beta': 0, 'gamma': 0,  # Cell angles
                  'sg': None,                         # Space group
                  'strong': 0,                        # Strong reflections
                  'res': 0, 'lres': 0,                # Resolution, low res.
                  'mos': 0, 'epv': 0,                 # Mosaicity, EPV
                  'wavelength': 0,                    # Wavelength
                  'distance': 0,                      # Distance
                  'beamX': 0, 'beamY': 0,             # Beam XY (mm)
                  'img_size':(0, 0),                  # Image size (pixels)
                  'pixel_size':0,                     # Pixel Size
                  'gain':0                            # Detector gain
                  }

  def write_to_file(self):
    pass


class ImageImporterBase():
  ''' Base class for image importer, which will:
        1. read an image file and extract data and info
        2. apply any modifications if requested / necessary
        3. output a datablock or file for processing
  '''

  def __init__(self, init):
    self.init = init
    self.img_type = None
    self.iparams = init.params
    self.filepath = None
    self.modify = False

  def instantiate_image_object(self, filepath, idx=None):
    ''' Override to instantiate a SingleImage object for current backend
    :param filepath: path to image file
    :return: an image object
    '''
    self.img_object = SingleImageBase(imgpath=filepath, idx=idx)

  def load_image_file(self, filepath):
    ''' Loads datablock and populates image information dictionary (can
    override to load images for old-timey HA14 processing)

    :param filepath: path to raw diffraction image (or pickle!)
    :return: datablock, error message (if any)
    '''
    # Create datablock from file
    try:
      datablock = db.DataBlockFactory.from_filenames(filenames=[filepath])[0]
    except Exception as e:
      error = 'IOTA IMPORTER ERROR: Import failed! {}'.format(e)
      print (error)
      return None, error

    # Load image infromation from datablock
    try:
      imgset = datablock.extract_imagesets()[0]
      beam = imgset.get_beam()
      s0 = beam.get_s0()
      detector = imgset.get_detector()[0]

      self.img_object.final['pixel_size'] = detector.get_pixel_size()[0]
      self.img_object.final['img_size'] = detector.get_image_size()
      self.img_object.final['beamX'] = detector.get_beam_centre(s0)[0]
      self.img_object.final['beamY'] = detector.get_beam_centre(s0)[1]
      self.img_object.final['gain'] = detector.get_gain()
      self.img_object.final['distance'] = detector.get_distance()
      self.img_object.final['wavelength'] = beam.get_wavelength()

    except Exception as e:
      error = 'IOTA IMPORTER ERROR: Information extraction failed! {}'.format(e)
      print (error)
      return datablock, error

    return datablock, None

  def modify_image(self, data=None):
    ''' Override for specific backend needs (i.e. convert, mask, and/or square
    images and output pickles for HA14'''
    return data, None

  def calculate_parameters(self, datablock=None):
    ''' Override to perform calculations of image-specific parameters
    :param datablock: Datablock with image info
    :return: datablock, error message
    '''
    return datablock, None

  def write_to_file(self, datablock, filepath):
    filename = os.path.basename(filepath)
    db_dump = db.DataBlockDumper(datablock)
    db_dump.as_file(filename=filename)

  def update_log(self, data=None, status='imported', msg=None):
    if not data:
      self.img_object.log_info.append('\n{}'.format(msg))
      self.img_object.status = 'failed import'
      self.img_object.fail = 'failed import'
    else:
      beamx = 'BEAM_X = {:<4.2f}, '.format(self.img_object.final['beamX'])
      beamy = 'BEAM_Y = {:<4.2f}, '.format(self.img_object.final['beamY'])
      pixel = 'PIXEL_SIZE = {:<8.6f}, ' \
              ''.format(self.img_object.final['pixel_size'])
      size  = 'IMG_SIZE = {:<4} X {:<4}, ' \
              ''.format(self.img_object.final['img_size'][0],
                        self.img_object.final['img_size'][1])
      dist  = 'DIST = {}'.format(self.img_object.final['distance'])
      info  = ['Parameters      :', beamx, beamy, pixel, size, dist]
      self.img_object.log_info.append(''.join(info))
      self.img_object.status = status
      self.img_object.fail = None

  def prep_output(self):
    ''' Assign output paths to image object (will be used by various modules
    to write out files in appropriate locations) '''

    # Object path (may not need)
    self.img_object.obj_path = util.make_image_path(self.img_object.img_path,
                                                    self.init.input_base,
                                                    self.init.obj_base)
    fname = util.make_filename(self.img_object.img_path, new_ext='int')
    self.img_object.obj_file = os.path.join(self.img_object.obj_path, fname)

    # Final integration pickle path
    self.img_object.int_path = util.make_image_path(self.img_object.img_path,
                                                    self.init.input_base,
                                                    self.init.fin_base)
    fname = util.make_filename(self.img_object.img_path, prefix='int',
                               new_ext='pickle')
    self.img_object.int_file = os.path.join(self.img_object.int_path, fname)

    # Processing log path for image
    self.img_object.log_path = util.make_image_path(self.img_object.img_path,
                                                    self.init.input_base,
                                                    self.init.log_base)
    fname = util.make_filename(self.img_object.img_path, new_ext='tmp')
    self.img_object.int_log = os.path.join(self.img_object.log_path, fname)

    # Visualization path (may need to deprecate)
    self.img_object.viz_path = util.make_image_path(self.img_object.img_path,
                                                    self.init.input_base,
                                                    self.init.viz_base)
    fname = util.make_filename(self.img_object.img_path, prefix='int',
                               new_ext='png')
    self.viz_file = os.path.join(self.img_object.viz_path, fname)

    # Populate the 'final' dictionary
    self.img_object.final['final'] = self.img_object.int_file
    self.img_object.final['img'] = self.img_object.img_path

  def import_image(self, input_entry):
    ''' Image importing: read file, modify if requested, make datablock
    :return: An image object with info and datablock
    '''

    # Interpret input
    if type(input_entry) == list:
      idx      = input_entry[0]
      filepath = input_entry[2]
    elif type(input_entry) == str:
      idx = None
      filepath = input_entry
    else:
      raise util.InputError('IOTA IMPORT ERROR: Unrecognized input -- {}'
                            ''.format(input_entry))

    # Instantiate image object and prep requisite paths
    self.instantiate_image_object(filepath=filepath, idx=idx)
    self.prep_output()

    # Load image (default is datablock, override for HA14-style pickling)
    self.datablock, error = self.load_image_file(filepath=filepath)

    # Log initial image information
    self.img_object.log_info.append('\n{:-^100}\n'.format(filepath))
    self.update_log(data=self.datablock, msg=error)

    # Stop there if data did not load
    if not self.datablock:
      return self.img_object, error

    # Modify image as per backend (or not)
    if self.modify:
      self.datablock, error = self.modify_image(data=self.datablock)
      if not self.datablock:
        self.update_log(data=None, msg=error)
        return self.img_object, error
      else:
        self.update_log(data=self.datablock, msg=error, status='converted')

    # Calculate image-specific parameters (work in progress)
    self.datablock, error = self.calculate_parameters(datablock=self.datablock)
    if error:
      self.update_log(data=self.datablock, msg=error)

    # Finalize and output
    self.img_object.datablock = self.datablock
    self.img_object.status = 'imported'
    return self.img_object, None

  def make_image_object(self, input_entry):
    '''Run image importer (override as needed)'''
    img_object, error = self.import_image(input_entry=input_entry)
    if error:
      print(error)
    return img_object

  def run(self, input_entry):
    return self.make_image_object(input_entry)


# ------------------------------ Processing Base ----------------------------- #

class ProcessingThreadBase(Thread):
  ''' Base class for submitting processing jobs to multiple cores '''
  def __init__(self,
               init=None,
               iterable=None,
               stage='all'):
    Thread.__init__(self, name='iota_proc')
    self.init = init
    self.iterable = iterable
    self.stage = stage
    self.importer = None
    self.integrator = None

  def import_and_process(self, input_entry):
    if self.stage != 'process':
      img_object = self.importer.run(input_entry)
    else:
      img_object = input_entry[2]
    if self.stage != 'import' and img_object.status == 'imported':
      img_object = self.integrator.run(img_object)
    return img_object

  def callback(self, result):
    """ Override for custom callback behavior (or not) """
    return result

  def create_image_iterable(self):
    # This is mostly for backward compatibility; may need to remove
    if self.stage == 'process':
      self.iterable = [[i, len(self.img_objects) + 1, j] for i, j in
                       enumerate(self.img_objects, 1)]
    else:
      if self.init.params.cctbx_ha14.selection.select_only.flag_on:
        self.iterable = [[i, len(self.init.gs_img_objects) + 1, j] for i, j in
                         enumerate(self.init.gs_img_objects, 1)]
      else:
        self.iterable = [[i, len(self.init.input_list) + 1, j] for i, j in
                         enumerate(self.init.input_list, 1)]

  def run_process(self):
    self.img_objects = parallel_map(iterable=self.iterable,
                                    func=self.import_and_process,
                                    callback=self.callback,
                                    processes=self.init.params.mp.n_processors)

  def run_analysis(self):
    """ Run analysis of integrated images """
    from iota.components.iota_analysis import Analyzer
    analysis = Analyzer(init=self.init, all_objects=self.img_objects)
    analysis.print_results()
    analysis.unit_cell_analysis()
    analysis.print_summary()
    analysis.make_prime_input()

  def process(self):
    """ Run importer and/or processor """

    # Create iterable if needed
    if self.stage != 'process':
      if not self.iterable:
        self.create_image_iterable()

    # Make sure all the important stuff has been provided
    assert self.init
    assert self.iterable
    if self.stage == 'all':
      assert (self.importer and self.integrator)
    elif self.stage == 'import':
      assert self.importer
    elif self.stage == 'process':
      assert self.integrator

    self.run_process()

  def run(self):
    self.process()

# -------------------------------- Init Base --------------------------------- #

class ProcessInfo(object):
  """ Object with all the processing info UI needs to plot results """
  def __init__(self, tmp_base, int_base):

    # Image tracking
    self.bookmark = 0
    self.finished_objects = []
    self.img_list = []
    self.unread_files = []
    self.unprocessed_images = 0
    self.int_pickles = []

    # Processing stats
    self.nref_list = []
    self.nref_xaxis = []
    self.nsref_x = None
    self.nsref_y = None
    self.nsref_median = None
    self.res_list = None
    self.res_x = None
    self.res_y = None
    self.res_median = None

    # HKL / bfactor stats
    self.user_sg = 'P1'
    self.idx_array = None
    self.merged_indices = None
    self.b_factors = []

    # Summary
    self.status_summary = {
      'nonzero': [],
      'names': [],
      'patches': []
    }

    # Clustering
    self.cluster_iterable = []
    self.cluster_info = {}

    # PRIME
    self.best_pg = None
    self.best_uc = None
    self.prime_info = {}

    # Runtime files
    self.obj_list_file         = os.path.join(tmp_base, 'finished_objects.lst')
    self.finished_pickles_file = os.path.join(int_base, 'finished_pickles.lst')
    self.cluster_info_file     = os.path.join(int_base, 'cluster_info.pickle')

    # Miscellaneous
    self.msg = ''
    self.test_attribute = 0


class InitBase(object):
  """ Base class to initialize an IOTA run """

  def __init__(self):
    """ Constructor  """
    from iota import iota_version, now
    self.ginp = util.InputFinder()
    self.pid = os.getpid()

    try:
      self.user_id = os.getlogin()
    except OSError:
      self.user_id = 'iota'

    self.iver = iota_version
    self.now = now

    self.input_base = None
    self.conv_base = None
    self.obj_base = None
    self.int_base = None

    self.params = None
    self.target_phil = None
    self.input_list = None

  def make_input_list(self):
    """ Reads input directory or directory tree and makes lists of input images.
        Optional selection of a random subset
    """

    # Read input from provided folder(s) or file(s)
    input_entries = [i for i in self.params.input if i is not None]
    input_list = self.ginp.make_input_list(input_entries, filter=True,
                                           filter_type='image')

    return input_list

  def select_image_range(self, full_list):
    """ Selects a range of images (can be complex) """
    img_range_string = str(self.params.advanced.image_range.range)
    img_range_elements = img_range_string.split(',')
    img_list = []
    for n in img_range_elements:
      if '-' in n:
        img_limits = [int(i) for i in n.split('-')]
        start = min(img_limits)
        end = max(img_limits)
        if start <= len(full_list) and end <= len(full_list):
          img_list.extend(full_list[start:end])
      else:
        if int(n) <= len(full_list):
         img_list.append(full_list[int(n)])

    if len(img_list) > 0:
      return img_list
    else:
      return full_list

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

    for i in range(random_sample_number):
      random_number = random.randrange(0, len(input_list))
      if input_list[random_number] in random_inp_list:
        while input_list[random_number] in random_inp_list:
          random_number = random.randrange(0, len(input_list))
        random_inp_list.append(input_list[random_number])
      else:
        random_inp_list.append(input_list[random_number])

    return random_inp_list

  def make_int_object_list(self):
    """ Generates list of image objects from previous grid search """
    from libtbx import easy_pickle as ep

    if self.params.cctbx_ha14.selection.select_only.grid_search_path is None:
      int_dir = util.set_base_dir('integration', True)
    else:
      int_dir = self.params.cctbx_ha14.selection.select_only.grid_search_path

    img_objects = []

    for root, dirs, files in os.walk(int_dir):
      for filename in files:
        found_file = os.path.join(root, filename)
        if found_file.endswith(('int')):
          obj = ep.load(found_file)
          img_objects.append(obj)

    return img_objects

  def initialize_interface(self):
    """ Override with interface-specific particulars """
    return True, 'BASE INTERFACE OKAY'

  def initialize_output(self):
    try:
      self.int_base = util.set_base_dir('integration',
                                        out_dir=self.params.output)
      self.obj_base = os.path.join(self.int_base, 'image_objects')
      self.fin_base = os.path.join(self.int_base, 'final')
      self.log_base = os.path.join(self.int_base, 'logs')
      self.viz_base = os.path.join(self.int_base, 'visualization')
      if not self.params.advanced.temporary_output_folder:
        self.tmp_base = os.path.join(self.int_base, 'tmp')
      else:
        self.tmp_base = os.path.join(self.params.advanced.temporary_output_folder)

      # Determine input base
      common_pfx = os.path.abspath(os.path.dirname(os.path.commonprefix(self.input_list)))
      if len(self.params.input) == 1:
        self.input_base = os.path.commonprefix([self.params.input[0], common_pfx])
      else:
        self.input_base = common_pfx

      # Generate base folders
      os.makedirs(self.int_base)
      os.makedirs(self.obj_base)
      os.makedirs(self.fin_base)
      os.makedirs(self.log_base)
      try:
        if not os.path.isdir(self.tmp_base):
          os.makedirs(self.tmp_base)
      except OSError:
        self.tmp_base = os.path.join(self.int_base, 'tmp')
        if not os.path.isdir(self.tmp_base):
          os.makedirs(self.tmp_base)

      # Designate files for init and image list iterable
      self.init_file = os.path.join(self.int_base, 'init.cfg')
      self.iter_file = os.path.join(self.int_base, 'iter.cfg')

      # Initialize info object
      self.info = ProcessInfo(int_base=self.int_base, tmp_base=self.tmp_base)
      self.info_file = os.path.join(self.int_base, 'proc.info')

      return True, 'IOTA_INIT: OUTPUT INITIALIZED'
    except Exception as e:
      return False, 'IOTA_INIT_ERROR: OUTPUT NOT INITIALIZED: {}'.format(e)


  def initialize_parameters(self):
    """ Initialize IOTA and backend parameters, write to files
    :return: True if successful, False if something failed
    """

    try:
      # If fewer images than requested processors are supplied, set the number of
      # processors to the number of images
      if self.params.mp.n_processors > len(self.input_list):
        self.params.mp.n_processors = len(self.input_list)

      # Read in backend parameters
      if self.target_phil is None:
        if self.params.advanced.processing_backend == 'ha14':
          target_file = self.params.cctbx_ha14.target
        elif self.params.advanced.processing_backend == 'cctbx.xfel':
          target_file = self.params.cctbx_xfel.target
        else:
          target_file = None

        if target_file:
          with open(target_file, 'r') as phil_file:
            self.target_phil = phil_file.readlines()
        else:
          return False, 'IOTA_INIT_ERROR: TARGET FILE NOT FOUND!'

      # Create local target file so that backend parameters stay with the run
      local_target_file = os.path.join(self.int_base, 'target.phil')
      if type(self.target_phil) == list:
        self.target_phil = '\n'.join(self.target_phil)
      with open(local_target_file, 'w') as tf:
        tf.write(self.target_phil)

      # Point IOTA parameters at local target file
      if self.params.advanced.processing_backend == 'ha14':
        self.params.cctbx_ha14.target = local_target_file
      elif self.params.advanced.processing_backend == 'cctbx.xfel':
        self.params.cctbx_xfel.target = local_target_file

      # Collect final params and convert to PHIL object
      from iota.components.iota_input import master_phil
      self.iota_phil = master_phil.format(python_object=self.params)

      return True, 'IOTA_INIT: PARAMS INITIALIZED'

    except Exception as e:
      return False, 'IOTA_INIT_ERROR: PARAMS NOT INITIALIZED, {}'.format(e)

  def initialize_main_log(self):
    """ Initialize main log (iota.log) and record starting parameters

    :return: True if successful, False if fails
    """
    try:
      # Generate text of params
      self.iota_phil_string = util.convert_phil_to_text(self.iota_phil)

      # Initialize main log
      self.logfile = os.path.abspath(os.path.join(self.int_base, 'iota.log'))

      # Log starting info
      util.main_log(self.logfile, '{:*^80} \n'.format(' IOTA MAIN LOG '))
      util.main_log(self.logfile, '{:-^80} \n'.format(' SETTINGS FOR THIS RUN '))
      util.main_log(self.logfile, self.iota_phil_string)

      # Log cctbx.xfel / DIALS settings
      util.main_log(self.logfile, '{:-^80} \n'.format('BACKEND SETTINGS'))
      util.main_log(self.logfile, self.target_phil)

      return True, 'IOTA_INIT: MAIN LOG INITIALIZED'
    except Exception as e:
      return False, 'IOTA_INIT_ERROR: MAIN LOG NOT INITIALIZED, {}'.format(e)

  def run(self):
    """ Run initialization functions (overwrite for customization)

    :return: True if successful, False if failed
    """

    # Create output file structure
    init_out, msg = self.initialize_output()
    if not init_out:
      return False, msg

    # Interface-specific options - override in subclass
    ui_init_good, msg = self.initialize_interface()
    if msg:
      print (msg)
    if not ui_init_good:
      return False, msg

    # Initialize IOTA and backend parameters
    init_param, msg = self.initialize_parameters()
    if not init_param:
      return False, msg

    # Initalize main log (iota.log)
    init_log, msg = self.initialize_main_log()
    if not init_log:
      return False, msg

    # Return True for successful initialization
    return True, msg

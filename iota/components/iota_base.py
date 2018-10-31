from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/18/2018
Last Changed: 10/30/2018
Description : IOTA base classes
'''

import os

from libtbx.easy_mp import parallel_map

from iota.components.iota_image import SingleImage
from threading import Thread
import iota.components.iota_utils as util

# ------------------------------ Processing Base ----------------------------- #

class ProcessImage():
  ''' Wrapper class to do full processing of an image '''
  def __init__(self, init, input_entry, stage='all'):
    self.init = init
    self.input_entry = input_entry
    self.stage = stage

  def run(self):
    """ Wrapper for processing function using the image object
    @param input_entry: [image_number, total_images, image_object]
    @return: image object
    """
    try:
      if type(self.input_entry[2]) is str:
        img_object = SingleImage(self.input_entry, self.init)
      else:
        img_object = self.input_entry[2]
    except Exception:
      return None

    try:
      if self.stage == 'import':
        if self.init.params.cctbx_ha14.selection.select_only.flag_on:
          img_object.import_int_file(self.init)
        else:
          img_object.import_image()
      elif self.stage == 'process':
        img_object.process()
      elif self.stage == 'all':
        img_object.import_image()
        img_object.process()
    except Exception as e:
      print ('\nDEBUG: PROCESSING ERROR!', e)
      return None

    return img_object

class ProcessGeneral(Thread):
  def __init__(self,
               init=None,
               iterable=None,
               stage='all',
               abort_file=None):
    Thread.__init__(self, name='iota_proc')
    self.init = init
    self.iterable = iterable
    self.stage = stage
    self.abort = os.path.isfile(abort_file)

  def proc_wrapper(self, input_entry):
    proc_image_instance = ProcessImage(init=self.init,
                                       input_entry=input_entry,
                                       stage=self.stage)
    return proc_image_instance.run()

  def callback(self, result):
    """ Override for custom callback behavior (or not) """
    return result

  def create_image_iterable(self):
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
                                    func=self.proc_wrapper,
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

  def run(self):
    """ Run IOTA in full-processing mode (i.e. process image from import to
    integration; allows real-time tracking of output """

    # Make sure IOTA has been initialized and a list of images provided
    assert self.init

    if not self.iterable:
      self.create_image_iterable()

    self.run_process()


# -------------------------------- Init Base --------------------------------- #

class InitGeneral(object):
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
      self.conv_base = util.set_base_dir('converted_pickles',
                                         out_dir=self.params.output)
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
        pass

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

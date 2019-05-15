from __future__ import absolute_import, division, print_function
from past.builtins import range
from six.moves import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/18/2018
Last Changed: 03/06/2019
Description : IOTA base classes
'''

import os
import json
try:  # for Py3 compatibility
    import itertools.ifilter as filter
except ImportError:
    pass

from dxtbx.model import experiment_list as exp
from libtbx.easy_mp import parallel_map
from libtbx import easy_pickle as ep

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
    self.experiments = None

    if idx:
      self.img_index = idx
    else:
      self.img_index = 0


    # Processing properties
    self.status = None
    self.fail = None
    self.log_info = []
    self.errors = []

    # Final statistics (may add stuff to dictionary later, as needed)
    self.final = {'img': None,                        # Image filepath
                  'final': None,                      # Integrated pickle
                  'observations':None,                # Miller array
                  'info': '',                         # Information (general)
                  'a': 0, 'b': 0, 'c': 0,             # Cell edges
                  'alpha': 0, 'beta': 0, 'gamma': 0,  # Cell angles
                  'sg': None,                         # Space group
                  'spots': 0,                         # Found spots
                  'indexed': 0,                       # Indexed reflections
                  'integrated': 0,                    # Integrated reflections
                  'strong': 0,                        # Strong int. reflections
                  'res': 0, 'lres': 0,                # Resolution, low res.
                  'mos': 0, 'epv': 0,                 # Mosaicity, EPV
                  'wavelength': 0,                    # Wavelength
                  'distance': 0,                      # Distance
                  'beamX': 0, 'beamY': 0,             # Beam XY (mm)
                  'img_size':(0, 0),                  # Image size (pixels)
                  'pixel_size':0,                     # Pixel Size
                  'gain':0                            # Detector gain
                  }

class ImageImporterBase():
  ''' Base class for image importer, which will:
        1. read an image file and extract data and info
        2. apply any modifications if requested / necessary
        3. output an experiment list or file for processing
  '''

  def __init__(self, init=None, info=None, write_output=True):

    self.init = init if init else None
    self.info = info if info else None
    assert (init or info)

    self.img_type = None
    self.filepath = None
    self.modify = False
    self.write_output = write_output

  def instantiate_image_object(self, filepath, idx=None):
    ''' Override to instantiate a SingleImage object for current backend
    :param filepath: path to image file
    :return: an image object
    '''
    self.img_object = SingleImageBase(imgpath=filepath, idx=idx)

  def load_image_file(self, filepath):
    ''' Loads experiment list and populates image information dictionary (can
    override to load images for old-timey HA14 processing)

    :param filepath: path to raw diffraction image (or pickle!)
    :return: experiment list, error message (if any)
    '''
    # Create experiment list from file
    try:
      experiments = exp.ExperimentListFactory.from_filenames(filenames=[filepath])
    except Exception as e:
      error = 'IOTA IMPORTER ERROR: Import failed! {}'.format(e)
      print (error)
      return None, error

    # Load image information from experiment list
    try:
      imgset = experiments.imagesets()[0]
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
      return experiments, error

    return experiments, None

  def modify_image(self, data=None):
    ''' Override for specific backend needs (i.e. convert, mask, and/or square
    images and output pickles for HA14'''
    return data, None

  def calculate_parameters(self, experiments=None):
    ''' Override to perform calculations of image-specific parameters
    :param experiments: Experiment list with image info
    :return: experiment list, error message
    '''
    return experiments, None

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

    # Get bases
    input_base = self.info.input_base if self.info else self.init.input_base
    obj_base = self.info.obj_base if self.info else self.init.obj_base
    fin_base = self.info.fin_base if self.info else self.init.fin_base
    log_base = self.info.log_base if self.info else self.init.log_base
    viz_base = self.info.viz_base if self.info else self.init.viz_base

    # Object path (may not need)
    self.img_object.obj_path = util.make_image_path(self.img_object.img_path,
                                                    input_base, obj_base)
    fname = util.make_filename(self.img_object.img_path, new_ext='int')
    self.img_object.obj_file = os.path.join(self.img_object.obj_path, fname)

    # Final integration pickle path
    self.img_object.int_path = util.make_image_path(self.img_object.img_path,
                                                    input_base, fin_base)
    fname = util.make_filename(self.img_object.img_path, prefix='int',
                               new_ext='pickle')
    self.img_object.int_file = os.path.join(self.img_object.int_path, fname)

    # Processing log path for image
    self.img_object.log_path = util.make_image_path(self.img_object.img_path,
                                                    input_base, log_base)
    fname = util.make_filename(self.img_object.img_path, new_ext='log')
    self.img_object.int_log = os.path.join(self.img_object.log_path, fname)

    # Visualization path (may need to deprecate)
    self.img_object.viz_path = util.make_image_path(self.img_object.img_path,
                                                    input_base, viz_base)
    fname = util.make_filename(self.img_object.img_path, prefix='int',
                               new_ext='png')
    self.viz_file = os.path.join(self.img_object.viz_path, fname)

    # Make paths if they don't exist already
    for path in [self.img_object.obj_path, self.img_object.int_path,
                 self.img_object.log_path, self.img_object.viz_path]:
      try:
        os.makedirs(path)
      except OSError:
        pass

    # Populate the 'final' dictionary
    self.img_object.final['final'] = self.img_object.int_file
    self.img_object.final['img'] = self.img_object.img_path

  def import_image(self, input_entry):
    ''' Image importing: read file, modify if requested, make experiment list
    :return: An image object with info and experiment list
    '''

    # Interpret input
    if type(input_entry) in (list, tuple):
      idx      = int(input_entry[0])
      filepath = str(input_entry[1])
    elif type(input_entry) == str:
      idx = None
      filepath = input_entry
    else:
      raise util.InputError('IOTA IMPORT ERROR: Unrecognized input -- {}'
                            ''.format(input_entry))

    # Instantiate image object
    self.instantiate_image_object(filepath=filepath, idx=idx)

    # Generate output paths
    if self.write_output:
      self.prep_output()

    # Load image (default is experiment list, override for HA14-style pickling)
    self.experiments, error = self.load_image_file(filepath=filepath)

    # Log initial image information
    self.img_object.log_info.append('\n{:-^100}\n'.format(filepath))
    self.update_log(data=self.experiments, msg=error)

    # Stop there if data did not load
    if not self.experiments:
      return self.img_object, error

    # Modify image as per backend (or not)
    if self.modify:
      self.experiments, error = self.modify_image(data=self.experiments)
      if not self.experiments:
        self.update_log(data=None, msg=error)
        return self.img_object, error
      else:
        self.update_log(data=self.experiments, msg=error, status='converted')

    # Calculate image-specific parameters (work in progress)
    self.experiments, error = self.calculate_parameters(
      experiments=self.experiments)
    if error:
      self.update_log(data=self.experiments, msg=error)

    # Finalize and output
    self.img_object.experiments = self.experiments
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


class ProcessingBase(Thread):
  ''' Base class for submitting processing jobs to multiple cores '''
  def __init__(self, *args, **kwargs):
    Thread.__init__(self, *args, **kwargs)

    # Set attributes from remaining kwargs
    for key, value in kwargs.items():
      setattr(self, name=key, value=value)

  def create_iterable(self, input_list):
    return [[i, len(input_list) + 1, str(j)] for i, j in enumerate(input_list, 1)]

  def import_and_process(self, input_entry):
    img_object = self.importer.run(input_entry)
    if img_object.status == 'imported':
      img_object = self.integrator.run(img_object)

    # Update main log
    if hasattr(self.info, 'logfile'):
      main_log_entry = "\n".join(img_object.log_info)
      util.main_log(self.info.logfile, main_log_entry)
      util.main_log(self.info.logfile, '\n{:-^100}\n'.format(''))

    # Set status to final
    img_object.status = 'final'
    return img_object

  def callback(self, result):
    """ Override for custom callback behavior (or not) """
    return result

  def run_process(self, iterable):
    img_objects = parallel_map(iterable=iterable,
                               func=self.import_and_process,
                               callback=self.callback,
                               processes=self.params.mp.n_processors)
    return img_objects

  def run_analysis(self):
    """ Run analysis of integrated images """
    from iota.components.iota_analysis import Analyzer
    analysis = Analyzer(info=self.info, params=self.params)
    self.info = analysis.run_all()
    self.info.export_json()

  def process(self):
    """ Run importer and/or processor """

    img_objects = self.run_process(iterable=self.info.unprocessed)
    self.info.finished_objects = img_objects
    self.run_analysis()

  def run(self):
    self.process()

  @classmethod
  def for_new_run(cls, paramfile, run_no, *args, **kwargs):
    from iota.components.iota_image import ImageImporter as Importer
    from iota.components.iota_processing import IOTAImageProcessor as Integrator

    # Initialize processing parameters
    from iota.components.iota_init import initialize_processing
    cls.info, cls.params = initialize_processing(paramfile, run_no)

    cls.importer = Importer(info=cls.info)
    cls.integrator = Integrator(iparams=cls.params)

    return cls(*args, **kwargs)

  @classmethod
  def for_existing_run(cls, info, *args, **kwargs):

    from iota.components.iota_image import ImageImporter as Importer
    from iota.components.iota_processing import IOTAImageProcessor as Integrator

    # Initialize processing parameters
    from iota.components.iota_init import resume_processing
    cls.info, cls.params = resume_processing(info)

    cls.importer = Importer(info=cls.info)
    cls.integrator = Integrator(iparams=cls.params)

    return cls(*args, **kwargs)

  @classmethod
  def for_single_image(cls, info, params, action_code='spotfinding',
                       verbose=False, *args,  **kwargs):

    from iota.components.iota_image import ImageImporter
    from iota.components.iota_processing import IOTAImageProcessor

    # Initialize processing parameters
    cls.info = info
    cls.params = params
    cls.action_code = action_code
    cls.verbose = verbose

    cls.importer = ImageImporter(info=cls.info, write_output=False)
    cls.integrator = IOTAImageProcessor(iparams=cls.params,
                                        write_pickle=False,
                                        write_logs=False,
                                        last_stage=action_code)

    return cls(*args, **kwargs)


class ProcInfo(object):
  ''' Container for all the processing info. Updated by Object Sentinel
      during processing. Stored as JSON dict and can be used to recover a
      previous run. '''

  def __init__(self, info_dict=None):
    ''' Constructor
    :param dict: dictionary of attributes, likely from JSON file
    '''

    if info_dict:
      self.update(info_dict)

  def update(self, info_dict=None, **kwargs):
    ''' Add / overwrite attributes from dictionary and/or set of kwargs '''
    for dictionary in (info_dict, kwargs):
      if dictionary:
        for key, value in dictionary.items():
          setattr(self, key, value)

  def _select_image_range(self, full_list, range_str):
    """ Selects a range of images (can be complex) """
    img_range_elements = range_str.split(',')
    img_list = []
    for n in img_range_elements:
      if '-' in n:
        img_limits = [int(i) for i in n.split('-')]
        start = min(img_limits) - 1
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

  def _select_random_subset(self, full_list, number=0):
    """ Selects random subset of input entries """
    import random

    random_inp_list = []
    if number == 0:
      if len(full_list) <= 5:
        number = len(full_list)
      elif len(full_list) <= 50:
        number = 5
      else:
        number = int(len(full_list) * 0.1)

    for i in range(number):
      random_number = random.randrange(0, len(full_list))
      if full_list[random_number] in random_inp_list:
        while full_list[random_number] in random_inp_list:
          random_number = random.randrange(0, len(full_list))
        random_inp_list.append(full_list[random_number])
      else:
        random_inp_list.append(full_list[random_number])

    return random_inp_list

  def flatten_input(self, **kw):
    ''' Generates a flat list of absolute imagefile paths
    :param kw: Typical kwargs may include a param file, path to images,
    whether to select a random subset of imagefiles and/or a range, etc.
    :return: imagefile path iterator
    '''

    if 'params' in kw:
      prm = kw['params']
    elif hasattr(self, 'params'):
      prm = self.params
    else:
      prm = None

    input_list = None
    if hasattr(self, 'input_list_file'):
      input_filepath = self.input_list_file
    elif 'filepath' in kw:
      input_filepath = kw['filepath']
    else:
      input_filepath = None

    if input_filepath:
      with open(input_filepath, 'r') as inpf:
        input_list = [i.replace('\n', '') for i in inpf.readlines()]
    else:
      if 'input' in kw:
        inputs = kw['input']
      elif prm:
        inputs = [i for i in prm.input if i is not None]
      else:
        inputs = None

      if inputs:
        input_list = util.ginp.make_input_list(inputs, filter=True,
                                                 filter_type='image')

    if prm and input_list:
      if prm.advanced.image_range.flag_on:
        input_list = self._select_image_range(input_list,
                                              prm.advanced.image_range.range)
      if prm.advanced.random_sample.flag_on:
        input_list = self._select_random_subset(input_list,
                                              prm.advanced.random_sample.number)

    return [str(i) for i in input_list]

  def generate_input_list(self, **kw):
    assert not hasattr(self, 'input_list')
    self.input_list = self.flatten_input(**kw)
    self.n_input_images = len(self.input_list)
    self.unprocessed = list(enumerate(self.input_list, 1))

  def update_input_list(self, new_input=None):
    assert hasattr(self, 'input_list') and hasattr(self, 'categories')
    if new_input:
      self.input_list = [str(i) for i in new_input]
      self.unprocessed = list(enumerate(self.input_list, self.n_input_images + 1))
      self.categories['not_processed'][0].extend(new_input)
      self.categories['total'][0].extend(new_input)
      self.n_input_images += len(new_input)

  def reset_input_list(self):
    assert hasattr(self, 'input_list') and hasattr(self, 'categories')
    self.input_list = [str(i[1]) for i in self.unprocessed]

  def get_finished_objects(self):
    if hasattr(self, 'finished_objects') and self.finished_objects:
      return (ep.load(o) for o in self.finished_objects)

  def get_finished_objects_from_filelist(self, filelist):
    assert filelist
    return (ep.load(o) for o in filelist)

  def get_finished_objects_from_file(self):
    if hasattr(self, 'obj_list_file') and os.path.isfile(self.obj_list_file):
      with open(self.obj_list_file, 'r') as objf:
        objf.seek(self.bookmark)
        obj_paths = [i.rstrip('\n') for i in objf.readlines()]
        self.bookmark = objf.tell()
      if obj_paths:
        self.finished_objects.extend(obj_paths)
        return (ep.load(o) for o in obj_paths)
      else:
        return None

  def get_final_objects(self):
    if hasattr(self, 'final_objects') and self.final_objects:
      return (ep.load(o) for o in self.final_objects)

  # def get_observations(self, to_p1=True):
  #   if hasattr(self, 'final_objects') and self.final_objects:
  #
  #     # Read final objects into a generator
  #     fin = (ep.load(o) for o in self.final_objects)
  #
  #     # Extract miller_index objects from final_objects or integrated pickles
  #     all_obs = []
  #     for obj in fin:
  #       if 'observations' in obj.final:
  #         obs = obj.final['observations'].as_non_anomalous_array()
  #       else:
  #         obs = ep.load(obj.final['final'])['observations'][0].as_non_anomalous_array()
  #       all_obs.append(obs)
  #
  #     # Combine miller_index objects into a single miller_index object
  #     with util.Capturing():
  #       observations = None
  #       for o in all_obs[1:]:
  #         if to_p1:
  #           o = o.change_symmetry('P1')
  #         if observations:
  #           observations = observations.concatenate(o, assert_is_similar_symmetry=False)
  #         else:
  #           observations = o
  #     return observations

  # def update_indices(self, filepath=None, obs=None):
  #   if not filepath:
  #     filepath = self.idx_file
  #   if not (hasattr(self, 'merged_indices')):
  #     self.merged_indices = {}
  #   if not obs:
  #     obs = ep.load(filepath)
  #
  #   with util.Capturing():
  #     p1_mrg = obs.change_symmetry('P1').merge_equivalents()
  #     p1_red = p1_mrg.redundancies()
  #   for i in p1_red:
  #     hkl = ' '.join([str(j) for j in i[0]])
  #     red = i[1]
  #     if hkl in self.merged_indices:
  #       self.merged_indices[hkl] += red
  #     else:
  #       self.merged_indices[hkl] = red
  #
  #   with open(os.path.join(self.int_base, 'merged.json'), 'w') as mjf:
  #     json.dump(self.merged_indices, mjf)

    # dict to list:
    # idx_list = [tuple([tuple(int(i) for i in k.split(' ')), v])
    #             for k, v in self.merged_indices.items()]

  def get_hkl_slice(self, sg='P1', axis='l'):

    try:
      all_obs = ep.load(self.idx_file)
    except IOError:
      return None

    with util.Capturing():
      ext_sg = str(all_obs.space_group_info()).replace(' ', '').lower()
      sg = sg.replace(' ', '').lower()
      if ext_sg != sg:
        all_obs = all_obs.change_symmetry(sg, merge_non_unique=False)

    # Recalculate redundancies and slices from updated indices
    red = all_obs.merge_equivalents().redundancies()
    return red.slice(axis=axis, slice_start=0, slice_end=0)

  def export_json(self, **kw):
    ''' Export contents as JSON dict '''

    if 'filepath' in kw:
      json_file = kw['filepath']
    elif hasattr(self, 'info_file'):
      json_file = self.info_file
    else:
      if hasattr(self, 'int_base'):
        int_base = self.int_base
      else:
        int_base = os.path.abspath(os.curdir)
      json_file = os.path.join(int_base, 'proc.info')

    try:
      with open(json_file, 'w') as jf:
        json.dump(self.__dict__, jf)
    except TypeError as e:
      raise Exception('IOTA JSON ERROR: {}'.format(e))

  @classmethod
  def from_json(cls, filepath, **kwargs):
    ''' Generate INFO object from a JSON file'''

    # Check for file
    if not os.path.isfile(filepath):
      return None

    # Read in JSON file
    with open(filepath, 'r') as json_file:
      info_dict = json.load(json_file)

    # Allow to override param values with in-code args, etc.
    info_dict.update(kwargs)

    return cls(info_dict)

  @classmethod
  def from_pickle(cls, filepath):
    ''' To recover an old pickled info object '''
    pass

  @classmethod
  def from_dict(cls, info_dict):
    return cls(info_dict)

  @classmethod
  def from_args(cls, **kwargs):
    return cls(kwargs)

  @classmethod
  def from_folder(cls, path, **kwargs):
    ''' Generate INFO object from an integration folder
    :param path: path to folder with integration results
    :param kwargs: additional keyword args
    :return: ProcInfo class generated with attributes
    '''

    # Check for folder
    if not os.path.isdir(path):
      return None

    # Check for proc.info file
    info_path = os.path.isfile(os.path.join(path, 'proc.info'))
    if info_path:
      try:
        with open(os.path.join(path, 'proc.info'), 'r') as json_file:
          info_dict = json.load(json_file)
          return cls(info_dict)
      except Exception:
        pass

    # If file not there, reconstruct from contents of folder
    # Create info dictionary
    obj_base = os.path.join(path, 'image_objects')
    fin_base = os.path.join(path, 'final')
    log_base = os.path.join(path, 'logs')
    viz_base = os.path.join(path, 'visualization')
    logfile = os.path.join(path, 'iota.log')
    info_dict = dict(int_base=path,
                     obj_base=obj_base,
                     fin_base=fin_base,
                     log_base=log_base,
                     viz_base=viz_base,
                     logfile=logfile)

    # Logfile is pretty much the most key element of an IOTA run
    if not os.path.isfile(logfile):
      return None

    # Look for the IOTA paramfile
    prm_list = util.ginp.get_file_list(path=path, ext_only='param')
    if prm_list:                             # Read from paramfile
      with open(prm_list[0], 'r') as prmf:
        info_dict['iota_phil'] = prmf.read()
    else:                                    # Last ditch: extract from log
      with open(logfile, 'r') as logf:
        lines = []
        while True:
          line = next(logf)
          if '-----' in line:
            break
        while True:
          line = next(logf)
          if '-----' in line:
            break
          else:
            lines.append(line)
      if lines:
        info_dict['iota_phil'] = ''.join(lines)

    # Look for the target PHIL file
    target_file = os.path.join(path, 'target.phil')
    if os.path.isfile(target_file):           # Read from target file
      with open(target_file, 'r') as tarf:
        info_dict['target_phil'] = tarf.read()
    else:                                     # Last ditch: extract from log
      with open(logfile, 'r') as logf:
        lines = []
        while True:
          line = next(logf)
          if 'BACKEND SETTINGS' in line:
            break
        while True:
          line = next(logf)
          if '-----' in line:
            break
          else:
            lines.append(line)
      if lines:
        info_dict['target_phil'] = ''.join(lines)

    # Generate object list
    if os.path.isdir(obj_base):
      info_dict['finished_objects'] = util.ginp.get_file_list(path, ext_only='int')

    return cls(info_dict)

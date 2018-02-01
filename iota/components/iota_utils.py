from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 12/19/2016
Last Changed: 11/03/2017
Description : Module with basic utilities of broad applications in IOTA
'''

import os
from cctbx import miller
assert miller

from libtbx import easy_pickle as ep
from libtbx import easy_run

from collections import Counter


class InputFinder():
  def __init__(self):
    self.images = ['cbf', 'img', 'corr', 'mccd', 'marccd', 'mar']
    self.datafiles = ['mtz', 'hkl']
    self.sequences = ['seq', 'fasta']
    self.pickles = ['pickle', 'int']
    self.texts = ['txt', 'lst', 'phil', 'param']
    self.pickle_type = None

  def get_file_list(self, path,
                    as_string=False,
                    ignore_ext=None,
                    ext_only=None,
                    last=None):
    ''' Runs the 'find' command to recursively get a list of filepaths. Has
    a few advangages over os.walk():
      1. Faster (by quite a lot when lots of files involved)
      2. Automatically recursive
      3. Automatically returns absolute paths
      4. Can be further modified in command line (size of file, wildcards, etc.)
    :param as_string: boolean, if true will return file list as a string, if false, as list
    :param ignore_ext:  will ignore extensions as supplied
    :param ext_only: will only find files with these extensions
    :param path: path to all data (top of folder tree)
    :param last: path to last file in a previously-generated input list (
    useful when using this to look for new files in the same folder)
    :return filepaths: list of absolute file paths
    '''
    if last is not None:
      newer_than = '-newer {}'.format(last)
    else:
      newer_than = ''
    command = 'find {} -type f {}'.format(path, newer_than)
    filepaths = easy_run.fully_buffered(command).stdout_lines

    if ignore_ext is not None:
      filepaths = [path for path in filepaths if not path.endswith(ignore_ext)]
    elif ext_only is not None:
      filepaths = [path for path in filepaths if path.endswith(ext_only)]
    filepaths = [path for path in filepaths if not
                 os.path.basename(path).startswith('.')]

    if as_string:
      return '\n'.join(filepaths)

    return filepaths

  def read_files_with_oswalk(self, path):
    '''os.walk() is typically slower than using 'find', but I am keeping it
    here just in case (and for easy switching)
    :param path:
    :return filelist: a list of absolute paths of files
    '''
    filelist = []
    for root, folder, files in os.walk(path):
      filepaths = [os.path.join(os.path.abspath(path), f) for f in files]
      filelist.extend(filepaths)

    return filelist

  def identify_file_type(self, filepath):
    ''' This will attempt to identify the filetype using several consequtive
    methods

    :param filepath: input filepath
    :return filetype: identified type of file
    '''
    filetype = self.test_extension(filepath)
    if filetype == 'unidentified':
      filetype = self.test_file(filepath)
    if filetype == 'text':
      filetype = self.test_text(filepath)

    return filetype


  def test_extension(self, filepath):

    # Check extensions
    filetype = 'unidentified'
    filename = os.path.basename(filepath)
    ext = filename.split(os.extsep)[1:]  # NOTE: may end up with multiple extensions

    for e in ext:
      e = e.lower().replace(' ', '')
      if e in self.images or e.isdigit():
        filetype = 'raw image'
      elif e in self.sequences:
        filetype = 'sequence'
      elif e in self.pickles:
        if self.pickle_type is None:
          self.pickle_type = self.test_pickle(filepath)
        filetype = self.pickle_type
      elif e == 'int':
        filetype = 'image object'
      elif e in self.texts:
        filetype = 'text'
      elif e == 'mtz':
        filetype = 'data (MTZ)'
      elif e == 'hkl':
        filetype = 'data (HKL)'
      elif e == 'pdb':
        filetype = 'coordinates'
      else:
        filetype = 'unidentified'

    return filetype

  def test_file(self, filepath):
    filename = os.path.basename(filepath)
    ext = filename.split(os.extsep)[1:]  # NOTE: may end up with multiple extensions

    # Check type using 'file' command (slow but do once per extension)
    raw_type = \
    easy_run.fully_buffered('file {}'.format(filepath)).stdout_lines
    raw_type = raw_type[0].split(':')[1]
    if '8086 relocatable' in raw_type.lower():
      self.pickle_type = self.test_pickle(filepath)
      filetype = self.pickle_type
      self.pickles.extend('.'.join(ext))
    elif 'mar area detector image' in raw_type.lower():
      filetype = 'raw image'
      self.images.append('.'.join(ext))
    elif raw_type.lower() == ' data':
      filetype = 'raw image'
      self.images.append('.'.join(ext))
    elif 'text' in raw_type.lower():
      filetype = 'text'
      self.texts.append('.'.join(ext))
    else:
      filetype = 'unidentified'

    return filetype

  def test_text(self, path):
    ''' Test the contents of a text file for it being a
          1. List of paths (i.e. an input file)
          2. A PHIL type file
    :param path: path to file
    :return: filetype determined from testing a text file
    '''
    with open(path, 'r') as tf:
      contents = tf.readlines()

    # TODO: Find a better way to catch an old-timey format
    # if '\r' in contents[0]:
    #   contents = contents[0].split('\r')

    contents = [i.replace('\n', '') for i in contents][:1000]
    content_test = [os.path.isfile(i) for i in contents]

    if Counter(content_test).most_common(1)[0][0]:
      return 'file list'
    else:
      return self.test_phil(path)

  def test_pickle(self, path):
    # Test if pickle, and if so, if it's an image or processed pickle
    pickle = ep.load(path)
    try:
      if 'DATA' in pickle:
        return 'image pickle'
      elif 'observations' in pickle:
        return 'processed pickle'
      else:
        return 'pickle'
    except TypeError:
      if hasattr(pickle, 'process'):
        return 'image object'
      else:
        return 'pickle'

  def test_phil(self, filepath):
    ''' Tests incoming PHIL file to try and determine what it's for '''

    import iotbx.phil as ip
    try:
      test_phil = ip.parse(open(filepath).read())

      # Test if IOTA parameter file
      from iota.components.iota_input import master_phil as iota_phil
      new_phil, unused = iota_phil.fetch(sources=[test_phil],
                                         track_unused_definitions=True)
      if len(unused) == 0:
        return 'IOTA settings'

      # Test if PRIME parameter file
      from prime.postrefine.mod_input import master_phil as prime_phil
      new_phil, unused = prime_phil.fetch(sources=[test_phil],
                                          track_unused_definitions=True)
      if len(unused) == 0:
        return 'PRIME settings'

      # Test if LABELIT target file
      from labelit.phil_preferences import iotbx_defs, libtbx_defs
      labelit_phil = ip.parse(input_string=iotbx_defs + libtbx_defs,
                              process_includes=True)
      new_phil, unused = labelit_phil.fetch(sources=[test_phil],
                                            track_unused_definitions=True)
      if len(unused) == 0:
        return 'LABELIT target'

      # Test if DIALS target file
      from dials.command_line.stills_process import control_phil_str, \
        dials_phil_str
      dials_phil = ip.parse(control_phil_str + dials_phil_str,
                            process_includes=True)
      new_phil, unused = dials_phil.fetch(sources=[test_phil],
                                          track_unused_definitions=True)
      if len(unused) == 0:
        return 'DIALS target'
      else:
        return 'text'
    except Exception:
      return 'text'

  def get_input(self, path, filter=True, filter_type='image', last=None):
    ''' Obtain list of files (or single file) from any input; obtain file type in input
    :param path: path to input file(s) or folder(s)
    :return: input_list: list of input file(s) (could be just one file)
             input_type: type of input file(s)
    '''

    input_list = None
    input_type = None
    suffix = 'file'

    if os.path.isfile(path):
      input_type = self.identify_file_type(path)
      if input_type == 'file list':
        with open(path, 'r') as f:
          input_list = [i.rstrip('\n') for i in f.readlines()]
          suffix = 'list'
      else:
        input_list = [os.path.abspath(path)]
    elif os.path.isdir(path):
      input_list = self.get_file_list(path, last=last)
      suffix = "folder"

    if input_list is None:
      return [], None

    if len(input_list) > 0:
      input_pairs = [(filepath, self.identify_file_type(filepath)) for
        filepath in input_list]

      if filter:
        input_pairs = [i for i in input_pairs if filter_type in i[1]]
        input_pairs = [i for i in input_pairs if not '_tmp' in i[0]]

      if len(input_pairs) > 0:
        input_list = [i[0] for i in input_pairs]
        input_types = [i[1] for i in input_pairs]
        consensus_type = Counter(input_types).most_common(1)[0][0]
        input_type = '{} {}'.format(consensus_type, suffix)
      else:
        return [], None

    return input_list, input_type

  def get_folder_type(self, path):
    if os.path.isdir(path):
      file_list = self.get_file_list(path)
      input_types = [self.identify_file_type(f) for f in file_list]
      folder_type = '{} folder'.format(Counter(input_types).most_common(1)[0][0])
      return folder_type
    else:
      return 'unknown'

  def get_file_type(self, path):
    if os.path.isfile(path):
      file_type = self.identify_file_type(path)
      if file_type == 'file list':
        with open(path, 'r') as f:
          input_list = [i.rstrip('\n') for i in f.readlines()]
          input_types = [self.identify_file_type(f) for f in input_list]
          consensus_type = Counter(input_types).most_common(1)[0][0]
          input_type = '{} list'.format(consensus_type)
      else:
        input_type = file_type
      return input_type
    else:
      return 'unknown'

  def make_input_list(self, input_entries,
                      filter=False,
                      filter_type=None,
                      last=None):
    ''' Makes input list from multiple entries
    :param input_entries: a list of input paths
    :return: input list: a list of input files
    '''

    input_list = []
    for path in input_entries:
      if path is not None:
        filepaths, _ = self.get_input(path, filter=filter,
                                      filter_type=filter_type,
                                      last=last)
        input_list.extend(filepaths)
    return input_list

class RadAverageCalculator(object):
  def __init__(self, image=None, datablock=None):
    if (image is None and datablock is None):
      print 'ERROR: Need image or datablock for Radial Average Calculator'
      return
    if datablock is None:
      from dxtbx.datablock import DataBlockFactory
      self.datablock = DataBlockFactory.from_filenames([image])[0]
    else:
      self.datablock = datablock

  def make_radial_average(self, num_bins=None, hires=None, lowres=None):
    from dials.algorithms.background import RadialAverage

    imageset = self.datablock.extract_imagesets()[0]
    beam = imageset.get_beam()
    detector = imageset.get_detector()
    scan_range = (0, len(imageset))

    summed_data = None
    summed_mask = None
    for i in range(*scan_range):
      data = imageset.get_raw_data(i)
      mask = imageset.get_mask(i)
      assert isinstance(data, tuple)
      assert isinstance(mask, tuple)
      if summed_data is None:
        summed_mask = mask
        summed_data = data
      else:
        summed_data = [ sd + d for sd, d in zip(summed_data, data) ]
        summed_mask = [ sm & m for sm, m in zip(summed_mask, mask) ]

    if num_bins is None:
      num_bins = int(sum(sum(p.get_image_size()) for p in detector) / 50)
    if lowres is not None:
      vmin = (1 / lowres) ** 2
    else:
      vmin = 0
    if hires is not None:
      vmax = (1 / hires) ** 2
    else:
      vmax = (1 / detector.get_max_resolution(beam.get_s0())) ** 2

    # Compute the radial average
    radial_average = RadialAverage(beam, detector, vmin, vmax, num_bins)
    for d, m in zip(summed_data, summed_mask):
      radial_average.add(d.as_double() / (scan_range[1] - scan_range[0]), m)
    mean = radial_average.mean()
    reso = radial_average.inv_d2()
    return mean, reso

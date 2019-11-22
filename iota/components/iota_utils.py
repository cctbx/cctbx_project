from __future__ import absolute_import, division, print_function
from six.moves import range, zip

'''
Author      : Lyubimov, A.Y.
Created     : 12/19/2016
Last Changed: 11/21/2019
Description : Module with basic utilities of broad applications in IOTA
'''

import os
import sys
from collections import Counter
import wx
from glob import glob

from cctbx import miller
assert miller
from libtbx import easy_pickle as ep, easy_run

# for Py3 compatibility
from six.moves import StringIO   # This, not io.BytesIO, works with DIALS logger

# For testing
import time
assert time

# Platform-specific stuff
# TODO: Will need to test this on Windows at some point
if wx.Platform == '__WXGTK__':
  plot_font_size = 10
  norm_font_size = 10
  button_font_size = 12
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = 'python'
elif wx.Platform == '__WXMAC__':
  plot_font_size = 9
  norm_font_size = 12
  button_font_size = 14
  LABEL_SIZE = 14
  CAPTION_SIZE = 12
  python = "Python"
elif (wx.Platform == '__WXMSW__'):
  plot_font_size = 9
  norm_font_size = 9
  button_font_size = 11
  LABEL_SIZE = 11
  CAPTION_SIZE = 9
  python = "Python"  # TODO: make sure it's right!

# --------------------------- Miscellaneous Utils ---------------------------- #

def noneset(value):
  if value == '':
    return 'None'
  elif 'none' in str(value).lower():
    return "None"
  elif value is None:
    return "None"
  else:
    return value

def makenone(value):
  if str(value).lower() in ('none', ''):
    return None
  else:
    return str(value)

class UnicodeCharacters():
  def __init__(self):
    self.alpha = u'\N{GREEK SMALL LETTER ALPHA}'.encode('utf-8')
    self.beta = u'\N{GREEK SMALL LETTER BETA}'.encode('utf-8')
    self.gamma = u'\N{GREEK SMALL LETTER GAMMA}'.encode('utf-8')
    self.sigma = u'\N{GREEK SMALL LETTER SIGMA}'.encode('utf-8')

class WxFlags():
  def __init__(self):
    self.stack = wx.TOP | wx.RIGHT | wx.LEFT
    self.expand = wx.TOP | wx.RIGHT | wx.LEFT | wx.EXPAND

class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """

  def __enter__(self):
    self._stdout = sys.stdout
    self._stderr = sys.stderr
    sys.stdout = self._io_stdout = StringIO()
    sys.stderr = self._io_stderr = StringIO()
    return self
  def __exit__(self, *args):
    bytes_str = self._io_stdout.getvalue()
    stdout_lines = bytes_str.decode('UTF-8').splitlines()
    self.extend(stdout_lines)
    sys.stdout = self._stdout

    bytes_str = self._io_stderr.getvalue()
    stderr_lines = bytes_str.decode('UTF-8').splitlines()
    self.extend(stderr_lines)
    sys.stderr = self._stderr


def convert_phil_to_text(phil, phil_file=None, att_level=0):
  """ Reads in a PHIL object and converts it to plain text; optionally writes
      out to text file if filepath is provided
  :param phil: PHIL object
  :param phil_file: absolute filepath for text file with parameters
  :return: PHIL text string
  """

  txt_out = phil.as_str(attributes_level=att_level)

  if phil_file:
    with open(phil_file, 'w') as pf:
      pf.write(txt_out)

  return txt_out

def get_mpi_rank_and_size():
  from mpi4py import MPI
  comm = MPI.COMM_WORLD
  rank = comm.Get_rank() # each process in MPI has a unique id, 0-indexed
  size = comm.Get_size() # size: number of processes running in this job
  return rank, size

def main_log(logfile, entry, print_tag=False):
  """ Write main log (so that I don't have to repeat this every time). All this
      is necessary so that I don't have to use the Python logger module, which
      creates a lot of annoying crosstalk with other cctbx.xfel modules.
  """
  if logfile is not None:
    with open(logfile, 'a') as lf:
      lf.write('{}\n'.format(entry))

  if print_tag:
    print (entry)

def set_base_dir(dirname=None, sel_flag=False, out_dir=None, get_run_no=False):
  """ Generates a base folder for converted pickles and/or integration results;
      creates subfolder numbered one more than existing
  """
  if out_dir is None and dirname is not None:
    path = os.path.abspath(os.path.join(os.curdir, dirname))
  elif out_dir is not None and dirname is None:
    path = os.path.abspath(out_dir)
  elif out_dir is None and dirname is None:
    path = os.path.abspath(os.curdir)
  else:
    path = os.path.join(os.path.abspath(out_dir), dirname)

  run_no = 1
  if os.path.isdir(path):
    dirs = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
    dirnums = [int(d) for d in dirs if d.isdigit()]
    if len(dirnums) > 0:
      n_top = max(dirnums)
      if sel_flag:
        run_no = n_top
      else:
        run_no = n_top + 1

  new_path = "{}/{:03d}".format(path, run_no)

  if get_run_no:
    return new_path, run_no
  else:
    return new_path

def find_base_dir(dirname):
  """ Function to determine the current folder name """
  def check_dirname(path, subdirname):
    if os.path.isdir(os.path.join(path, subdirname)):
      try:
        int(subdirname)
        return True
      except ValueError:
        return False
    else:
      return False

  path = os.path.abspath(os.path.join(os.curdir, dirname))
  if os.path.isdir(path):
    if len(os.listdir(path)) > 0:
      dirs = [int(i) for i in os.listdir(path) if check_dirname(path, i)]
      found_path = "{}/{:03d}".format(path, max(dirs))
    else:
      found_path = path
  else:
    found_path = os.curdir
  return found_path


def make_image_path(raw_img, input_base, base_path):
  """ Makes path for output images """
  path = os.path.dirname(raw_img)
  relpath = os.path.relpath(path, input_base)
  if relpath == '.':
    dest_folder = base_path
  else:
    dest_folder = os.path.join(base_path, relpath)
  return os.path.normpath(dest_folder)
  # return dest_folder

def make_filename(path, prefix=None, suffix=None, new_ext=None):
  bname = os.path.basename(path)
  ext = bname.split(os.extsep)[-1]
  if ext.isdigit():
    filename = bname
  else:
    fn_list = bname.split('.')
    filename = '.'.join(fn_list[0:-1])
  if prefix is not None:
    filename = "{}_{}".format(prefix, filename)
  if suffix is not None:
    filename = "{}_{}".format(filename, suffix)
  if new_ext is not None:
    filename = "{}.{}".format(filename, new_ext)
  return filename

def iota_exit(silent=False, msg=None):
  if not silent:
    from iota import iota_version, now
    if msg:
      print (msg)
    print ('\n\nIOTA version {0}'.format(iota_version))
    print ('{}\n'.format(now))
  sys.exit()

def get_size(obj, seen=None):
  """Recursively finds size of objects (by Wissam Jarjoui at SHippo,
     https://goshippo.com/blog/measure-real-size-any-python-object/)"""
  size = sys.getsizeof(obj)
  if seen is None:
    seen = set()
  obj_id = id(obj)
  if obj_id in seen:
    return 0
  seen.add(obj_id)
  if isinstance(obj, dict):
    size += sum([get_size(v, seen) for v in obj.values()])
    size += sum([get_size(k, seen) for k in obj.keys()])
  elif hasattr(obj, '__dict__'):
    size += get_size(obj.__dict__, seen)
  elif hasattr(obj, '__iter__') and not isinstance(obj,(str, bytes, bytearray)):
    size += sum([get_size(i, seen) for i in obj])

  return size


# ------------------------------ Input Finder -------------------------------- #

class InputFinder(object):
  def __init__(self):
    self.images = ['cbf', 'img', 'corr', 'mccd', 'marccd', 'mar']
    self.multi_image = ['h5', 'xtc']
    self.datafiles = ['mtz', 'hkl']
    self.sequences = ['seq', 'fasta']
    self.pickles = ['pickle', 'int']
    self.texts = ['txt', 'lst', 'phil', 'param']
    self.pickle_type = None

  def get_list_type(self, path=None, file_list=None):
    if path and os.path.isdir(path):
      file_list = self.get_input_from_folder(path=path)
      is_folder = True
    else:
      is_folder = False

    if file_list:
      input_types = [self.identify_file_type(f) for f in file_list]

      # Images can be outnumbered by other files in a folder; for that
      # reason, search for even one occurence of an image
      choices = ['raw image', 'hdf5 image', 'image pickle', 'processed pickle']
      input_type = next((s for s in input_types if s in choices), None)

      # Test for whatever's most common
      if not input_type:
        input_type = '{}'.format(Counter(input_types).most_common(1)[0][0])

      if is_folder:
        input_type += ' folder'

      return input_type
    else:
      return 'unknown'

  def get_file_type(self, path):
    if os.path.isfile(path):
      file_type = self.identify_file_type(path)
      if file_type == 'file list':
        with open(path, 'r') as f:
          input_list = [i.rstrip('\n') for i in f.readlines()]
          consensus_type = self.get_list_type(file_list=input_list)
          input_type = consensus_type.replace('folder', 'list')
      else:
        input_type = file_type
      return input_type
    else:
      return 'unknown'

  def identify_file_type(self, filepath):
    """ This will attempt to identify the filetype using several consequtive
    methods

    :param filepath: input filepath
    :return filetype: identified type of file
    """
    if isinstance(filepath, tuple) or isinstance(filepath, list):
      filepath = [i for i in filepath if os.path.exists(str(i))][0]

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
    ext = filename.split(os.extsep)[1:]  # NOTE: may get multiple extensions
    for e in ext:
      e = e.lower().replace(' ', '')
      if e in self.images or e.isdigit():
        filetype = 'raw image'
      elif e in self.multi_image:
        fmt = self.test_multi_image(path=filepath)
        filetype = '{} image'.format(fmt)
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
      elif e in ['png', 'tif', 'tiff', 'jpg', 'jpeg', 'raw', 'bmp']:
        filetype = 'picture'
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
    elif 'data' in raw_type.lower():
      try:
        img_format = self.test_multi_image(filepath)
        if img_format is not None:
          filetype = '{} image'.format(img_format)
          self.multi_image.append('.'.join(ext))
        else:
          filetype = 'raw image'
          self.images.append('.'.join(ext))
      except Exception:
        filetype = 'binary'
    elif 'text' in raw_type.lower():
      filetype = 'text'
      self.texts.append('.'.join(ext))
    else:
      filetype = 'unidentified'

    return filetype

  def test_multi_image(self, path):
    from dxtbx.format.Registry import get_format_class_dag, get_format_class_for
    dag = get_format_class_dag()

    img_format = None
    for subformat in dag.get('Format', []):
      try:
        understood = get_format_class_for(subformat).understand(path)
      except Exception as e:
        pass
      else:
        if understood:
          img_format = subformat
    if 'xtc' in img_format.lower():
      img_format = 'xtc'
    elif 'hdf5' in img_format.lower():
      img_format = 'hdf5'
    else:
      img_format = None

    return img_format

  def test_text(self, path):
    """ Test the contents of a text file for it being a
          1. List of paths (i.e. an input file)
          2. A PHIL type file
    :param path: path to file
    :return: filetype determined from testing a text file
    """
    with open(path, 'r') as tf:
      contents = tf.readlines()
    contents = [i.replace('\n', '') for i in contents][:1000]
    content_test = [os.path.isfile(i.replace(' ', '')) for item in contents
                    for i in item.split(',')]

    try:
      if (
              Counter(content_test).most_common(1)[0][0] or
              content_test.count(True) / content_test.count(False) >= 0.5
      ):
        return 'file list'
      else:
        return self.test_phil(path)
    except IndexError:
      return 'unidentified'

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
    """ Tests incoming PHIL file to try and determine what it's for """

    from iotbx import phil as ip
    from iotbx.file_reader import any_file as af

    try:
      if af(filepath).file_type == 'phil':
        test_phil = ip.parse(open(filepath).read())
      else:
        test_phil = None
    except RuntimeError:  # If not a PHIL file or a bad PHIL file
      return 'text'
    else:
      if test_phil:
        # Test if IOTA parameter file
        from iota.components.iota_input import master_phil as iota_phil
        new_phil, unused = iota_phil.fetch(sources=[test_phil],
                                           track_unused_definitions=True)
        len_test = len(test_phil.all_definitions(suppress_multiple=True))
        percent_fit = (1 - len(unused) / len_test) * 100
        if percent_fit >= 50:
          return 'IOTA settings'

        # Test if PRIME parameter file
        from prime.postrefine.mod_input import master_phil as prime_phil
        new_phil, unused = prime_phil.fetch(sources=[test_phil],
                                            track_unused_definitions=True)
        len_test = len(test_phil.all_definitions(suppress_multiple=True))
        percent_fit = (1 - len(unused) / len_test) * 100
        if percent_fit >= 50:
          return 'PRIME settings'

        # Test if LABELIT target file (LABELIT not always available)
        try:
          from labelit.phil_preferences import iotbx_defs, libtbx_defs
        except ImportError:
          pass
        else:
          labelit_phil = ip.parse(input_string=iotbx_defs + libtbx_defs,
                                  process_includes=True)
          new_phil, unused = labelit_phil.fetch(sources=[test_phil],
                                                track_unused_definitions=True)
          len_test = len(test_phil.all_definitions(suppress_multiple=True))
          percent_fit = (1 - len(unused) / len_test) * 100
          if percent_fit >= 50:
            return 'LABELIT target'

        # Test if DIALS target file
        from dials.command_line.stills_process import control_phil_str, \
          dials_phil_str
        dials_phil = ip.parse(control_phil_str + dials_phil_str,
                              process_includes=True)
        new_phil, unused = dials_phil.fetch(sources=[test_phil],
                                            track_unused_definitions=True)
        len_test = len(test_phil.all_definitions(suppress_multiple=True))
        percent_fit = (1 - len(unused) / len_test) * 100
        if percent_fit >= 50:
          return 'DIALS target'
        else:
          return 'text'
      else:
        return 'text'

  def _expand_multiple(self, path):
    import h5py

    if isinstance(path, tuple) or isinstance(path, list):
      try:
        path = [i for i in path if os.path.exists(str(i))][0]
      except IndexError:
        return []

    exp_input_list = []
    with h5py.File(path, 'r') as f:
      data_keys = f['entry']['data'].keys()
    for i in range(len(data_keys)):
      exp_input_list.append((path, i))
    return [(str(il[0]), il[1]) for il in sorted(exp_input_list)]

  def _check_for_master_hdf5(self, path):
    import h5py

    if isinstance(path, tuple) or isinstance(path, list):
      try:
        path = [i for i in path if os.path.exists(str(i))][0]
      except IndexError:
        return False

    with h5py.File(path, 'r') as f:
      if 'instrument' in f['entry'].keys():
        return True
      else:
         return False

  def _get_hdf5_entry_count(self, path):
    import h5py
    if os.path.exists(path):
      with h5py.File(path, 'r') as f:
        n_images = len(f['entry']['data'])
    else:
      n_images = 0
    return n_images

  def _pair_with_types(self, input_list,
                       filter_results=True,
                       filter_type='image'):
    # Create list of filepaths paired with file types; remove any filenames
    # with "_tmp" in them (#TODO: need a better check for incomplete images)
    input_pairs = [(fp, self.identify_file_type(fp)) for fp in input_list
                   if not '_tmp' in fp]

    # No matter what, check filter non-master hdf5 files
    checked_pairs = []
    for ip in input_pairs:
      if 'hdf5' in ip[1]:
        if self._check_for_master_hdf5(path=ip[0]):
          checked_pairs.append(ip)
      else:
        checked_pairs.append(ip)
    input_pairs = checked_pairs

    if filter_results:
      if filter_type == 'self' or filter_type is None:
        filter_type = self.get_list_type(file_list=input_list)
      input_pairs = [i for i in input_pairs if filter_type in i[1]]
    return input_pairs

  def get_input_from_list_file(self, path, double_check_type=False):
    # Double check file type if requested (only 'file list' accepted)
    if double_check_type:
      filetype = self.identify_file_type(filepath=path)
      if filetype != 'file list':
        raise InputError('IOTA INPUT ERROR: {} is not an image list file!'
                         ''.format(path))

    # Open and read file
    with open(path, 'r') as f:
      input_list = [i.rstrip('\n') for i in f.readlines()]
      if len(input_list) == 0:
        raise InputError('IOTA INPUT ERROR: {} is an empty file!'.format(path))

    # If list is comma-delimited, return as list of tuples
    if ',' in input_list[0]:
      input_list = [i.split(',') for i in input_list]
    return input_list

  def get_input_from_multi_image(self, path, double_check_type=False,
                                 expand_multiple=True):
    if double_check_type:
      filetype = self.identify_file_type(filepath=path)
      if 'hdf5' not in filetype:
        raise InputError('IOTA INPUT ERROR: {} is not an HDF5 file!'
                         ''.format(path))
    if expand_multiple:
      input_list = self._expand_multiple(path)
    else:
      input_list = [path]
    return input_list

  def get_input_from_file(self, path, expand_multiple=True):
    filetype = self.identify_file_type(filepath=path)
    if filetype == 'file list':
      input_list = self.get_input_from_list_file(path)
    elif 'hdf5' in filetype:
      input_list = self.get_input_from_multi_image(
        path=path,
        expand_multiple=expand_multiple)
    else:
      input_list = [os.path.abspath(path)]

    return input_list, filetype

  def get_input_from_folder(self, path,
                            as_string=False,
                            ignore_ext=None,
                            ext_only=None,
                            last=None,
                            min_back=None):
    """ Runs the 'find' command to recursively get a list of filepaths. Has
    a few advangages over os.walk():
      1. Faster (by quite a lot when lots of files involved)
      2. Automatically recursive
      3. Automatically returns absolute paths
      4. Can be further modified in command line (size of file, wildcards, etc.)
    :param min_back: return only files last modified this many minutes ago
    :param as_string: boolean, if true will return file list as a string,
                      if false, as list
    :param ignore_ext:  will ignore extensions as supplied
    :param ext_only: will only find files with these extensions
    :param path: path to all data (top of folder tree)
    :param last: path to last file in a previously-generated input list (
    useful when using this to look for new files in the same folder)
    :return filepaths: list (or string) with absolute file paths
    """
    if last is not None:
      newer_than = '-newer {}'.format(last)
    else:
      newer_than = ''

    if min_back is not None:
      mmin = '-mmin {}'.format(min_back)
    else:
      mmin = ''

    command = 'find {} -type f {} {}'.format(path, newer_than, mmin)
    filepaths = easy_run.fully_buffered(command).stdout_lines
    if ignore_ext is not None:
      filepaths = [path for path in filepaths if not path.endswith(ignore_ext)]
    elif ext_only is not None:
      filepaths = [path for path in filepaths if path.endswith(ext_only)]
    filepaths = [path for path in filepaths if not
                 os.path.basename(path).startswith(('.', '_'))]
    filepaths = [os.path.abspath(p) for p in filepaths]

    if as_string:
      return '\n'.join(filepaths)

    return filepaths

  def get_input(self, path,
                filter_results=True,
                filter_type='image',
                last=None,
                min_back=None,
                expand_multiple=False):
    """ Obtain list of files (or single file) from any input; obtain file type
        in input
    :param filter_results: Set to True to filter the input to include only
                         files of particular type
    :param filter_type: type of file to accept (e.g. 'image' or 'list')
    :param last: Last file to be found (search for files newer than this one)
    :param min_back: Search for files created this many minutes in the past
    :param path: path to input file(s) or folder(s)
    :return: input_list: list of input file(s) (could be just one file)
             input_type: type of input file(s)
    """

    input_list = []
    input_type = None
    suffix = 'file'
    expanded = False

    # Get input list from file
    if os.path.isfile(path):
      input_list, filetype = self.get_input_from_file(
        path=path,
        expand_multiple=expand_multiple)
      expanded = expand_multiple
      if 'list' in filetype:
        suffix = 'list'

    # Get input list from folder
    elif os.path.isdir(path):
      input_list = self.get_input_from_folder(path,
                                              last=last,
                                              min_back=min_back)
      suffix = "folder"

    if len(input_list) > 1:
      input_pairs = self._pair_with_types(input_list,
                                          filter_results=filter_results,
                                          filter_type=filter_type)
      input_list, input_types = list(zip(*input_pairs))
      input_types = list(set(input_types))
      if len(input_types) == 1:
        input_type = '{} {}'.format(input_types[0], suffix)
      else:
        if filter_type == 'image':
          input_type = 'image {}'.format(suffix)
        else:
          input_type = 'mixed input {}'.format(suffix)

      # sort input by filename and then by image index; also convert all
      # filenames to string, in case any are in unicode. Expand multiples if
      # necessary
      if expand_multiple:
        if not expanded:
          exp_input_list = []
          for item in input_pairs:
            if 'hdf5' in item[1]:
              exp_item = self._expand_multiple(path=item[0])
              if len(exp_item) > 0:
                exp_input_list.extend(exp_item)
            else:
              exp_input_list.append((item[0], 0))
          input_list = [(str(il[0]), il[1]) for il in sorted(exp_input_list)]
        else:
          input_list = [(str(il[0]), il[1]) for il in sorted(input_list)]
        input_count = len(input_list)
      else:
        input_list = [str(il) for il in sorted(input_list)]
        input_count = 0
        for ip in input_pairs:
          if 'hdf5' in ip[1].lower():
            input_count += self._get_hdf5_entry_count(path=ip[0])
          else:
            input_count += 1

    elif len(input_list) == 1:
      # This works for both a single file and a folder containing one file
      input_type = '{} {}'.format(self.identify_file_type(input_list[0]),
                                  suffix)

      # Expand multi-image file if requested
      if 'hdf5' in input_type.lower():
        if expand_multiple:
          input_list = self._expand_multiple(path=input_list[0])
          input_count = len(input_list)
        else:
          input_count = self._get_hdf5_entry_count(path=input_list[0])
      else:
        input_count = 1

    else:
      input_count = 0

    return input_list, input_type, input_count

  def make_input_list(self, paths,
                      filter_results=False,
                      filter_type=None,
                      last=None,
                      min_back=None,
                      expand_multiple=False,
                      as_tuple=False):
    """ Makes input list from multiple entries
    :param filter_results: Set to True to filter found files / folders to
                           include only those of specific type
    :param filter_type: Type of file/folder to accept
    :param last: last file/folder found
    :param min_back: how far back (in minutes) to start the search
    :param paths: one or multiple input paths
    :return: input list: a list of input files
    """

    if 'scope_extract' in type(paths).__name__:
      paths = [str(i) for i in paths if i is not None]
    elif type(paths) == str:
      paths = glob(paths)
    elif isinstance(paths, list) or isinstance(paths, tuple):
      exp_paths = []
      for p in paths:
        exp_paths.extend(glob(p))
      paths = exp_paths

    input_list = []
    total_count = 0
    for path in paths:
      if path is not None:
        filepaths, input_type, input_count = self.get_input(
          path,
          filter_results=filter_results,
          filter_type=filter_type,
          last=last,
          min_back=min_back,
          expand_multiple=expand_multiple)
        input_list.extend(filepaths)
        total_count += input_count

    if not expand_multiple and as_tuple:
      input_list = [(i, 0) for i in input_list]

    return sorted(input_list), total_count

  def process_mixed_input(self, paths):
    input_dict = dict(paramfile=None,
                      imagefiles=[], imagepaths=[],
                      objectfiles=[], objectpaths=[],
                      neither=[], badpaths=[])

    if type(paths) == str:
      raw_paths = [paths]
    elif isinstance(paths, list) or isinstance(paths, tuple):
      raw_paths = paths

    paths = []
    for path in raw_paths:
      paths.extend(glob(path))

    for path in paths:
      path = os.path.abspath(path)
      if os.path.exists(path):
        if 'IOTA settings' in self.get_file_type(path):
          input_dict['paramfile'] = path

          # If there's a paramfile, get data from it first (will always be
          # imagefiles, never objects!)
          from iota.components.iota_input import get_input_phil
          phil, _ = get_input_phil(paramfile=input_dict['paramfile'])
          prm = phil.extract()
          input_dict['imagefiles'].extend(self.make_input_list(prm.input)[0])
          input_dict['imagepaths'].extend(prm.input)
        else:
          contents, ctype, _ = self.get_input(path, filter_type='self')
          if not contents:
            continue
          if ctype is None:
            ctype = ''
          if type(contents) == str:
            contents = [contents]
          if 'object' in ctype:
            input_dict['objectfiles'].extend(contents)
            if 'folder' in ctype:
              input_dict['objectpaths'].append(os.path.abspath(path))
            else:
              input_dict['objectpaths'].append(os.path.abspath(os.path.dirname(path)))
          elif 'image' in ctype:
            input_dict['imagefiles'].extend(contents)
            if 'folder' in ctype:
              input_dict['imagepaths'].append(os.path.abspath(path))
            else:
              input_dict['imagepaths'].append(os.path.abspath(os.path.dirname(path)))
          else:
            input_dict['neither'].extend(contents)
            input_dict['badpaths'].append(path)
      else:
        input_dict['badpaths'].append(path)

      # Make paths unique
      if input_dict['imagefiles']:
        input_dict['imagefiles'] = list(set(input_dict['imagefiles']))
      if input_dict['objectfiles']:
        input_dict['objectfiles'] = list(set(input_dict['objectfiles']))
      if input_dict['imagepaths']:
        input_dict['imagepaths'] = list(set(input_dict['imagepaths']))
      if input_dict['objectpaths']:
        input_dict['objectpaths'] = list(set(input_dict['objectpaths']))

    return input_dict

ginp = InputFinder()

class ObjectFinder(object):
  """ A class for finding pickled IOTA image objects and reading in their
  contents; outputs a list of Python objects containing information about
  individual images, including a list of integrated intensities """

  def __init__(self):
    """ Constructor """
    pass

  def find_objects(self, obj_folder, read_object_files=None,
                   find_old=False, finished_only=False):
    """ Seek and import IOTA image objects

    :param finished_only:
    :param obj_folder: path to objects (which can be in subfolders)
    :param read_object_files: list of already-read-in objects
    :param find_old: find all objects in folder, regardless of other settings
    :return: list of image objects
    """
    if find_old:
      min_back = None
    else:
      min_back = -1

    # Find objects and filter out already-read objects if any
    object_files = ginp.get_input_from_folder(obj_folder,
                                              ext_only='int',
                                              min_back=min_back)

    if read_object_files is not None:
      new_object_files = list(set(object_files) - set(read_object_files))
    else:
      new_object_files = object_files

    # For backwards compatibility, read and append observations to objects
    new_objects = [self.read_object_file(i) for i in new_object_files]
    new_finished_objects = [i for i in new_objects if
                            i is not None and i.status == 'final']

    if finished_only:
      return new_finished_objects
    else:
      return new_objects

  def read_object_file(self, filepath):
    """ Load pickled image object; if necessary, extract observations from
    the image pickle associated with object, and append to object

    :param filepath: path to image object file
    :return: read-in (and modified) image object
    """
    try:
      object = ep.load(filepath)
      if object.final['final'] is not None:
        pickle_path = object.final['final']
        if os.path.isfile(pickle_path):
          pickle = ep.load(pickle_path)
          object.final['observations'] = pickle['observations'][0]
      return object
    except Exception as e:
      print ('OBJECT_IMPORT_ERROR for {}: {}'.format(filepath, e))
      return None

gobj = ObjectFinder()

# ---------------------------------- Other ----------------------------------- #

class RadAverageCalculator(object):
  def __init__(self, image=None, experiments=None):
    if (image is None and experiments is None):
      print ('ERROR: Need image or experiments for Radial Average Calculator')
      return
    if experiments is None:
      from dxtbx.model.experiment_list import ExperimentListFactory
      self.experiments = ExperimentListFactory.from_filenames([image])[0]
    else:
      self.experiments = experiments

  def make_radial_average(self, num_bins=None, hires=None, lowres=None):
    from dials.algorithms.background import RadialAverage

    imageset = self.experiments.extract_imagesets()[0]
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
        summed_data = [ sd + d for sd, d in list(zip(summed_data, data)) ]
        summed_mask = [ sm & m for sm, m in list(zip(summed_mask, mask)) ]

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
    for d, m in list(zip(summed_data, summed_mask)):
      radial_average.add(d.as_double() / (scan_range[1] - scan_range[0]), m)
    mean = radial_average.mean()
    reso = radial_average.inv_d2()
    return mean, reso

class IOTATermination(Exception):
  def __init__(self, termination):
    Exception.__init__(self, termination)

class InputError(Exception):
  def __init__(self, termination):
    Exception.__init__(self, termination)

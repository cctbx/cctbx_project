from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 12/19/2016
Last Changed: 01/06/2017
Description : Module with basic utilities of broad applications in IOTA
'''

import os
import string
from cctbx import miller
assert miller
from libtbx import easy_pickle as ep
import dxtbx
from dxtbx.format.FormatPY import FormatPY as FPy
from collections import Counter
from iota.components.iota_misc import Capturing

class GenerateInput():
  ''' Class will try to determine the type of the input file/folder  '''

  def __init__(self):
    pass

  def make_input_list(self, input_entries):
    """ Reads input directory or directory tree and makes lists of input images
        (in pickle format) using absolute path for each file. If a separate file
        with list of images is provided, parses that file and uses that as the
        input list. If random input option is selected, pulls a specified number
        of random images from the list and outputs that subset as the input list.
    """
    input_list = []

    # run through the list of multiple input entries (or just the one) and
    # concatenate the input list (right now GUI only supplies folder, but
    # this will change in future)
    for input_entry in input_entries:
      if os.path.isfile(input_entry):
        if input_entry.endswith('.lst'):  # read from file list
          with open(input_entry, 'r') as listfile:
            listfile_contents = listfile.read()
          input_list.extend(listfile_contents.splitlines())
        elif input_entry.endswith(('pickle', 'mccd', 'cbf', 'img')):
          input_list.append(input_entry)  # read in image directly

      elif os.path.isdir(input_entry):  # read from folder
        abs_inp_path = os.path.abspath(input_entry)
        for root, dirs, files in os.walk(abs_inp_path):
          for filename in files:
            input_list.append(os.path.join(root, filename))

    return input_list


  def is_text(self, filename):
    '''Determines whether the file is text or binary '''

    s = open(filename).read(512)

    text_characters = "".join(map(chr, range(32, 127)) + list("\n\r\t\b"))
    _null_trans = string.maketrans("", "")

    if not s:
      # Empty files are considered text
      return True
    if "\0" in s:
      # Files with null bytes are likely binary
      return False

    # Get the non-text characters (maps a character to itself then
    # use the 'remove' option to get rid of the text characters.)
    t = s.translate(_null_trans, text_characters)

    # If more than 30% non-text characters, then
    # this is considered a binary file
    if float(len(t)) / float(len(s)) > 0.30:
      return False

    return True

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


  def get_folder_type(self, path):
    ''' (Crudely) determine type of files inside the folder '''

    file_list = []
    for root, dirs, files in os.walk(path):
      for filename in files:
        if not filename.endswith(('log', 'tmp', 'txt')):
          found_file = os.path.join(root, filename)
          file_list.append(found_file)

    extensions = [os.path.splitext(i)[1][1:].strip().lower() for i in file_list]
    dext = Counter(extensions).most_common(1)[0][0]

    if dext == 'pickle':
      pickle = ep.load(file_list[0])
      try:
        if 'TIMESTAMP' in pickle:
          return 'converted pickles folder'
        elif 'ewald_proximal_volume' in pickle:
          return 'processed pickles folder'
        else:
          return 'raw images folder'
      except Exception:
        pass
    elif dext == 'int':
      pickle = ep.load(file_list[0])
      try:
        if hasattr(pickle, 'generate_grid'):
          return 'image objects folder'
        else:
          return 'raw images folder'
      except Exception:
        pass
    else:
      import dxtbx
      for single_file in file_list:
        try:
          with Capturing() as junk_output:
            loaded_img = dxtbx.load(single_file)
          break
        except IOError, e:
          loaded_img = None
          pass
      if loaded_img is not None:
        return 'raw images folder'
      else:
        return 'unknown folder'


  def test_file_list(self, filelist):
    # Go through list of files, test each one for type and record time
    import time

    all_start_time = time.time()
    for item in filelist:
      start_time = time.time()
      file_type = self.get_file_type(item)
      run_time = time.time() - start_time
      print 'File {0}, type "{1}"   .... {2} sec'.format(os.path.basename(item),
                                                     file_type, run_time)
    print
    print '{0} FILES INSPECTED, TOTAL TIME: {1} sec' \
          ''.format(len(filelist), time.time() - all_start_time)

  def get_file_type(self, path):

    # Test if pickle
    if FPy.understand(path):
      pickle = ep.load(path)
      try:
        if 'TIMESTAMP' in pickle:
          return 'image pickle'
        else:
          return 'pickle'
      except TypeError:
        if hasattr(pickle, 'generate_grid'):
          return 'image object'
        else:
          return 'pickle'

    # Try to load it as an image (some img files - e.g. cbf - register as
    # text files, others - e.g. mccd - register as binaries)
    else:
      try:
        with Capturing() as junk_output:
          loaded_img = dxtbx.load(path)
      except IOError, e:
        loaded_img = None
        pass

      if loaded_img is not None:
        return 'raw image'
      else:
        return 'not image'

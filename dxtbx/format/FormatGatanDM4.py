#!/usr/bin/env python
# FormatGatanDM4.py
#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#
# Experimental format for Gatan Digital Micrograph DM4 files. See
# http://www.er-c.org/cbb/info/dmformat/

from __future__ import absolute_import, division
from dxtbx.format.Format import Format
from dxtbx.format.FormatMultiImage import FormatMultiImage
import struct
import os

class FormatGatanDM4(FormatMultiImage, Format):
  '''An image reading class for images in Gatan Digital Micrograph DM4
  format, which may be multi image stacks.

  The header does not contain useful information about the geometry, therefore
  we will construct dummy objects and expect to override on import using
  site.phil.
'''

  # set up tag data types and lengths
  _tag_data_type = {
    2:('h',2),
    3:('l',4),
    4:('H',2),
    5:('L',4),
    6:('f',4),
    7:('d',8),
    8:('?',1),
    9:('c',1),
    10:('b',1),
    11:('q',8),
    15:('struct',None),
    18:('string',None),
    20:('array',None)}

  # set up image data types and lengths
  _image_data_type = {
     0:('null',None),
     1:('h',2),
     2:('f',4),
     3:('complex',8),
     4:('obsolete',None),
     5:('complex',4),
     6:('B',1),
     7:('l',4),
     8:('_rgb',4),
     9:('b',1),
    10:('H',2),
    11:('I',4),
    12:('d',8),
    13:('complex',16),
    14:('?',),
    23:('rgba',4)}

  def __init__(self, image_file, **kwargs):

    from dxtbx import IncorrectFormatError
    if not self.understand(image_file):
      raise IncorrectFormatError(self, image_file)
    FormatMultiImage.__init__(self, **kwargs)
    Format.__init__(self, image_file, **kwargs)

  @staticmethod
  def understand(image_file):

    header = FormatGatanDM4._read_header(image_file)

    if header['version'] != 4: return False
    if header['byteord'] not in [0,1]: return False
    filesize = os.stat(image_file).st_size
    if header['rootlen'] != filesize - 24: return False
    if header['sortf'] not in [0,1]: return False
    if header['closef'] not in [0,1]: return False

    return True

  @staticmethod
  def _read_header(image_file):
    '''read the minimal header metadata and structure of the root tag
    directory'''

    hd = {}

    with FormatGatanDM4.open_file(image_file, 'rb') as f:

      hd['version'] = struct.unpack('>I', f.read(4))[0]

      if hd['version'] == 3:
        hd['rootlen'] = struct.unpack('>I', f.read(4))[0]
      else:
        hd['rootlen'] = struct.unpack('>Q', f.read(8))[0]

      hd['byteord'] = struct.unpack('>I', f.read(4))[0]

      hd['sortf'] = struct.unpack('B', f.read(1))[0]

      hd['closef'] = struct.unpack('B', f.read(1))[0]

      if hd['version'] == 3:
        hd['ntags'] = struct.unpack('>I', f.read(4))[0]
      else:
        hd['ntags'] = struct.unpack('>Q', f.read(8))[0]

      hd['offset'] = f.tell()

    return hd

  def _read_tag_directory(self, f):
    '''read a tag directory from an open file. The first byte that identifies
    the tag directory is expected to have been read already'''

    tag_dir = {}

    ltname = struct.unpack('>H', f.read(2))[0]
    tag_dir['tname'] = f.read(ltname)
    tag_dir['tlen'] = struct.unpack('>Q', f.read(8))[0]
    tag_dir['sortf'] = struct.unpack('B', f.read(1))[0]
    tag_dir['closef'] = struct.unpack('B', f.read(1))[0]
    tag_dir['ntags'] = struct.unpack('>Q', f.read(8))[0]
    tag_dir['tags'] = []

    for i in range(tag_dir['ntags']):
      tag = self._read_tag(f)
      if tag is None: break # Give up if EOF reached
      tag_dir['tags'].append(tag)

    return tag_dir

  def _read_tag(self, f):
    '''read a tag from an open file, omitting tag data for groups and arrays
    (to be fulfilled later as needed)'''

    tag = {}
    tag_type = struct.unpack('B', f.read(1))[0]
    if tag_type == 0: # EOF
      return None
    elif tag_type == 20:
      return self._read_tag_directory(f)
    elif tag_type != 21:
      # unrecognized format
      raise ValueError("Unrecognised tag type while reading data")

    ltname = struct.unpack('>H', f.read(2))[0]
    tag['tname'] = f.read(ltname)
    tag['tlen'] = struct.unpack('>Q', f.read(8))[0]
    tag['tend'] = f.tell() + tag['tlen']
    assert f.read(4) == "%%%%"
    tag['ninfo'] = struct.unpack('>Q', f.read(8))[0]
    tag['info'] = struct.unpack('>' + 'Q' * tag['ninfo'], f.read(8 * tag['ninfo']))

    if tag['ninfo'] == 1: # single entry tag
      lookup = self._tag_data_type.get(tag['info'][0])

      if lookup is None: # unknown tag data type - skip this tag
        tag['value'] = None
        f.seek(tag['tend'])
        return tag

      # read the tag value
      data_type, nbytes = lookup
      tag['value'] = struct.unpack(self._byteord + data_type, f.read(nbytes))[0]

    else: # group or array, don't read now, just record the file pos
      tag['offset'] = f.tell()
      f.seek(tag['tend'])

    return tag

  def _start(self):
    '''Open the image file, read useful metadata into an internal dictionary
    and read the tags to discover where the image data are and how to read it'''

    self._header_dictionary = self._read_header(
        self._image_file)

    self._byteord = '<' if self._header_dictionary['byteord'] == 1 else '>'

    # read tags
    with FormatGatanDM4.open_file(self._image_file, 'rb') as f:

      f.seek(self._header_dictionary['offset'])

      root = []
      for i in range(self._header_dictionary['ntags']):
        root.append(self._read_tag(f))

      #self._print_tag_hierarchy(root)

    # Get indices of images to ignore because they are thumbnails
    ignore_idx = []
    for tag in root:
      if tag['tname'] == 'Thumbnails':
        vals = [v['value'] for v in self._search_tag_hierarchy('ImageIndex', [tag])]
        ignore_idx.extend(vals)

    # Get the ImageList and remove thumbnails
    imlist = []
    for tag in root:
      if tag['tname'] == 'ImageList':
        data = [t for i, t in enumerate(tag['tags']) if i not in ignore_idx]
        imlist.extend(data)

    # Should only be one item left now! Extract ImageData (leaving ImageTags,
    # Name and UniqueID aside)
    assert len(imlist) == 1
    image_data = self._search_tag_hierarchy('ImageData', imlist)
    assert len(image_data) == 1

    # data type and size
    lookup = self._search_tag_hierarchy('DataType', image_data)[0]['value']
    lookup = self._image_data_type.get(lookup)
    if lookup is None:
      raise ValueError("Unrecognised data type for image data")
    else:
      self._data_type, self._data_size = lookup
    if self._data_type not in ['h','f','B','l','b','H','I','d','?']:
      raise ValueError("Unrecognised data type for image data")

    # image dimensions
    dimensions = self._search_tag_hierarchy('Dimensions', image_data)[0]['tags']
    self._image_size = dimensions[0]['value'], dimensions[1]['value']

    # work out how to read the data
    data = self._search_tag_hierarchy('Data', image_data)
    assert len(data) == 1
    data = data[0]
    self._image_num_elements = self._image_size[0] * self._image_size[1]
    self._total_num_elements = data['info'][2]
    self._num_images = self._total_num_elements // self._image_num_elements
    assert self._total_num_elements % self._image_num_elements == 0
    self._data_offset = data['offset']

    return

  def _print_tag_hierarchy(self, tag_list, prefix=""):
    '''for debugging'''

    for tag in tag_list:
      name = tag['tname']
      if name == '': name = 'NONAME'
      if 'tags' in tag: # tag directory
        print prefix + "DIR: {0}".format(name)
        prefix += '  '
        self._print_tag_hierarchy(tag['tags'], prefix)
        prefix = prefix[:-2]
      else: # standard tag
        val = tag.get('value')
        if val is not None:
          val = ': ' + str(val)
        else:
          val = ': ' + str(tag['info'])
        print prefix + name + val

  def _search_tag_hierarchy(self, name, tag_list):
    '''Search a list of tags for any with the given name.'''
    result = []
    for tag in tag_list:
      if tag['tname'] == name:
        result.extend([tag])
      else:
        if 'tags' in tag:
          result.extend(self._search_tag_hierarchy(name, tag['tags']))
    return result

  def get_num_images(self):
    return self._num_images

  def get_goniometer(self, index=None):
    return Format.get_goniometer(self)

  def get_detector(self, index=None):
    return Format.get_detector(self)

  def get_beam(self, index=None):
    return Format.get_beam(self)

  def get_scan(self, index=None):
    if index == None:
      return Format.get_scan(self)
    else:
      scan = Format.get_scan(self)
      return scan[index]

  def get_image_file(self, index=None):
    return Format.get_image_file(self)

  def _goniometer(self):
    '''Dummy goniometer, 'vertical' as the images are viewed'''

    return self._goniometer_factory.known_axis((0,-1,0))

  def _detector(self):
    '''Dummy detector. It may be possible to extract some of the required
    details from the image metadata, but not sure how consistent they are.'''

    # Assuming Gatan OneView camera, with 15um pixels, but trusted range is
    # (probably) wrong.
    pixel_size = 0.015, 0.015
    image_size = self._image_size
    trusted_range = (-1, 65535)
    distance=2000
    beam_centre = [(p * i) / 2 for p, i in zip(pixel_size, image_size)]
    d = self._detector_factory.simple('PAD', distance, beam_centre, '+x', '-y',
                              pixel_size, image_size, trusted_range)
    # Not sure what the gain is. A Gatan OneView has a scintillator and fibre
    # optic, so is not a direct detection camera and probably has gain != 1
    #for p in d: p.set_gain(8)
    return d

  def _beam(self):
    '''Dummy unpolarized beam, energy 200 keV'''

    wavelength = 0.02508
    return self._beam_factory.make_polarized_beam(
        sample_to_source=(0.0, 0.0, 1.0),
        wavelength=wavelength,
        polarization=(0, 1, 0),
        polarization_fraction=0.5)

  def _scan(self):
    '''Dummy scan for this stack'''

    nframes = self.get_num_images()
    image_range = (1, nframes)
    exposure_times = 0.0
    oscillation = (0, 0.5)
    epochs = [0] * nframes

    return self._scan_factory.make_scan(image_range, exposure_times,
      oscillation, epochs, deg=True)

  def get_raw_data(self, index):

    from boost.python import streambuf
    from scitbx.array_family import flex

    # is this image a type we can read?
    assert self._data_type in ['h', 'f', 'B', 'l', 'b', 'H', 'I', 'd']

    with FormatGatanDM4.open_file(self._image_file, 'rb') as f:
      f.seek(self._data_offset)

      skip_bytes = index * self._image_num_elements * self._data_size
      f.seek(skip_bytes, whence=1)

      if self._data_type == 'f':
        from dxtbx import read_float32
        raw_data = read_float32(streambuf(f), self._image_num_elements)
      elif self._data_type == 'B':
        from dxtbx import read_uint8
        raw_data = read_uint8(streambuf(f), self._image_num_elements)
      elif self._data_type == 'h':
        from dxtbx import read_int16
        raw_data = read_int16(streambuf(f), self._image_num_elements)
      elif self._data_type == 'H':
        from dxtbx import read_uint16
        raw_data = read_uint16(streambuf(f), self._image_num_elements)
      elif self._data_type == 'l':
        from dxtbx import read_int32
        raw_data = read_int32(streambuf(f), self._image_num_elements)
      elif self._data_type == 'I':
        from dxtbx import read_uint32
        raw_data = read_uint32(streambuf(f), self._image_num_elements)

      # no C++ reader for remaining types (should be unusual anyway)
      else:
        vals = struct.unpack(self._byteord + self._data_type * self._image_num_elements,
          f.read(self._data_size * self._image_num_elements))
        if self._data_type == 'd':
          raw_data = flex.double(vals)
        if self._data_type in ['b', '?']:
          raw_data = flex.int(vals)

    raw_data.reshape(flex.grid(self._image_size[1], self._image_size[0]))
    return raw_data


if __name__ == '__main__':
  import sys
  for arg in sys.argv[1:]:
    print FormatGatanDM4.understand(arg)
